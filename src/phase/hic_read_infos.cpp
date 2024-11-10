#include "hic_read_infos.hpp"

#include "../read_store.hpp"
#include "../overlap_store.hpp"
#include "../utils/project_file.hpp"

namespace fsa {


void HicReadInfos::Build(const std::string& fn_hic1, const std::string& fn_paf1, 
                         const std::string& fn_hic2, const std::string& fn_paf2, 
                         const std::string &fn_vars, size_t thread_size) {
    SnpStore vars(string_pool_);
    //vars.LoadFromVcf(fn_vars);
    vars.Load(fn_vars);
    

    size_t sp_offset = string_pool_.Size();
    LOG(INFO)("Scan HiC1");
    BuildOne(fn_hic1, fn_paf1, 0, sp_offset, vars, thread_size);
    LOG(INFO)("Scan HiC2");
    BuildOne(fn_hic2, fn_paf2, 1, sp_offset, vars, thread_size);
}

void HicReadInfos::BuildOne(const std::string& fn_hic, const std::string& fn_paf, size_t ihic, size_t sp_offset, const SnpStore &vars, size_t thread_size) {
    assert(ihic == 0 || ihic == 1);
    std::mutex mutex;
    ReadStore hic_store(string_pool_, sp_offset);
    hic_store.Load(fn_hic);

    LOG(INFO)("Id range: %zd, %zd",hic_store.GetIdRange()[0], hic_store.GetIdRange()[1]);
    
    size_t hiclen = hic_store.GetSeqLength(hic_store.GetIdRange()[0]);
    LOG(INFO)("hiclen: %zd",hiclen);
    
    
    auto filter = [](Overlap &o) {
        return o.Identity()*o.AlignedLength() > 0.95*o.a_.len;
    };
    
    OverlapStore ol_store(string_pool_);
    ol_store.Load(fn_paf, "", thread_size, filter);
    LOG(INFO)("ols_store.Load: %zd",ol_store.Size());

    auto combine_func = [this, &mutex, ihic](const std::vector<std::pair<Seq::Id, HicReadMapping>>& hic_maps) {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto &&i : hic_maps) {
            infos_[i.first].hic[ihic].push_back(std::move(i.second));
        }
    };

    std::atomic<size_t> index { 0 };
    auto work_func = [&](int tid) {

        std::vector<std::pair<Seq::Id, HicReadMapping>> hic_maps;

        for (size_t i = index.fetch_add(1); i < ol_store.Size(); i = index.fetch_add(1)) {
            const auto &ol = ol_store.Get(i);

            // filter low-quality mapping
            if (ol.AlignedLength() < hiclen*0.9) continue;
            hic_maps.push_back(std::make_pair(ol.a_.id, GetMappingInfo(ol, hic_store, vars)));

            if (hic_maps.size() > 10000) {
                combine_func(hic_maps);
                hic_maps.clear();
            }
            //printf("XXX(%zd): tid=%d, hic_maps.size()=%zd\n", i, tid, hic_maps.size());
        }
        if (hic_maps.size() > 0) {
            combine_func(hic_maps);
            hic_maps.clear();
        }
    };
    
    MultiThreadRun(thread_size, work_func);
}


HicReadMapping HicReadInfos::GetMappingInfo(const Overlap& ol, const ReadStore& hic_store, const class SnpStore &vars) {
    HicReadMapping info;
    info.ctg = ol.b_.id;
    info.offset = ol.b_.start;

    size_t ctg_off = 0;
    size_t rd_off = 0;
    const auto &rd = hic_store.GetSeq(ol.a_.id); 
    const auto& ctgvar = vars.Get(ol.b_.id);
    const int C = 3;

    if (!ctgvar.empty()) {
        
        for (const auto &d : ol.detail_) {
            switch (d.type) {
            case 'M':
                if (d.len >= C) {
                    for (int i=C/2; i<d.len-C/2; ++i) {
                        size_t ctg_i = ol.b_.strand == 0 ? ol.b_.start+ctg_off+i : ol.b_.end-ctg_off-i-1;
                        size_t rd_i = ol.a_.strand == 0 ? ol.a_.start+rd_off+i : ol.a_.end-rd_off-i-1;

                        auto snp = ctgvar.find(ctg_i);
                        if (snp != ctgvar.end()) {
                            uint8_t rd_b = ol.SameDirect() ? rd[rd_i] : 3 - rd[rd_i];
                            if (rd_b == snp->second[0] || rd_b == snp->second[1]) {
                                info.alleles.push_back({{(uint32_t)ol.b_.id, (uint32_t)ctg_i}, rd_b});
                            }
                        }
                    }
                }
                
                ctg_off += d.len;
                rd_off += d.len;
                break;

            case 'D':
                ctg_off += d.len;
                break;

            case 'I':
                rd_off += d.len;
                break;

            case '=':
            default:
                LOG(ERROR)("Not support cigar type '%c'.", d.type);
            }
        }
    
    }
    return info;
}

void HicReadInfos::Save(const std::string& fname) const {
    GzFileWriter writer(fname);

    auto dump_HicReadMapping = [this](std::ostringstream &oss, const std::vector<HicReadMapping>& maps) {
        for (const auto &m : maps) {
            oss << string_pool_.QueryStringById(m.ctg) << "_" << m.offset;
            for (const auto &s : m.alleles) {
                oss << "_" << s.site.offset << "-" << (size_t)(s.base);
            }
            if (&m != &maps.back()) {
                oss << ",";
            }
        }
    };

    if (writer.Valid()) {
        std::ostringstream oss;
        for (auto &i : infos_) {
            if (i.second.hic[0].size() > 0 && i.second.hic[1].size() > 0) {
                oss << string_pool_.QueryStringById(i.first) << " ";
                dump_HicReadMapping(oss, i.second.hic[0]);
                oss << "|";
                dump_HicReadMapping(oss, i.second.hic[1]);
                oss << "\n";
                writer.Flush(oss);
                oss.str("");
            }
        }
    } else {
        LOG(WARNING)("Failed to write the file: %s", fname.c_str());
    }
}

void HicReadInfos::Load(const std::string& fname) {
    GzFileReader reader(fname);

    if (reader.Valid()) {
        std::string line = reader.GetNoEmptyLine();
        while (!line.empty()) {
            LoadLine(line);
            line = reader.GetNoEmptyLine();
        }
        if (!reader.IsEnd()) {
            LOG(WARNING)("Not all data is loaded: %s", fname.c_str());
        }
    } else {
        LOG(WARNING)("Failed to read the file: %s", fname.c_str());
    }
}

void HicReadInfos::LoadLine(const std::string &line) {
    auto load_mapping = [this](HicReadMapping& m, const std::string& str) {
        auto its0 = SplitStringByChar(str, '_');
        if (its0.size() >= 2) {
            m.ctg = string_pool_.GetIdByString(its0[0]);
            m.offset = std::stoul(its0[1]);

            for (size_t i = 2; i < its0.size(); ++i) {
                m.alleles.push_back(SnpAllele());

                auto its1 = SplitStringByChar(its0[i], '-');
                m.alleles.back().site.ctg = m.ctg;
                m.alleles.back().site.offset = std::stoul(its1[0]);
                m.alleles.back().base = std::stoul(its1[1]);
            }
        }
    };

    auto its0 = SplitStringBySpace(line);
    if (its0.size() == 2) {
        auto rid = string_pool_.GetIdByString(its0[0]);
        auto& info = infos_[rid].hic;

        auto its1 = SplitStringByChar(its0[1], '|');
        for (size_t i = 0; i < 2; ++i) {
            auto its2 = SplitStringByChar(its1[i], ',');
            for (size_t ii = 0; ii < its2.size(); ++ii) {
                info[i].push_back(HicReadMapping());
                load_mapping(info[i].back(), its2[ii]);
            }
        }
    }

}

}  // namespace fsa
