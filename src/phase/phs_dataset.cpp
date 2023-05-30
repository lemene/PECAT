
#include "phs_dataset.hpp"

#include "../utils/logger.hpp"

#include "phs_options.hpp"


namespace fsa {

void PhsDataset::Load() {


    LOG(INFO)("Preload reads and contigs");
    PreloadReads();

    if (!opts_.vcf_fname_.empty()) {
        LOG(INFO)("Load VCF");
        LoadSnpFromVcf(opts_.vcf_fname_, rd_store_.GetStringPool());
        opts_.phase_opts_.using_vcf = true;
        LOG(INFO)("Set using vcf");
    }

    LOG(INFO)("Load overlaps");
    LoadOverlaps(opts_.ol_fname_);

    if (!opts_.rd2rd.empty()) {
        LoadAva(opts_.rd2rd);
    }

    LOG(INFO)("Load reads and contigs");
    LoadReads();
    
}   


void PhsDataset::LoadOverlaps(const std::string &fname) {
    auto filter_simple = [&](Overlap& o) {
        if (!opts_.filter_opts_.Valid(o)) {
            return false;
        }

        if (contig_ids_.find(o.b_.id) == contig_ids_.end()) {
            return false;
        }

//        // TODO for debug
//        if (o.b_.start < 66229153 - 1000000 || o.b_.start > 66332061 + 1000000) return false;

        return true;
    };

    ol_store_.LoadFast(fname, "", (size_t)opts_.thread_size_, filter_simple);

    if (ol_store_.Size() > 0) {
        LOG(INFO)("Overlap size: %zd", ol_store_.Size());
    } else {
        LOG(ERROR)("No overlap was loaded");
    }

    ol_store_.Save(opts_.OutputPath("load.paf"));
}


void PhsDataset::LoadAva(const std::string &fname) {
 

    ol_ava_.LoadFast(fname, "", (size_t)opts_.thread_size_);

    if (ol_ava_.Size() > 0) {
        LOG(INFO)("Overlap size: %zd", ol_ava_.Size());
    } else {
        LOG(ERROR)("No overlap was loaded");
    }
    ol_ava_.Group(ava_groups_, opts_.thread_size_);
}

void PhsDataset::LoadSnpFromVcf(const std::string &fname, StringPool& string_pool) {
    GzFileReader vcf(fname);

    size_t count = 0;
    auto line = vcf.GetNoEmptyLine();
    for (auto line = vcf.GetNoEmptyLine(); !line.empty(); line = vcf.GetNoEmptyLine()) {
        if (line[0] == '#') continue;
        const int CHROM = 0;
        const int POS = 1;
        const int REF = 3;
        const int ALT = 4;
        const int FILTER = 6;
        const int SAMPLE = 9;
        auto its = SplitStringBySpace(line);
        assert(its.size() >= 10);
        std::array<uint8_t,2> bases = {(uint8_t)-1, (uint8_t)-1};
        if (its[FILTER] == "PASS" && its[SAMPLE].size() >= 3) {

            if (its[SAMPLE].find("0/1") == 0 ) {
                if (its[REF].size() == 1 && its[ALT].size() == 1) {
                    bases[0] = DnaSeq::Serial(its[REF][0]);
                    bases[1] = DnaSeq::Serial(its[ALT][0]);
                    count ++;
                }

            } else if (its[SAMPLE].find("1/2") == 0) {
                if (its[ALT].size() == 3) { //A,C
                    bases[0] = DnaSeq::Serial(its[ALT][0]);
                    bases[1] = DnaSeq::Serial(its[ALT][2]);
                    count ++;
                }
            } else {
                // pass
            }
        }

        if (bases[0] != (uint8_t)-1) {
            auto id = string_pool.QueryIdByString(its[CHROM]);
            auto pos = std::stoi(its[POS]) - 1;
            if (bases[0] > bases[1]) std::swap(bases[0], bases[1]);

            snps_[id][pos] = bases;
        }

        
    }
    LOG(INFO)("Load SNPs: count = %zd", count);
}

void PhsDataset::LoadSnpFromVariants(const std::string &fname, StringPool& string_pool) {
    GzFileReader vcf(fname);

    size_t count = 0;
    auto line = vcf.GetNoEmptyLine();
    for (auto line = vcf.GetNoEmptyLine(); !line.empty(); line = vcf.GetNoEmptyLine()) {
        auto its = SplitStringBySpace(line);
        
        if (its.size() == 13) {
            //if (std::stoi(its[12]) == 1) {
            if (its[12] == "1") {

                int counts[4] = {0, 0, 0, 0};
                counts[0] = std::stoi(its[2]);
                counts[1] = std::stoi(its[3]);
                counts[2] = std::stoi(its[4]);
                counts[3] = std::stoi(its[5]);

                std::array<int8_t, 4> argmax = {0, 1, 2, 3};
                std::sort(argmax.begin(), argmax.end(), [&counts](int8_t a, int8_t b) {
                    return counts[a] > counts[b];
                });

                auto id = string_pool.QueryIdByString(its[0]);
                auto pos = std::stoi(its[1]);
                if (argmax[0] < argmax[1]) {
                    snps_[id][pos][0] = argmax[0];
                    snps_[id][pos][1] = argmax[1];
                } else {
                    snps_[id][pos][0] = argmax[1];
                    snps_[id][pos][1] = argmax[0];
                }
                assert(snps_[id][pos][0] <4 && snps_[id][pos][1] <4 && snps_[id][pos][0] < snps_[id][pos][1]);
                assert(snps_[id][pos][0] != snps_[id][pos][1]);
                count++;
                
            } 
        } else {
            LOG(WARNING)("Load variants: wrong number of columns");
        }
        
    }
    LOG(INFO)("Load SNPs: count = %zd", count);
}

std::unordered_set<int> PhsDataset::QueryGroup(int id) {
    auto iter = ava_groups_.find(id);
    if (iter != ava_groups_.end()) {
        std::unordered_set<int> sss;
        for (auto ii : iter->second) {
            sss.insert(ii.first);
        }
        return sss;
    } else {
        return {};
    }
}

void PhsDataset::LoadReads() {
    std::unordered_set<Seq::Id> ids;

    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        ids.insert(ol_store_.Get(i).a_.id);
        read_ids2_[ol_store_.Get(i).a_.id]++;
    }

    rd_store_.Load(opts_.rd_fname_, "", true, ids);
    rd_store_.Load(opts_.ctg_fname_, "", true, contig_ids_);
}

void PhsDataset::PreloadReads() {
    rd_store_.Load(opts_.rd_fname_, "", false, std::unordered_set<Seq::Id>());
    std::array<size_t, 2> ctg_id_range { rd_store_.GetIdRange()[1], 0};
    rd_store_.Load(opts_.ctg_fname_, "", false, std::unordered_set<Seq::Id>());
    ctg_id_range[1] = rd_store_.GetIdRange()[1];
    
    if (!opts_.ctgname_fname_.empty()) {
        GzFileReader reader(opts_.ctgname_fname_);
        if (reader.Valid()) {
            std::string name = reader.GetNoEmptyLine();
            while (!name.empty()) {
                contig_names_.insert(name);
                contig_ids_.insert(rd_store_.QueryIdByName(name));
                name = reader.GetNoEmptyLine();
            }
        } else {
            LOG(ERROR)("Failed to load file: %s", opts_.ctgname_fname_.c_str());
        }
    } else {
        for (size_t i = ctg_id_range[0]; i < ctg_id_range[1]; ++i) {
            contig_names_.insert(rd_store_.QueryNameById(i));
            contig_ids_.insert(i);
        }
    }

    if (opts_.loglevel >= 4) rd_store_.SaveIdToName(opts_.OutputPath("id2name.gz"));

}

std::vector<Seq::Id> PhsDataset::GetSortedContigs() const {
    std::vector<Seq::Id> sorted(contig_ids_.begin(), contig_ids_.end());

    std::unordered_map<Seq::Id, size_t> counts;

    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);
        counts[o.b_.id] += 1;
    }

    for (auto &c : counts) {
        c.second = c.second * c.second  / std::max<size_t>(1, (rd_store_.GetSeqLength(c.first) / 1000000)); 
    }

    std::sort(sorted.begin(), sorted.end(), [&](Seq::Id a, Seq::Id b) {
        return counts[a] > counts[b];
    });

    return sorted;
}

}