#include "pol_dataset.hpp"

#include <random>

#include "pol_options.hpp"
#include "../correct/align/tool_aligner.hpp"

namespace fsa {

void PolDataset::Load() {

    LoadReadIds();

    LoadOverlaps();
    
    // load mapping
    if (!opts_.rd_2_ref_fname_.empty()) {
        LoadMappings();
    }

    LoadReads();

    if (opts_.debug) {
        read_store_.SaveIdToName("id_2_name");
    }
    
    // 保持乱序，避免多线程同时纠错超长读数，使得内存报表
    // std::sort(read_ids_.begin(), read_ids_.end(), [this](int a, int b) { 
    //     return read_store_.GetSeqLength(a) > read_store_.GetSeqLength(b); 
    // });

    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>(read_ids_.begin(), read_ids_.end()));
    EstimateParameters();
}

std::unique_ptr<Dispatcher> PolDataset::GetDispatcher() {
    return std::unique_ptr<Dispatcher>(opts_.use_cache ? 
        (Dispatcher*)new GroupDispatcher(*this) : 
        (Dispatcher*)new SimpleDispatcher(*this));
}

void PolDataset::LoadReadIds() {
    std::unordered_set<Seq::Id> ids;    // Remove duplicate names
    if (!opts_.read_name_.empty()) {
        ids.insert(string_pool_.GetIdByStringUnsafe(opts_.read_name_));
    } else if (!opts_.read_name_fname_.empty()) {
        std::ifstream file(opts_.read_name_fname_);
        std::string line;
        while (std::getline(file, line)) {
            ids.insert(string_pool_.GetIdByStringUnsafe(line));
        }
    } else {
        // correct all reads in read file
    }
    read_ids_.assign(ids.begin(), ids.end());
    LOG(INFO)("read ids: %zd", read_ids_.size());
}


    


std::vector<Bed> PolDataset::CollectBedFromBam(const std::vector<Seq::Id>& read_ids) {
    assert(!opts_.rd_2_ref_fname_.empty());

    std::unordered_set<std::string> read_names;
    for (auto i : read_ids) {
        read_names.insert(string_pool_.QueryStringById(i));
    }

    std::vector<Bed> beds;

    std::mutex mutex;
    std::atomic<size_t> index { 0 };

    auto combine_func = [&mutex, &beds](const std::string& target, const std::vector<std::array<size_t,2>>& ranges) {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto& r : ranges) {
            beds.push_back({target, r[0], r[1]});
        }
    };

    auto work_func = [&read_names, combine_func, &index, this](int tid) {

        samFile *in = hts_open(opts_.rd_2_ref_fname_.c_str(), "r");
        assert(in != nullptr);
        
        bam_hdr_t *header = sam_hdr_read(in);
        assert(header != nullptr);
        auto idx = sam_index_load2(in, opts_.rd_2_ref_fname_.c_str(), (opts_.rd_2_ref_fname_+".bai").c_str());
        assert(idx != nullptr);

        std::vector<std::array<size_t, 2>> ranges;
        for (size_t tgtid = index.fetch_add(1); tgtid < header->n_targets; tgtid = index.fetch_add(1)) {

            hts_itr_t* itr = sam_itr_queryi(idx, tgtid, 0, header->target_len[tgtid]);
            bam1_t *b = bam_init1();
            
            while (sam_itr_next(in, itr, b) >= 0) {
                if (read_names.find(bam_get_qname(b)) != read_names.end()) {                   
                    size_t start = b->core.pos;
                    size_t end = b->core.pos + bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
                    
                    bool done = false;
                    for (auto &r : ranges) {
                        if (end > r[0] && start < r[1]) {
                            r[0] = std::min(start, r[0]);
                            r[1] = std::max(end, r[1]);
                            done = true;
                            break;
                        }
                    }
                    if (!done) {
                        ranges.push_back({start, end});
                    }
                }

            }
            bam_destroy1(b);

            std::sort(ranges.begin(), ranges.end(), [](const std::array<size_t, 2>& a, const std::array<size_t, 2> &b) {
                return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
            });

            std::vector<std::array<size_t, 2>> new_ranges;
            for (const auto& r0 : ranges) {
                
                bool done = false;
                for (auto &r : new_ranges) {
                    if (r0[1] > r[0] && r0[0] < r[1]) {
                        r[0] = std::min(r0[0], r[0]);
                        r[1] = std::max(r0[1], r[1]);
                        done = true;
                        break;
                    }
                }
                if (!done) {
                    new_ranges.push_back(r0);
                }

            }
            
            combine_func(header->target_name[tgtid], new_ranges);
            
        }
        bam_hdr_destroy(header);
        hts_close(in);
        
    };

        
    MultiThreadRun(opts_.thread_size, work_func);
    return beds;
}

void PolDataset::LoadOverlaps() {
    const std::string& fname = opts_.overlap_fname_;
    std::unordered_set<Seq::Id> ids(read_ids_.begin(), read_ids_.end());
    
    ol_store_.Load(fname, "", (size_t)opts_.thread_size, [this, &ids](Overlap &o) {
        bool rel = ids.empty() || ids.find(o.a_.id) != ids.end() || ids.find(o.b_.id) != ids.end();
        return rel && opts_.filter0_.Valid(o);
    });

    LOG(INFO)("Load %zd overlaps from file %s", ol_store_.Size(), fname.c_str());
   
}

void PolDataset::LoadMappings() {
    assert(!opts_.rd_2_ref_fname_.empty());

    if (!opts_.read_name_.empty() || !opts_.read_name_fname_.empty()) {
        LOG(INFO)("Start collecting bed");
        
        auto beds = CollectBedFromBam(read_ids_);
        LOG(INFO)("BED: %zd", beds.size());
        rd_2_ref_.LoadFileBam(opts_.rd_2_ref_fname_, [](const Overlap &o){return true;}, opts_.thread_size, beds);
    } else {
        rd_2_ref_.LoadFileBam(opts_.rd_2_ref_fname_, [](const Overlap &o){return true;}, opts_.thread_size);
    }

    mapping_.BuildIndex();
    LOG(INFO)("Load rd_2_ref size = %zd", rd_2_ref_.Size());
    

}

void PolDataset::LoadReads() {

    std::unordered_set<Seq::Id> ids;
    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        const Overlap& o = ol_store_.Get(i);
        ids.insert(o.a_.id);
        ids.insert(o.b_.id);
    }

    for (size_t i = 0; i < rd_2_ref_.Size(); ++i)  {
        const Overlap& o = rd_2_ref_.Get(i);
        ids.insert(o.a_.id);
    }
    read_store_.Load(opts_.rread_fname_, "", false, ids);
    LOG(INFO)("Load reads: %zd / %zd %zd", ids.size(), read_store_.Size(), rd_2_ref_.Size());
    
    if (read_ids_.empty()) {
        read_ids_.reserve(read_store_.Size());
        auto rs = read_store_.GetIdRange();
        for (Seq::Id i = rs[0]; i < (Seq::Id)rs[1]; ++i) {
            read_ids_.push_back(i);
        }
    }
}


std::vector<std::vector<Seq::Id>> PolDataset::GroupReadIds() const {
    
    std::vector<std::vector<Seq::Id>> clu_ids_;
    const float GOOD_ALIGNED_RATE = 0.6;
    std::unordered_map<Seq::Id, bool> done;

    for (auto i : read_ids_) {
        done[i] = false;
    }

    for (auto i : read_ids_) {
        if (done[i]) continue;
        
        if (clu_ids_.size() == 0 || clu_ids_.back().size() > 100) {
            clu_ids_.push_back(std::vector<Seq::Id>());
        }
        
        auto& curr = clu_ids_.back();
        curr.push_back(i);
        done[i] = true;

        auto gp = grouper_.Get(i);
        for (size_t ii = 0; ii < gp.Size(); ++ii) {
            auto ol = gp.Get(ii, 0);
            if (ol->AlignedLength() >= GOOD_ALIGNED_RATE*ol->TargetLength() || 
                ol->AlignedLength() >= GOOD_ALIGNED_RATE*ol->QueryLength()) {

                auto d = done.find(ol->GetOtherRead(gp.id).id);
                if (d != done.end() && !d->second) {
                    curr.push_back(d->first);
                    d->second = true;
                }
            }
        }
    }

    LOG(INFO)("Cluster reads: %zd", clu_ids_.size());
    return clu_ids_;
}

void PolDataset::EstimateParameters() {

    size_t count = std::min<size_t>(10, read_ids_.size());

    std::unordered_set<int> tests;

    std::default_random_engine e;
    std::uniform_int_distribution<int> u(0, read_ids_.size()-1);
    e.seed(time(0));
    
    while (tests.size() < count) {
        tests.insert(read_ids_[u(e)]);
    }

    auto aligner = ToolAligner::Create("edlib");

    double max_idt = 0.0;
    
    for (auto id : tests) {
        auto group = grouper_.Get(id);
        if (group.Empty()) continue;

        for (size_t i = 0; i < group.Size() && i < 10; ++i) {
            auto o = group.Get(i, 0);
            const auto& tread = o->GetRead(id);
            const auto& qread = o->GetOtherRead(id);
            std::vector<uint8_t> tseq = read_store_.GetSeq(tread.id).ToUInt8(tread.start, tread.end, false);
            std::vector<uint8_t> qseq = read_store_.GetSeq(qread.id).ToUInt8(qread.start, qread.end, !o->SameDirect());
            Alignment al;
            auto r = aligner->Align((const char*)&qseq[0], qseq.size(), (const char*)&tseq[0], tseq.size(), {0, qseq.size()}, {0, tseq.size()}, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
            if (r && al.AlignSize() > al.TargetSize() / 2) {
                if (al.Identity() > max_idt) {
                    max_idt = al.Identity();
                }

            }
        }
    }

    double min_idt = 0.0, min_lc_idt = 0.0;
    if (max_idt >= 0.99) {
        min_idt = 95;
        min_lc_idt = 90;
    } else if (max_idt >= 0.95) {
        min_idt = 90;
        min_lc_idt = 80;
    } else if (max_idt >= 0.85) {
        min_idt = 75;
        min_lc_idt = 65;
    } else {
        min_idt = 60;
        min_lc_idt = 50;
    }
    opts_.min_identity_ = opts_.min_identity_ < 0 ? min_idt : opts_.min_identity_;
    opts_.min_local_identity_ = opts_.min_local_identity_ < 0 ? min_lc_idt : opts_.min_local_identity_;
    LOG(INFO)("Estimate parameters(%0.02f): min_identity = %f min_local_identity = %f", max_idt, opts_.min_identity_, opts_.min_local_identity_);
}


PolDataset::OlGroup PolDataset::GetOverlaps(Seq::Id id) const {
    OlGroup group(id); 
    group.map = mapping_.QueryOverlaps(id);
    //LOG(INFO)("getols0: map %zd", group.map.size());
    
    
    group.ava = grouper_.GetRelatedOverlaps(id);
    DEBUG_printf("get_ols: map=%zd, ava=%zd\n", group.map.size(), group.ava.size());
    //LOG(INFO)("getols0: ava %zd", group.ava.size());

    group.ols.reserve(group.map.size() + group.ava.size());
    for (size_t i = 0; i < group.map.size(); ++i) {
        OlGroup::SetType(group.map[i], OlGroup::Type::MAP);
        if (opts_.filter0_.Valid(group.map[i])) {
            group.ols.push_back({1, i});
            //LOG(INFO)("ols(%s,%s): %s", read_store_.QueryNameById(group.map[i].a_.id).c_str(), 
            //    read_store_.QueryNameById(group.map[i].b_.id).c_str(), group.map[i].ToM4Line().c_str());

        }
    }
    for (size_t i = 0; i < group.ava.size(); ++i) {
        OlGroup::SetType(*group.ava[i], OlGroup::Type::AVA);
        group.ols.push_back({0, i});
    }
    
    DEBUG_printf("get_ols: size=%zd\n", group.ols.size());
    group.BuildIndex();
    return group;
}

void PolDataset::OlGroup::BuildIndex() {
    std::sort(ols.begin(), ols.end(), [this](const Index &ia, const Index &ib) { 

        const Overlap* a = Get(ia);
        const Overlap* b = Get(ib);
        const auto& r0 = a->GetOtherRead(id);
        const auto& r1 = b->GetOtherRead(id);

        return (r0.id < r1.id) ||
               (r0.id == r1.id && a->AlignedSize() > b->AlignedSize()) ||
               (r0.id == r1.id && a->AlignedSize() == b->AlignedSize() && a->SameDirect() && !b->SameDirect());
    });

    index.push_back({0, ols.size()});
    for (size_t i = 0; i < ols.size(); ++i) {
        const auto& r0 = Get(ols[index.back()[0]])->GetOtherRead(id);
        const auto& r1 = Get(ols[i])->GetOtherRead(id);

        if (r0.id != r1.id) {
            index.back()[1] = i;
            index.push_back({i, ols.size()});
        }
    }  

    for (size_t i = 0; i<index.size(); ++i) {
        DEBUG_printf("INDEX(%zd): %zd %zd\n", i, index[i][0], index[i][1]);
    }
}


void PolDataset::OlGroup::Sort(double opt_ohwt) {
    assert(!Empty());

    auto weights = GetWeight(opt_ohwt);

    std::sort(index.begin(), index.end(), [&weights](const std::array<size_t, 2> &a, const std::array<size_t,2> &b) {
        return weights[a[0]] > weights[b[0]];
    });

    
    for (size_t i = 0; i<index.size(); ++i) {
        DEBUG_printf("INDEX2(%zd): %zd %zd\n", i, index[i][0], index[i][1]);
    }
}

std::vector<double> PolDataset::OlGroup::GetWeight(double opt_ohwt) {
    assert(!Empty());
    size_t target_length = Get(ols[0])->GetRead(id).len;
    std::vector<double> cand_cov_wts (target_length+1);

    double wtsum = 0.0;
    for (size_t i = 0; i < Size(); ++i) {
        auto o = Get(i, 0);
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        double ohwt = opt_ohwt * o->identity_ / 100;
        double olwt = o->identity_ / 100;

        auto mr = o->MappingTo<2>(t, {0, q.len});
        auto start = std::max(0, mr[0] < mr[1] ? mr[0] : mr[1]);
        auto end =   std::min(t.len, mr[0] >= mr[1] ? mr[0] : mr[1]);
        // start -- t.start -- t.end -- end
        assert(start <= t.start && t.end <= end);

        cand_cov_wts[start]   += ohwt;
        cand_cov_wts[t.start] += (olwt - ohwt);
        cand_cov_wts[t.end]   -= (olwt - ohwt);
        cand_cov_wts[end]     -= ohwt;

        wtsum += olwt;
    }

    for (size_t i=1; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] += cand_cov_wts[i-1];
    }
    assert(std::abs(cand_cov_wts.back()) < 0.0000001);  // cand_cov_wts.back() == 0

    for (size_t i=0; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] = wtsum - cand_cov_wts[i];
    }

    std::vector<double> weights(ols.size(), 0.0);
    for (size_t i = 0; i < index.size(); ++i) {
        auto o = Get(ols[index[i][0]]);
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        double ohwt = opt_ohwt * o->identity_ / 100;
        double olwt = o->identity_ / 100;

        auto mr = o->MappingTo<2>(t, {0, q.len});
        auto start = std::max(0, mr[0] < mr[1] ? mr[0] : mr[1]);
        auto end =   std::min(t.len, mr[0] >= mr[1] ? mr[0] : mr[1]);
        // start -- t.start -- t.end -- end
        assert(start <= t.start && t.end <= end);

        double wt = std::accumulate(cand_cov_wts.begin()+start, cand_cov_wts.begin()+t.start, 0.0) * ohwt +
                    std::accumulate(cand_cov_wts.begin()+t.start, cand_cov_wts.begin()+t.end, 0.0) * olwt + 
                    std::accumulate(cand_cov_wts.begin()+t.end, cand_cov_wts.begin()+end, 0.0) * ohwt;

        weights[index[i][0]] = wt;
    }
    return weights;
}

}   // namespace fsa