#include "crr_dataset.hpp"

#include <random>

#include "crr_options.hpp"
#include "align/tool_aligner.hpp"

namespace fsa {

void CrrDataset::Load() {

    LoadReadIds();
    LoadOverlaps();
    LoadReads();
    
    // load mapping
    if (!opts_.rd_2_ref_fname_.empty()) {
        rd_2_ref_.Load(opts_.rd_2_ref_fname_);
        mapping_.BuildIndex();
    }

    std::sort(read_ids_.begin(), read_ids_.end(), [this](int a, int b) { 
        return read_store_.GetSeqLength(a) > read_store_.GetSeqLength(b); 
    });

    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>(read_ids_.begin(), read_ids_.end()));
    EstimateParameters();
}

std::unique_ptr<Dispatcher> CrrDataset::GetDispatcher() {
    return std::unique_ptr<Dispatcher>(opts_.use_cache ? 
        (Dispatcher*)new GroupDispatcher(*this) : 
        (Dispatcher*)new SimpleDispatcher(*this));
}

void CrrDataset::LoadReadIds() {
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

void CrrDataset::LoadOverlaps() {
    const std::string& fname = opts_.overlap_fname_;
    std::unordered_set<Seq::Id> ids(read_ids_.begin(), read_ids_.end());
    
    ol_store_.Load(fname, "", (size_t)opts_.thread_size, [this, &ids](Overlap &o) {
        bool rel = ids.empty() || ids.find(o.a_.id) != ids.end() || ids.find(o.b_.id) != ids.end();
        return rel && opts_.filter0_.Valid(o);
    });

    LOG(INFO)("Load %zd overlaps from file %s", ol_store_.Size(), fname.c_str());
   
}

void CrrDataset::LoadReads() {

    std::unordered_set<Seq::Id> ids;
    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        const Overlap& o = ol_store_.Get(i);
        ids.insert(o.a_.id);
        ids.insert(o.b_.id);
    }
    
    if (!opts_.rd_2_ref_fname_.empty()) {
        read_store_.Load(opts_.rread_fname_, "", false, ids);
    } else {
        read_store_.Load(opts_.rread_fname_, "");

    }
    
    if (read_ids_.empty()) {
        read_ids_.reserve(read_store_.Size());
        auto rs = read_store_.GetIdRange();
        for (Seq::Id i = rs[0]; i < rs[1]; ++i) {
            read_ids_.push_back(i);
        }
    }
}


std::vector<std::vector<Seq::Id>> CrrDataset::GroupReadIds() const {
    
    std::vector<std::vector<Seq::Id>> clu_ids_;
    const float GOOD_ALIGNED_RATE = 0.6;
    std::unordered_map<Seq::Id, bool> done;

    for (auto i : read_ids_) {
        done[i] = false;
    }

    for (auto i : read_ids_) {
        if (done[i]) continue;
        
        clu_ids_.push_back(std::vector<Seq::Id>());
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

void CrrDataset::EstimateParameters() {

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


CrrDataset::OlGroup CrrDataset::GetOverlaps(Seq::Id id) const {
    OlGroup group(id); 
    group.map = mapping_.QueryOverlaps(id);
    group.ava = grouper_.GetRelatedOverlaps(id);
; 
    group.ols.reserve(group.map.size() + group.ava.size());
    for (size_t i = 0; i < group.map.size(); ++i) {
        if (opts_.filter0_.Valid(group.map[i]))
            group.ols.push_back({1, i});
    }
    for (size_t i = 0; i < group.ava.size(); ++i) {
        group.ols.push_back({0, i});
    }

    group.BuildIndex();
    return group;
}

void CrrDataset::OlGroup::BuildIndex() {
    std::sort(ols.begin(), ols.end(), [this](const Index &ia, const Index &ib) { 

        const Overlap* a = Get(ia);
        const Overlap* b = Get(ib);
        const auto& r0 = a->GetOtherRead(id);
        const auto& r1 = b->GetOtherRead(id);

        return (r0.id < r1.id) ||
               (r0.id == r1.id && a->AlignedLength() > b->AlignedLength()) ||
               (r0.id == r1.id && a->AlignedLength() == b->AlignedLength() && a->SameDirect() && !b->SameDirect());
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
}


void CrrDataset::OlGroup::Sort(double opt_ohwt) {
    assert(!Empty());

    auto weights = GetWeight(opt_ohwt);

    std::sort(index.begin(), index.end(), [&weights](const std::array<size_t, 2> &a, const std::array<size_t,2> &b) {
        return weights[a[0]] > weights[b[0]];
    });
}

std::vector<double> CrrDataset::OlGroup::GetWeight(double opt_ohwt) {
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