#include "crr_dataset.hpp"

#include <random>

#include "crr_options.hpp"
#include "align/tool_aligner.hpp"

namespace fsa {

void CrrDataset::Load() {
    LoadReadIds();
    LoadOverlaps();
    LoadReads();

    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>(read_ids_.begin(), read_ids_.end()));

    if (opts_.use_cache) GroupReadIds();
    EstimateParameters();
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
    read_store_.Load(opts_.rread_fname_, "", false, ids);
    
    if (read_ids_.empty()) {
        read_ids_.assign(ids.begin(), ids.end());
    }
}


void CrrDataset::GroupReadIds() {
    const int CLU_SIZE = 100;
    std::unordered_map<Seq::Id, bool> done;

    for (auto i : read_ids_) {
        done[i] = false;
    }

    std::sort(read_ids_.begin(), read_ids_.end(), [this](int a, int b) { 
        return read_store_.GetSeqLength(a) > read_store_.GetSeqLength(b); 
    });

    for (auto i : read_ids_) {
        if (done[i]) continue;

        if (clu_ids_.size() == 0 || clu_ids_.back().size() >= CLU_SIZE) {
            clu_ids_.push_back(std::vector<Seq::Id>());
        }

        auto& curr = clu_ids_.back();
        curr.push_back(i);
        done[i] = true;

        for (size_t idx = 0; idx < curr.size() && curr.size() <= CLU_SIZE; idx++) {
            auto gp = grouper_.Get(curr[idx]);
            for (size_t i = 0; i < gp.Size(); ++i) {
                auto ol = gp.Get(i, 0);
                auto d = done.find(ol->GetOtherRead(gp.id).id);
                if (d != done.end() && !d->second) {
                    curr.push_back(d->first);
                    d->second = true;
                }
            }
        }
    }

    LOG(INFO)("Cluster reads: %zd", clu_ids_.size());
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


}   // namespace fsa