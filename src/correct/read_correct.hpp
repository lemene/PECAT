#pragma once

#include <string>
#include <vector>
#include <mutex>

#include "alignment_graph.hpp"
#include "utils/program.hpp"
#include "align/alignment_cache.hpp"
#include "crr_options.hpp"
#include "crr_dataset.hpp"

namespace fsa {

class ReadCorrect : public Program {
public:
    ReadCorrect();
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();
    virtual void CheckArguments() { opts_.CheckArguments(); }

protected:
    void Correct();
    void SaveCRead(std::ostream &os, int tid, const std::string &cread, const std::array<size_t, 2> &range);

    // 
    struct StatInfo {
        void Merge(const StatInfo si) {
            aligns[0] += si.aligns[0];
            aligns[1] += si.aligns[1];
        }
        std::array<int, 4> aligns {{0,0,0,0}}; // 统计执行详细比对的测试 all, succ, fails
    };

    class Worker {
    public:
        Worker(ReadCorrect& owner) : owner_(owner), graph_(owner.opts_, owner_.dataset_) {
            graph_.SetParameter("score", owner.opts_.score_);
            aligner_.SetParameter("min_identity", owner.opts_.min_identity_);  
            aligner_.SetParameter("min_local_identity", owner.opts_.min_local_identity_);
            aligner_.SetParameter("aligner", owner_.opts_.aligner_);
        };
        ~Worker() {  }
        bool Correct(int id, bool uc=true);
        void CalculateWeight(Seq::Id tid,  const DnaSeq& target, std::vector<std::tuple<const Overlap*, double, size_t>> & cands, double opt_ohwt);
        bool IsCoverageEnough(const std::vector<int> &cov);
        bool ExactFilter(const Alignment& r);
        bool ExactFilter(const Alignment& r, const std::array<size_t,2>& trange);
        bool GetAlignment(Seq::Id id, const Overlap* o, bool uc, Alignment &al);
        void Clear() {graph_.Clear(); aligned_.clear(); corrected.clear(); }
        void ClearCache() { return cache_.Clear(); }
        void ResetCache(const std::vector<Seq::Id> &ids, size_t size) { return cache_.Reset(ids, size); }
        const std::string GetCorrected() const { return graph_.GetSequence(); }
        const std::array<size_t, 2> & GetTrueRange() const { return graph_.GetTrueRange(); }
        void SaveReadInfos(std::ostream& os, int tid, const ReadStore &rd) { graph_.SaveReadInfos(os, tid, rd); }
   
        StatInfo stat_info;
    protected:
        ReadCorrect& owner_;
        AlignmentGraph graph_;
        Aligner aligner_;
        std::vector<Alignment> aligned_;
        std::string corrected;
        AlignmentCache cache_;
    };   
    friend class Worker;

protected:
    void CollectWorkerInfo(const Worker &w) {
        stat_info_.Merge(w.stat_info);
    }
    void CollectWorkerInfo(const Worker &w, std::mutex& m) { 
        std::lock_guard<std::mutex> lock(m); 
        CollectWorkerInfo(w);
    }
    void Report() const {
        LOG(INFO)("alignment %d %d", stat_info_.aligns[0], stat_info_.aligns[1]);
    }
protected:
    CrrOptions opts_;
    CrrDataset dataset_ { opts_ };

    StatInfo stat_info_;

};

} // namespace fsa {


