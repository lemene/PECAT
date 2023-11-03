#pragma once

#include <string>
#include <vector>
#include <mutex>

#include "overlap_store.hpp"
#include "read_store.hpp"
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
    virtual void CheckArguments();
    
    // the options for selecting candiate overlaps
    struct CandidateOptions {
        CandidateOptions(const std::string &str) { From(str); }
        void From(const std::string &str);
        std::string ToString() const;

        bool IsEndCondition(const std::vector<int> &cov, size_t number, size_t fail) const {
            return (int)fail >= failures || IsEnough(cov);
        }
        bool IsEnough(const std::vector<int> &cov) const ;
        bool IsEnough(size_t number) const { return false; }

        double percent { 0.95 };             // p Percentage of filled matrix
        double overhang_weight   { 0.0 };                // w overhang的比重
        int failures { 10 };               // f 连续失败次数
        int max_number { 200 };             // 
        int coverage { 80 };                    // 需要多少层数据
    };

protected:

    void LoadReads();
    void LoadOverlaps(const std::string &fname);
    void LoadReadIds();
    void Correct();
    void SaveCRead(std::ostream &os, int tid, const std::string &cread, const std::array<size_t, 2> &range);
    

    void GroupReadIds();
    void EstimateParameters();


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

    struct Progress {
        Progress(ReadCorrect& rc) : owner(rc) {}
        Seq::Id Get() {
            auto curr = index.fetch_add(1);
            if (curr % log_block_size == 0) {
             LOG(INFO)("done %zd, all %zd %zd", curr, owner.read_ids_.size(), GetMemoryUsage());
            }
            return curr < owner.read_ids_.size() ? owner.read_ids_[curr] : Seq::NID;
        }
        
        size_t Get(std::vector<Seq::Id> &ids) {
            auto curr = index.fetch_add(1);
            if (curr < owner.group_ticks.size()-1) {
                if (owner.group_ticks[curr] - last_log >= log_block_size) {
                    last_log = owner.group_ticks[curr];
                    LOG(INFO)("done %zd/%zd %zd", owner.group_ticks[curr], owner.grouped_ids_.size(), GetMemoryUsage());
                }
                assert(owner.group_ticks[curr+1]-owner.group_ticks[curr] <= ids.size());
                for (size_t i=owner.group_ticks[curr]; i<owner.group_ticks[curr+1]; ++i) {
                    ids[i-owner.group_ticks[curr]] = owner.grouped_ids_[i];
                }
                return owner.group_ticks[curr+1] - owner.group_ticks[curr];
            } else {
                return (size_t)0;
            }
        };

        std::atomic<size_t> index {0};
        const ReadCorrect& owner;
        const size_t log_block_size { 5000 };
        size_t last_log { 0 };  
    };
    friend class Progress;

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
    Overlap::Filter filter0_;
    Overlap::Filter filter1_; 
    CandidateOptions cands_opts_ { opts_.cands_opts_str_ };

    std::vector<Seq::Id> read_ids_;
    std::vector<Seq::Id> grouped_ids_;
    std::vector<size_t> group_ticks;
    size_t  group_size  = 1000;

    std::unordered_set<int> reads_;
    StringPool string_pool_;
    ReadStore read_store_ {string_pool_};
    OverlapStore ol_store_{string_pool_ };

    OverlapGrouper grouper_ { ol_store_ };

    StatInfo stat_info_;

};

} // namespace fsa {


