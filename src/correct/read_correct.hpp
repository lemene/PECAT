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


    bool ExactFilter(const Alignment& r);
    bool ExactFilter(const Alignment& r, const std::array<size_t,2>& trange);

    // 
    struct StatInfo {
        void Merge(const StatInfo si) {
            total += si.total;
            cache += si.cache;
            succ += si.succ;
        }
        void Clear() {
            total = 0;
            cache = 0;
            succ = 0;
        }
        void Report() const {
            LOG(INFO)("alignment %d %d %d", total, cache, succ);
        }
        int total { 0 };
        int cache { 0 };
        int succ { 0 };
    };

    class Worker {
    public:
        Worker(ReadCorrect& owner) : owner_(owner), graph_(owner.opts_, owner_.dataset_) {
            graph_.SetParameter("score", owner.opts_.score_);
            aligner_.SetParameter("aligner", owner_.opts_.aligner_);
        };
        ~Worker() {  }
        bool GetAlignment(Seq::Id id, const Overlap* o, Alignment &al);
        void Clear() {graph_.Clear(); aligned_.clear(); corrected.clear(); }
        void ClearCache() { return cache_.Clear(); }
        void ResetCache(const std::vector<Seq::Id> &ids, size_t size) { return cache_.Reset(ids, size); }
        const std::string GetCorrected() const { return graph_.GetSequence(); }
        const std::array<size_t, 2> & GetTrueRange() const { return graph_.GetTrueRange(); }
        void SaveReadInfos(std::ostream& os, int tid, const ReadStore &rd) { graph_.SaveReadInfos(os, tid, rd); }
        std::vector<Alignment> CheckLocalDistance0(const std::vector<Alignment>& als);
        std::vector<std::array<size_t, 2>> GroupPositions(const std::vector<size_t> &sorted_positions);
        bool Cigar2Alignment(Seq::Id tid, const Overlap* ol, Alignment &al);
   
        StatInfo stat_info;

        ReadCorrect& owner_;
        AlignmentGraph graph_;
        Aligner aligner_;
        std::vector<Alignment> aligned_;
        std::string corrected;
        AlignmentCache cache_;
    };   
    friend class Worker;
    
    bool Correct(Seq::Id tid, Worker& wrk);
    Alignment GetAligmentBetweenTwoReads(Seq::Id tid, const CrrDataset::OlGroup& group, size_t ig, Worker& wrk);
protected:
    CrrOptions opts_;
    CrrDataset dataset_ { opts_ };
};

} // namespace fsa {


