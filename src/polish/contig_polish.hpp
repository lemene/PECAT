#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "../correct/alignment_graph.hpp"
#include "utils/program.hpp"
#include "pol_dataset.hpp"
#include "pol_options.hpp"

namespace fsa {
using ArrayGraph = AlignmentGraph;

class ContigPolish : public Program {
public:
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();
    virtual void CheckArguments() { opts_.CheckArguments(); }
    
protected:

    void LoadOverlaps(const std::string &fname);
    void Correct();
    void CalcCoverage();
    
    struct ContigJob;
    struct WindowJob {
        WindowJob(ContigJob *w, int s, int e) : owner(w), start(s), end(e) {}

        std::vector<const Overlap*> GetOverlaps() {
            std::vector<const Overlap*> ols;
            for (auto i : owner->overlaps) {
                for (auto o : i.second) {
                    auto& r = o->GetRead(owner->tid);
                    auto s = std::max(r.start, start);
                    auto e = std::min(r.end, end);
                    if (e > s + 2000) {
                        ols.push_back(o);
                    }
                }
            }
            return ols;
        }
        Seq::Id GetTId() { return owner->tid; }
        ContigJob *owner { nullptr};
        int start, end;
        std::string seq;
        std::atomic<bool> done { false};

    };

    struct ContigJob {
        ContigJob(Seq::Id id, size_t len, const std::unordered_map<int, std::vector<const Overlap*>>& ols, size_t wsize, size_t osize);
        // 将各个窗口的数据拼接起来
        std::string GetSeq() const;
        
        bool IsDone() const {
            for (const auto& w : windows) {
                if (!w->done) return false;
            }
            return true;
        }

        bool Savable() {
            return IsDone() && !saved.exchange(true);
        }


        Seq::Id tid;    // target id
        size_t tlen;
        const std::unordered_map<int, std::vector<const Overlap*>>& overlaps;
        std::vector<std::shared_ptr<WindowJob>> windows;
        size_t win_size;
        size_t ovl_size;
        std::atomic<bool> saved { false };
    };


    class Worker {
    public:
        Worker(ContigPolish& owner) : owner_(owner), graph_(owner.opts_.min_coverage, owner_.dataset_.GetStringPool()) {
            aligner_.SetParameter("aligner", owner_.opts_.aligner_);
        };
        ~Worker() {  }
        bool Correct(WindowJob &job);
        void CalculateWeight(Seq::Id tid,  const DnaSeq& target, const std::vector<const Overlap*> & cands, int offset, const std::array<int,2>& range);
        bool IsCoverageEnough(const std::vector<int> &cov);
        bool ExactFilter(const Alignment& r);
        bool GetAlignment(Seq::Id id, const Overlap& ol, Alignment &al, int ctgstart);
        void GetAlignmentFromCigar(Seq::Id tid, const Overlap &ol, Alignment &al);
        void Clear() {graph_.Clear(); aligned_.clear(); corrected.clear(), scores_.clear(); }
        const std::string GetCorrected() const { return graph_.GetSequence(); }
    protected:
        ContigPolish& owner_;
        ArrayGraph graph_;
        Aligner aligner_;
        std::vector<Alignment> aligned_;
        std::string corrected;
        std::vector<ArrayGraph::Score> scores_;
    public:
        std::array<int,3> counts_ {{0, 0, 0}}; // for debug
    };   
    friend class Worker;

protected:
 


    std::vector<ContigJob> jobs_;
    
    PolOptions opts_;
    PolDataset dataset_ { opts_ };
};

} // namespace fsa {
    
