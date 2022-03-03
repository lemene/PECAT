#ifndef FSA_CONTIG_CORRECT_HPP
#define FSA_CONTIG_CORRECT_HPP

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "alignment_graph.hpp"
#include "utils/program.hpp"

namespace fsa {
using ArrayGraph = AlignmentGraph;

class ContigCorrect : public Program {
public:
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();
    virtual void CheckArguments();
    
protected:

    void LoadOverlaps(const std::string &fname);
    void LoadReadIds();
    void Correct();
    
    std::string OutputPath(const std::string &fname) { return output_directory_+"/"+fname; }

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
        Worker(ContigCorrect& owner) : owner_(owner) {
            graph_.SetParameter("score", owner.score_);
            graph_.SetParameter("min_coverage", owner.min_coverage_);
            aligner_.SetParameter("min_identity", owner.min_identity_);  
            aligner_.SetParameter("min_local_identity", owner.min_local_identity_);
            aligner_.SetParameter("aligner", owner_.aligner_);
        };
        ~Worker() {  }
        bool Correct(WindowJob &job);
        void CalculateWeight(Seq::Id tid,  const DnaSeq& target, const std::vector<const Overlap*> & cands, int offset, const std::array<int,2>& range);
        bool IsCoverageEnough(const std::vector<int> &cov);
        bool ExactFilter(const Alignment& r);
        bool GetAlignment(Seq::Id id, const Overlap* o, Alignment &al);
        void Clear() {graph_.Clear(); aligned_.clear(); corrected.clear(), scores_.clear(); }
        const std::string GetCorrected() const { return graph_.GetSequence(); }
    protected:
        ContigCorrect& owner_;
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
    int min_length_ {2000};
    int min_aligned_length_ { 10000 };
    double min_aligned_ratio_ { 0.50 };
    int max_aligned_length_ { 3000 };
    //int min_acceptable_length_ { 4000 };
    int max_overhang_ { 3000 };
    double max_overhang_rate_ { 0.30 };
    int min_coverage_ { 4 };
    double min_identity_ { 60 };
    double min_local_identity_ { 50 };

    std::string filter0_opts_ {"l=2000:al=2000:alr=0.50"};
    std::string filter1_opts_ {"l=2000:al=3000:alr=0.50:aal=5000:oh=2000:ohr=0.2"};
    Overlap::Filter filter0_;
    Overlap::Filter filter1_; 
 

    int max_number_ { 400 };    // MAX_COV - 1
    int coverage_ { 50 };
    double branch_score_ { 0.3 };
    int window_size_ { 50000 };
    int overlap_size_ { 500 };

    std::string read_name_ {""};
    std::string read_name_fname_ { "" };
    int thread_size_{ 4 };

    std::string aligner_ { "diff" };
    std::string score_ { "count" };
    std::string output_directory_ {"."};

    std::string overlap_fname_;
    std::string rread_fname_;
    std::string ctg_fname_;
    std::string cread_fname_;

    std::vector<Seq::Id> read_ids_;
    std::vector<ContigJob> jobs_;

    std::unordered_set<int> reads_;
    ReadStore read_store_;
    OverlapStore ol_store_{read_store_.GetStringPool() };
    std::unordered_map<int, std::unordered_map<int, std::vector<const Overlap*>>> groups_;
};

} // namespace fsa {
    
#endif // FSA_CONTIG_CORRECT_HPP

