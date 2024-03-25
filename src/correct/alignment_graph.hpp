#pragma once

#include <cassert>
#include <vector>
#include <list>
#include <map>
#include <set>

#include "corrector.hpp"

#include "aligner.hpp"

namespace fsa {

class ReadStore;
class CrrOptions;
class CrrDataset;

class AlignmentGraph : public Corrector {
public:

    struct Node {
        Node() {}

        void SortLinks() {
            for (size_t i=0; i<links.size(); ++i) {
                for (size_t j= i+1; j<links.size(); ++j) {
                    if (links[j].count > links[i].count) {
                        std::swap(links[i], links[j]);
                    }
                }
            }
        }
        uint16_t count{ 0 };
        double score{ 0 };
        const Link* best_link{ nullptr };

        std::vector<Link> links;
        MyBitSet seqs;
    };

    struct NodeGroup {
        Node& operator [] (size_t i) {
            assert(0 <= i && i < 5);
            return base[i];
        }
        const Node& operator [] (size_t i) const {
            assert(0 <= i && i < 5);
            return base[i];
        }
        size_t Size() const { return 5; }
        Node base[5];
    };


    struct Column {
        Column() {  }
        NodeGroup& operator[] (size_t i) { return rows[i]; }
        size_t Size() const { return rows.size(); }

        std::vector<NodeGroup> rows;
        uint coverage {0};
        int selected {0};
        double weight {0};
        MyBitSet queries;

    };

    struct Score {
        Score(Seq::Id id, size_t ts, size_t te, size_t colsize) 
         : qid(id), tstart(ts), tend(te), branches(colsize, nullptr) {}
        
        int qid { -1 };             // query read id
        size_t tstart { 0 };        // alignment in target
        size_t tend { 0 };
        int cross { 0 };            // number of positions where have 
        int all_cross { 0 };
        int t_in_cross { 0 };       // 
        int q_in_cross { 0 };       //
        int q_t_one { 0 };          //
        int q_t_two { 0 };          //

        struct BlockScore {
            int cross {0};
            int all_cross { 0 };
            int q_t_one { 0 };
            int q_t_two { 0 };
            double Weight() const {
                return cross == 0 ? 0.0 : (q_t_one - q_t_two)*1.0 / std::max(10, cross) ;
            }
        };
        std::vector<BlockScore> block_scores;

        double Weight() const {
            return cross == 0 ? 0.0 : (q_t_one - q_t_two)*1.0 / std::max(10, cross) ;
        }

        double WeightInGraph() const {
            return cross == 0 ? 0.0 : (q_t_one - q_t_two)*1.0 / std::max(10, cross) ;
        }

        double WeightInGraph(const std::array<double, 2>& r, const std::array<double,2>& wr) const {
            auto s = WeightInGraph();
            return r[0] == r[1] ? 0.5 : wr[0] + (s - r[0])/(r[1]-r[0]) * (wr[1] - wr[0]);
        }

        std::vector<const Link*> branches;
    };

    struct QueryInfos {

        void Clear() {
            scores_.clear();
            selected_.clear();
            windows.clear();
        }
        void SelectReads(int minsel);
        void SelectReads2(size_t minsel, const std::array<size_t,2> &range);
        void SelectReads3(size_t minsel, const std::array<size_t,2> &range);
        void SaveReadInfos(std::ostream& os, int tid, const ReadStore &rs) const ;

        double FindScoreThreshold();
        double FindScoreThreshold1();
        double FindScoreThreshold2();
        double FindScoreThreshold3(const std::vector<double>& score) const;
        size_t GetBlockSize() const;
        std::array<size_t, 3> GetWindowSize(const std::array<size_t, 2> &range) const;

        void SplitWindows(const std::array<size_t, 2> &range);
        std::vector<Score> scores_;
        std::unordered_set<int> selected_;
        
        std::vector<std::array<size_t, 2>> windows;
    };


    struct Segment {
        Loc end;
        Loc begin;
        int type = 0;
    };

    struct Options {
        int BranchThreshold(int cov);
        std::array<int, 2>      range {{10, 200}};
        std::array<float, 2>    rate {{0.4, 0.25}};
    };

    struct BlockRange {

    };

    AlignmentGraph(const CrrOptions& opts, const CrrDataset& ds);

    void SetParameter(const std::string &name, const std::string &v);
    void SetParameter(const std::string &name, double v);

    void ParseScoreParamter(const std::string &opts);

    
    Column& operator[] (size_t i) { return cols[i]; }
    Node& operator[] (const Loc &loc) { return cols[loc.col][loc.row][loc.base]; }

    void Build(const DnaSeq& target, const std::array<size_t,2> &range, const std::vector<Alignment> &aligned);
    void BuildCol(size_t s, size_t e);
    void BuildNode(size_t s, size_t e);
    void Clear();

    Loc Locate(const Node& node);

    Node* Get(const Loc& loc) {
        if (loc.col < 0 || loc.col >= (int)cols.size()) return nullptr;
        if (loc.row < 0 || loc.row >= (int)cols[loc.col].Size()) return nullptr;
        if (loc.base < 0 || loc.base >= (int)cols[loc.col].rows[loc.row].Size()) return nullptr;
        return &cols[loc.col].rows[loc.row].base[loc.base];
    }

    Segment FindBestPathBasedOnCount();
    Segment FindBestPathBasedOnWeight();

    std::vector<Segment> SplitSegment(const Loc &end);
    std::vector<Segment> SplitSegment2(const Loc &end);
    void Reconstruct(const std::vector<Segment>& segs);
    void Consensus();

    std::string ReconstructSimple(const Segment& seg);
    std::string ReconstructComplex(const Segment& seg);

    // std::vector<std::string> RestoreSegment(const Segment &seg);

    const std::string& GetSequence() const { return sequence_; }
    const std::string& GetQuality() const { return quality_; }
    const std::array<size_t, 2>& GetTrueRange() const { return true_range_; }

    void ComputeSimilarity4();
    struct ImportantBranch {
        size_t c;   // column
        struct LinkCol{
            const Link* l;
            uint8_t     r;    // 0-5
        } ;
        std::array<LinkCol,2> links;
        bool valid { true };
    };

    std::vector<int> ExtendLeft(const Loc& start, size_t n, const MyBitSet& seqs);
    std::vector<int> ExtendRight(const Loc& right, size_t n, const MyBitSet& seqs);
    std::vector<ImportantBranch> CollectImportantBranches();
    void VerifyImportantBranches1(std::vector<ImportantBranch>& cands);
    void VerifyImportantBranchesByVariants(std::vector<ImportantBranch>& cands);
    bool VerifiyImportantBranch(const std::vector<ImportantBranch>& cands, size_t start, size_t end);
    void VerifyImportantBranchesByDensity(std::vector<ImportantBranch> &cands);
    void VerifyImportantSitesByConsistent(std::vector<ImportantBranch> &cands, double th, std::unordered_set<size_t>& removed);
    void ReactivateImportantSitesByConsistent(std::vector<ImportantBranch> &cands, double th, const std::unordered_set<size_t>& removed);
    
    std::vector<const Link*> CollectLinks(size_t i);

    void SaveGraph(const std::string& fname, size_t s, size_t e) const;
    void SaveReadInfos(std::ostream &os, int tid, const ReadStore& rs) const { query_infos_.SaveReadInfos(os, tid, rs); }
protected:
    void AddTarget(const DnaSeq &target, const std::array<size_t, 2> &range);
    void AddQuery(size_t qid, size_t query_start, const std::string &aligned_query,  size_t target_start, const std::string &aligned_target);
    double LinkScoreCount(size_t col, size_t row, const Link& link);
    double LinkScoreWeight(size_t col, size_t row, const Link& link);
    bool IsSimpleColumn(size_t col);
protected:
    static DnaSerialTable2 Base2Num;

    Seq::Id tid_;
    std::vector<Column> cols;
    const DnaSeq* target_;
    std::array<size_t,2> range_;
    std::array<size_t,2> true_range_ {{0, 0}};
    std::array<double, 2> score_range_;
    std::array<double, 2> weight_range_ { {0.4, 0.8 }};
    std::array<double, 3> branch_score_ { {0.4, 0.3, 0.8} };
    double reduction_ { 0.8 };
    int min_selected { 20 };
    int max_bubble_length_ { 100 };
    Options opts_;

    const CrrOptions &sopts_;
    const CrrDataset &dataset_;

    std::vector<Tag> tags_;

    std::string sequence_;
    std::string quality_;
    double (AlignmentGraph::*LinkScore)(size_t, size_t, const Link&) { nullptr };
    void (AlignmentGraph::*VerifyImportantBranches0)(std::vector<ImportantBranch>& cands) { nullptr };
    QueryInfos query_infos_;
};


} // namespace fsa {
