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
        std::bitset<MAX_COV> seqs;
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
        std::bitset<MAX_COV> queries;

    };

    struct Score {
        Score(Seq::Id id, size_t colsize) : qid(id), branches(colsize, nullptr) {}
        int qid { -1 };             // query read id
        int cross { 0 };            // number of positions where have 
        int t_in_cross { 0 };       // 
        int q_in_cross { 0 };       //
        int q_t_one { 0 };          //
        int q_t_two { 0 };          //

        double Weight() const {
            return cross == 0 ? 0.0 : (q_t_one - q_t_two)*1.0 / cross ;
        }

        double WeightInGraph() const {
            return cross == 0 ? 0.0 : (q_t_one - q_t_two)*1.0 / cross ;
        }

        double WeightInGraph(const std::array<double, 2>& r, const std::array<double,2>& wr) const {
            auto s = WeightInGraph();
            return r[0] == r[1] ? 0.5 : 0.1 + (s - r[0])/(r[1]-r[0]) * (wr[1] - wr[0]);
        }

        std::vector<const Link*> branches;
    };

    struct QueryInfos {

        void SelectReads(int minsel);
        void SaveReadInfos(std::ostream& os, int tid, const ReadStore &rs) const ;

        double FindScoreThreshold();
        double FindScoreThreshold1();
        double FindScoreThreshold2();
        std::vector<Score> scores_;
        std::unordered_set<int> selected_;
    };


    struct Segment {
        Loc end;
        Loc begin;
        int type = 0;
    };

    struct Options {
        int BranchThreshold(int cov);
        std::array<int, 2>      range {{10, 200}};
        std::array<float, 2>    rate {{0.5, 0.2}};
    };

    AlignmentGraph();

    void SetParameter(const std::string &name, const std::string &v);
    void SetParameter(const std::string &name, int v);
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
        if (loc.col == -1) return nullptr;
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

    std::vector<std::string> RestoreSegment(const Segment &seg);

    const std::string& GetSequence() const { return sequence_; }
    const std::string& GetQuality() const { return quality_; }


    void ComputeSimilarity();
    void SelectReads();

    std::vector<const Link*> CollectLinks(size_t i);
    std::vector<const Link*> CollectLinks1(size_t i);
    std::vector<const Link*> CollectLinks2(size_t i);

    Link* AllocLink(size_t sz);
    NodeGroup* AllocNodeGroup(size_t sz);
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

    std::vector<Column> cols;
    std::vector<const Link*> branches_;     // target每个base所在的分支，nullptr表示该base不在分支上
    std::vector<bool> branch_flags_;        // base是否有分支
    const DnaSeq* target_;
    std::array<size_t,2> range_;
    std::array<double, 2> score_range_;
    std::array<double, 2> weight_range_ { {0.2, 0.8 }};
    std::array<double, 3> branch_score_ { {0.4, 0.3, 0.8} };
    double reduction_ { 0.8 };
    int min_selected { 20 };
    int max_bubble_length_ { 100 };
    int min_coverage_ { 4 };
    Options opts_;
    std::vector<Tag> tags_;

    std::string sequence_;
    std::string quality_;
    double (AlignmentGraph::*LinkScore)(size_t, size_t, const Link&) { nullptr };
    const Node* (AlignmentGraph::*FindBestPath)() {nullptr };

    QueryInfos query_infos_;
};


} // namespace fsa {
