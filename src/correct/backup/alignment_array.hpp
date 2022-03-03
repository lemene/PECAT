#ifndef FSA_ALIGNMENT_ARRAY_HPP
#define FSA_ALIGNMENT_ARRAY_HPP


#include <cassert>
#include <vector>
#include <list>
#include <bitset>
#include <map>
#include <set>

#include "graph_utils.hpp"
#include "aligner.hpp"

namespace fsa {


class AlignmentArray : public Corrector {
    static const int MAX_COV = 240;
public:
    struct Node {
        void AddSeq(int sid) { seqs.set(sid, true); }
        uint16_t count{ 0 };
        double score{ 0 };
        std::bitset<MAX_COV> seqs;
        Loc best {-1, -1, -1};
    };

    struct NodeGroup {
        Node& operator [] (size_t i) {
            return base[i];
        }
        static size_t Size() { return 5; }
        Node base[5];
    };


    struct Column {
        Column() { }
        NodeGroup& operator[] (size_t i) { return rows[i]; }
        size_t Size() const { return size; }
        void AddSeq(size_t id) { coverage++; seqs.set(id, true); }
        NodeGroup* rows { nullptr };
        size_t size { 0 };
        int index {0};
        int coverage {0};
        double weight {0};
        std::bitset<MAX_COV> seqs;
        //std::vector<NodeGroup> rows;
    };

    struct Score {
        int qid { -1 };
        int length { 0 };        
        int cross { 0 };         // number of positions where have 
        int t_in_cross { 0 };
        int q_in_cross { 0 };
        int q_t_one { 0 };
        int q_t_two { 0 };       //
        int distance { 0 };      // edit distance of two reads

        int t_start {0};
        int t_end { 0 };

        int range_q_t_one {0};
        int range_q_t_two {0};
        int range_q_t_total {0};
        int range_pos {-1};

    };


    AlignmentArray();
    
    void SetParameter(const std::string &name, const std::string &v);
    void SetParameter(const std::string &name, int v);
    void SetParameter(const std::string &name, double v);
    
    void Build(const DnaSeq& target, const std::array<size_t,2> &range, const std::vector<Alignment> &aligned, const std::vector<int> &qids);
    void Clear();


    Column& operator[] (size_t i) {
        return cols_[i];
    }


    size_t ColSize() { return cols_.size() - 1; };
    Column& GetCol(size_t i) { return cols_[i]; }

    size_t RowSize(size_t i) { return cols_[i+1].index - cols_[i].index; }
    NodeGroup& GetRow(size_t i, size_t j);

    size_t NodeSize() { return NodeGroup::Size(); }
    Node& GetNode(size_t c, size_t r, size_t n) { 
        return nodegroups_[cols_[c].index+r][n]; 
    }

    Node& operator[] (const Loc &loc) {
        return nodegroups_[cols_[loc.col].index+loc.row][loc.base];
    }

    void Save(const std::string &fname) const;

    Node* Get(const Loc& loc) {
        if (loc.col == -1) return nullptr;
        return &nodegroups_[cols_[loc.col].index+loc.row][loc.base];
    }

    Loc FindBestPath();
    Loc FindBestPath1();
    Loc FindBestPathWithCount();
    void Reconstruct(const Loc &end);
    void Consensus();

    const std::string& GetSequence() const { return sequence_; }

    void ComputeScore();

protected:
    void InitializeNodes(const DnaSeq &target, const std::vector<Alignment> &aligned);
    void AddTarget(const DnaSeq &target);
    void AddQuery(size_t qid, const Alignment &al);
    double LinkScoreCount(size_t col, size_t row) { return 1; } 
    double LinkScoreWeight(size_t col, size_t rowk);
protected:
    static DnaSerialTable2 Base2Num;

    std::vector<NodeGroup> nodegroups_;
    std::vector<Column> cols_;
    const DnaSeq* target_;
    std::array<size_t,2> range_;
    std::vector<std::array<int, 4>> ranges_;
    std::vector<Score> scores_;
    size_t window_size_ { 500 };
    int weight_base_ { 10 };

    double row_weight_ {0.5};

    std::string sequence_;
};


} // namespace fsa {
    
#endif // FSA_ALIGNMENT_ARRAY_HPP
