#pragma once

#include "minimizer_counter.hpp"

#include <numeric> // accumulate
#include <inttypes.h>

#include "mkseq_store.hpp"
namespace fsa {

class MkseqStore;

class MinimizerGraph {
public:

class RawNode {
    public:
        RawNode(const Minimizer* m, uint8_t d) : mkmer_(m),dir_(d) {}
        bool operator==(const RawNode& b) const {
            return  Hash() == b.Hash() && Direction()  == b.Direction();
        }
        bool operator!=(const RawNode& b) const {
            return !(*this == b);
        }
        bool operator<(const RawNode& b) const {
            return Hash() < b.Hash() || (Hash() == b.Hash() && Direction() < b.Direction());
        }

        std::string ToString() const {
            char buf[128];
            snprintf(buf, sizeof(buf), "%016" PRIX64 "_%d", Hash(), Direction());
            return std::string(buf);
        }
        uint8_t Direction() const { return dir_ ^ mkmer_->dir; }
        uint64_t Hash() const { return mkmer_->hash; }
        
        
    protected:
        const Minimizer* mkmer_;
        uint8_t dir_;
    };

    class RawEdge {
    public:
        RawEdge(size_t seqid, size_t pos, uint8_t dir) {
            seq_id_ = seqid;
            pos_ = pos;
            dir_ = dir;
        }
        const Minimizer& GetInNode(const MkseqStore& mkseqs) const {
            if (dir_ == 0) 
                return mkseqs.GetMinimizer(seq_id_, pos_);
            else  
                return mkseqs.GetMinimizer(seq_id_, pos_+1);
        }
        const Minimizer& GetOutNode(const MkseqStore& mkseqs) const {
            if (dir_ == 0) 
                return mkseqs.GetMinimizer(seq_id_, pos_+1);
            else
                return mkseqs.GetMinimizer(seq_id_, pos_);
        }

        RawNode InNode(const MkseqStore& mkseqs) const {
            if (dir_ == 0) 
                return RawNode(&mkseqs.GetMinimizer(seq_id_, pos_), dir_);
            else    
                return RawNode(&mkseqs.GetMinimizer(seq_id_, pos_+1), dir_);
        }

        RawNode OutNode(const MkseqStore& mkseqs) const {
            if (dir_ == 0) 
                return RawNode(&mkseqs.GetMinimizer(seq_id_, pos_+1), dir_);
            else    
                return RawNode(&mkseqs.GetMinimizer(seq_id_, pos_), dir_);
        }


        bool Equal(const MkseqStore& mkseqs, const RawEdge& b) const {
            const auto& a0 = InNode(mkseqs);
            const auto& a1 = OutNode(mkseqs);
            const auto& b0 = b.InNode(mkseqs);
            const auto& b1 = b.OutNode(mkseqs);
            return a0 == b0 && a1 == b1;
        }
        size_t seq_id_;
        size_t pos_;
        uint8_t dir_;
    };

    
    class Node {
    public:
        struct Hash {
            size_t operator()(const Node& n) const {
                std::hash<uint8_t> uint8_hash;
                return n.hash_ ^ (uint8_hash(n.dir_) << 1);
            }
        };
        Node(uint64_t h, uint8_t d, uint8_t p) : hash_(h),dir_(d), pos_(p) {}
        bool operator==(const Node& b) const {
            return hash_ == b.hash_ && dir_ == b.dir_;
        }
        bool operator!=(const Node& b) const {
            return !(*this == b);
        }
        bool operator<(const Node& b) const {
            return hash_ < b.hash_ || (hash_ == b.hash_ && dir_ < b.dir_);
        }
    protected:
        uint64_t hash_;
        uint8_t dir_;
        uint8_t pos_;
    };

    class Edge {
    public:
        Edge(size_t s, size_t e) : start_(s), end_(e) {

        }
        size_t Count() const { return end_ - start_; }
        // Node InNode(const MinimizerGraph &graph) const { 
        //     auto rnode = graph.redge_list_[start_]->GetInNode(*graph.mkseqs_);
        //     return Node(mkseqs_.GetMinimizer(mkseqs_.Get(redge_list_[start_]->seq_id_).Get(redge_list_[start_]->pos_ )));
        // }
        // Node OutNode(const MinimizerGraph &graph) const { 
        //     return Node(graph.redge_list_[start_]->GetOutNode(*graph.mkseqs_));
        // }
        size_t start_;
        size_t end_;
        // TODO add path information
        size_t position_ { 0 };  // which path this edge belongs to

    };


public:
    MinimizerGraph(){}
    ~MinimizerGraph() {
        for (auto e : redge_list_) {
            delete e;
        }
        redge_list_.clear();
        for (auto e : edge_list_) {
            delete e;
        }
        edge_list_.clear();
    }
    void Build(const MkseqStore& mkseqs);
    void Save(const std::string& fname) const;
protected:
    const MkseqStore *mkseqs_;
    std::vector<RawEdge*> redge_list_;
    std::vector<Edge*> edge_list_;
    std::unordered_map<Node, size_t, Node::Hash> node_2_edge;
};

}