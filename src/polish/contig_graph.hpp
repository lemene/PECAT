#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <cassert>

#include "contig_fragment.hpp"
#include "../utils/logger.hpp"

namespace fsa {

class MatchInfo;

class ContigGraph {
public:
    class Node;
    class Edge {
    public:
        Edge(Node* n0, Node* n1) : node0_(n0), node1_(n1) {}
        const Node* OtherNode(const Node* n) const { 
            if (n == node0_) {
                return node1_; 
            } else if (n == node1_) {
                return node0_; 
            } else {
                assert(!"Invalid node");
                return nullptr; // Should never reach here
            }
        }

        Node* OtherNode(Node* n) {
            return node1_; 
        }
        void AddMatch(const std::vector<const MatchInfo*> &mi0, const std::vector<const MatchInfo*>& mi1) {
            match0_.insert(match0_.end(), mi0.begin(), mi0.end());
            match1_.insert(match1_.end(), mi1.begin(), mi1.end());
        }

    protected:
        Node* node0_ {nullptr};
        Node* node1_ {nullptr};
        std::vector<const MatchInfo*> match0_;
        std::vector<const MatchInfo*> match1_;
    };
    class Node {
    public:
        Node(uint32_t id, const ContigFragment* frg, uint8_t end)
         : id_(id), frg_(frg), end_(end) {}

        void AddEdge(Edge* e) {
            edges_.push_back(e);
        }    

        uint32_t Id() const { return id_; }
        const ContigFragment* Fragment() const { return frg_; }
        const std::vector<Edge*>& Edges() const { return edges_; }
        uint8_t End() const { return end_; }
 
        Edge* FindEdge(const Node* n) {
            for (auto& e : edges_) {
                if (e->OtherNode(this) == n) return e;
            }
            return nullptr;
        }

    protected:
        uint32_t id_ { 0 }; // Node ID
        const ContigFragment* frg_;
        uint8_t end_ {0}; // 0: start, 1: end
        std::vector<Edge*> edges_; 
    };

    class Chain {
    public:
        Chain(std::vector<const Node*>&& paths) : paths_(std::move(paths)) {
            // Initialize other members if needed
        }
        std::string Polish();
    protected:
        std::vector<const Node*> paths_;
    };

public:
    ContigGraph() = default;
    ContigGraph(const ContigGraph&) = delete;
    ContigGraph& operator=(const ContigGraph&) = delete;
    ~ContigGraph();

    void AddFragment(const std::vector<ContigFragment>& frags) {
        fragments_.insert(fragments_.end(), frags.begin(), frags.end());
    }

    void Build();
    void Simplify();
    std::vector<Chain> GetChains();
    void Save(const std::string& fname) const ;

    Node* GetPairNode(Node* n) {
        assert(n != nullptr);
        return nodes_.at(GetPairNodeId(n->Id())); // flip last bit to get the pair node
    }
    const Node* GetPairNode(const Node* n) const {
        assert(n != nullptr);
        LOG(INFO)("PairNode %08X %08X", n->Id(), GetPairNodeId(n->Id()));
        return nodes_.at(GetPairNodeId(n->Id())); // flip last bit to get the pair node
    }
    Node* GetNode(uint32_t id) {
        assert(nodes_.find(id) != nodes_.end());
        return nodes_[id];
    }
    const Node* GetNode(uint32_t id) const {
        assert(nodes_.find(id) != nodes_.end());
        return nodes_.at(id);
    }
    
    uint32_t GetPairNodeId(uint32_t id) const {
        return id ^ 1; // flip last bit to get the pair node
    }

    uint32_t MakeNodeId(uint32_t id, uint8_t end) {
        assert(id >= 0 && id < fragments_.size() && (end == 0 || end == 1));
        return (id << 1) + end;
    }

    uint64_t MakeEdgeId(uint32_t n0, uint32_t n1) {
        if (n0 < n1) {
            return ((uint64_t)n0 << 32) + n1;
        } else {
            return ((uint64_t)n1 << 32) + n0;
        }
    }

protected:
    std::vector<ContigFragment> fragments_;
    std::unordered_map<uint32_t, Node*> nodes_;
    std::unordered_map<uint64_t, Edge*> edges_; 
};

}