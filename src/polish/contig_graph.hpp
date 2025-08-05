#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <cassert>

#include "contig_fragment.hpp"

namespace fsa {

class ContigGraph {
public:
    class Node;
    class Link {
    public:
        Link(Node* to) : to_(to) {}
        Node* To() const { return to_; }  
    protected:
        Node* to_ {nullptr};
    };
    class Node {
    public:
        Node(uint32_t id, const ContigFragment* frg, uint8_t end)
         : id_(id), frg_(frg), end_(end) {}

        void AddLink(Node* to) {
            link_to_.push_back(Link(to));
        }
        uint32_t Id() const { return id_; }
        const ContigFragment* Fragment() const { return frg_; }
        const std::vector<Link>& Links() const { return link_to_; }

    protected:
        uint32_t id_ {0}; // Node ID
        const ContigFragment* frg_;
        uint8_t end_ {0}; // 0: start, 1: end
        std::vector<Link> link_to_; // Nodes that this node connects to
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
    std::vector<Chain> GetChains();

    Node* GetPairNode(Node* n) {
        assert(n != nullptr);
        return nodes_.at(GetPairNodeId(n->Id())); // flip last bit to get the pair node
    }
    const Node* GetPairNode(const Node* n) const {
        assert(n != nullptr);
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
};

}