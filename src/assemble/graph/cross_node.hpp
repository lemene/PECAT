#pragma once

#include <list>

#include "../string_node.hpp"

namespace fsa {

class CrossNode : public PathNode {
public:
    CrossNode(const std::vector<SgNode*>& in_nodes, const std::vector<SgNode*>& out_nodes);
    CrossNode* Reverse(PathGraph& sg) const;
    bool IsRaw() const;

    size_t OriginInDegree() const { return in_nodes_[0]->InDegree() + in_nodes_[1]->InDegree(); }
    size_t OriginOutDegree() const { return  out_nodes_[0]->OutDegree() + out_nodes_[1]->OutDegree(); }

    const SgEdge* OriginOutEdge(size_t i) const { 
        if (i < out_nodes_[0]->OutDegree()) {
            return out_nodes_[0]->OutEdge(i);
        } else {
            return out_nodes_[1]->OutEdge(i-out_nodes_[0]->OutDegree());
        }
    }
    const SgEdge* OriginInEdge(size_t i) const {
        if (i < in_nodes_[0]->InDegree()) {
            return in_nodes_[0]->InEdge(i);
        } else {
            return in_nodes_[1]->InEdge(i-in_nodes_[0]->InDegree());
        }
    }
    virtual bool IsConsistOf(const SgNode* n) const { 
        assert(in_nodes_.size() == 2 && out_nodes_.size() == 2);

        return in_nodes_[0] == n || in_nodes_[1] == n || out_nodes_[0] == n || out_nodes_[1] == n;
    }

    const std::list<SgEdge*>& GetInternalPath(size_t i_in, size_t i_out) {
        assert(origin_edges_.size() == 4 && i_in < 2 && i_out < 2);
        return origin_edges_[2*i_in + i_out].edges;
    }

    void ReduceSub();
    static std::list<SgEdge*> GetOutLinearPath(SgEdge* e);
    static std::list<SgEdge*> GetInLinearPath(SgEdge* e);

protected:
    struct LinearPath {
        std::list<SgEdge*> edges;
        SgNode* OutNode() const;
        SgNode* InNode()  const;
    };
protected:
    SgNode::ID CreateID() const;

    std::vector<SgEdge*> origin_in_edges_;
    std::vector<SgEdge*> origin_out_edges_;
    std::vector<LinearPath> origin_edges_;
    std::vector<SgNode*> in_nodes_;
    std::vector<SgNode*> out_nodes_;
};

}

