#pragma once

#include "./simplifier.hpp"
namespace fsa {


class CrossSimplifier : public Simplifier {
public:
    CrossSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "cross";
        desc_ = "Detect cross structures";
    }

    virtual bool PreCondition();
    virtual void Running();
    
    bool IsInconsistent(const SgNode* n0, const SgNode* n1) const;
    

    CrossNode* DetectCross(SgNode* n) const;

    bool TrySoluteCrossNode(CrossNode* node, std::unordered_set<SgEdge*>& removed);
    bool IsLinkedReversedNode(const SgNode* front, const SgNode* back, size_t max_depth);
protected:
    PathGraph& graph_;   
    size_t max_node_size_ { 10 };
    size_t max_length_ { 1000000 };
};

}