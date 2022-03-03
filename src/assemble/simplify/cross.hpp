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
    
    bool IsCross(const SgNode* n) const;
    bool IsInconsistent(const SgNode* n0, const SgNode* n1) const;
    bool TestLength(const SgNode* in0, const SgNode* in1) const;
    bool TestExtends(const SgNode* in0, const SgNode* in1, const SgNode* out0, const SgNode* out1) const;
    CrossNode* DetectCross(SgNode* n) const;
    PathGraph& graph_;   
    size_t max_node_size { 5 };
};

}