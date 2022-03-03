#pragma once

#include "./simplifier.hpp"
namespace fsa {


class BubbleSimplifier : public Simplifier {
public:
    BubbleSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "Bubble";
        desc_ = "Detect bubble structures";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);

    virtual void Running();
    BubbleEdge*  Detect(PathNode* s, bool check, int depth_cutoff = 100, int width_cutoff = 1600);
    PathGraph& graph_;   
    bool simple { false };
};



}