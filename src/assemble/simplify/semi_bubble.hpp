#pragma once

#include "./simplifier.hpp"
namespace fsa {

class SemiBubbleSimplifier : public Simplifier {
public:
    SemiBubbleSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "semi";
        desc_ = "Detect semi-bubble structures";
    }

    
    virtual bool PreCondition();

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual void Running();
    SemiBubbleEdge* Detect(PathNode* start_node, int depth_cutoff = 100, int width_cutoff = 1600, int length_cutoff = 50000000);
    PathGraph& graph_;   
};

}