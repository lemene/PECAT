#pragma once

#include "./simplifier.hpp"
namespace fsa {


class LoopSimplifier : public Simplifier {
public:
    LoopSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "loop";
        desc_ = "Detect loop structures";
    }

    struct NodeOrEdge{
        NodeOrEdge() {}
        NodeOrEdge(LoopNode* n) { node = n;}
        NodeOrEdge(LoopEdge* e) { edge = e;}
        LoopNode* node { nullptr };
        LoopEdge* edge { nullptr };
    } ;

    virtual void Running();
    void FindLoops();

    NodeOrEdge DetectLoop(SgNode* start);


    PathGraph& graph_;   
    size_t max_length = 5000000;
};




}