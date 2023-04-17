#pragma once

#include "./simplifier.hpp"
namespace fsa {


class Spur2Simplifier : public Simplifier {
public:
    Spur2Simplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "spur2";
        desc_ = "Remove spur edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual void Running();

    void RemoveDeadEnds();
    bool TestInExtend(SgEdge* e, int minlen, int minnode);
    PathGraph& graph_;   

    size_t max_length { 1000000 };
    size_t max_nodesize { 20 };
};

}