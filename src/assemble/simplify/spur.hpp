#pragma once

#include "./simplifier.hpp"
namespace fsa {

class SpurSimplifier : public Simplifier {
public:
    SpurSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "spur";
        desc_ = "Remove spur edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual std::string GetParameters() const;
    virtual void Running();
protected:
    bool IsTrivialEdge(const SgEdge* e, int min_nodesize, int min_length);
protected:
    StringGraph& graph_;   
    int min_nodesize {3};
    int min_length {-1};
};



}