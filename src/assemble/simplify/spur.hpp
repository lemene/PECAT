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
    StringGraph& graph_;   
    int max_nodesize = 1;
    int max_length = -1;
};



}