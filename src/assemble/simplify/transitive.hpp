#pragma once

#include "./simplifier.hpp"
namespace fsa {


class TransitiveSimplifier : public Simplifier {
public:
    TransitiveSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "transitive";
        desc_ = "Remove transitive edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);

    virtual void Running();
    
    StringGraph& graph_;   
    size_t fuzz_ { 500 };
};


}