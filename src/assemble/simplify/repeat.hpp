#pragma once

#include "./simplifier.hpp"
namespace fsa {



class RepeatSimplifier : public Simplifier {
public:
    RepeatSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "repeat";
        desc_ = "Remove repeat edges";
    }

    virtual void Running();
    StringGraph& graph_;   
};



}