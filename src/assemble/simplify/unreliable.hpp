#pragma once

#include "./simplifier.hpp"
namespace fsa {


class UnreliableSimplifier : public Simplifier {
public:
    UnreliableSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "unreliable";
        desc_ = "Remove unreliable edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual void Running();
    virtual bool PreCondition() { return !graph_.Options().skip_purge; }
    StringGraph& graph_;   

    double max_cov_rate_ { 4.0 };
    size_t min_length_ { 10000 };
    double min_length_rate_ { 0.33 };
    size_t max_sub_length_ { 500 };
};



}