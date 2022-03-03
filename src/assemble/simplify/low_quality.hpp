#pragma once

#include "./simplifier.hpp"
namespace fsa {


class LowQualitySimplifier : public Simplifier {
public:
    LowQualitySimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "low-quailty";
        desc_ = "Remove low-quality edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual bool PreCondition() { return !graph_.Options().read_file.empty(); }
    virtual void Running();
    virtual void Clear() { nodes_.clear(); quals_.clear(); }

    double ComputeThreshold();
    void RemoveLowQuality(double threshold);
    void RepairRmoved(double threshold, std::unordered_set<BaseEdge*> &removed);
    void ReactiveEdges(double threshold);
    StringGraph& graph_;   
    std::vector<BaseNode*> nodes_;
    std::vector<std::vector<double>> quals_;
};


}