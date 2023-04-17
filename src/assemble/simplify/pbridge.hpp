
#pragma once

#include "./simplifier.hpp"
namespace fsa {

class PathBridgeSimplifier : public Simplifier {
public:
    PathBridgeSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "pbridge";
        desc_ = "Remove abnormal bridging edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual void Running();
    bool IsAbnormalBridge(const PathGraph::LinearPath &path);
    void CalcPathScore(PathGraph::LinearPath &path) const;
    bool TestOutExtend(PathEdge* e, int minlen, int minnode);
    bool TestInExtend(PathEdge* e, int minlen, int minnode);
    PathGraph& graph_;   

    size_t max_length { 1000000 };
    size_t max_nodesize { 20 };
    size_t round { 3 };
    size_t level {2};
};


} // namespace fsa {