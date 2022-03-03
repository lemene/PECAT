
#pragma once

#include "./simplifier.hpp"
namespace fsa {

class BridgeSimplifier : public Simplifier {
public:
    BridgeSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "bridge";
        desc_ = "Remove abnormal bridging edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual void Running();
    
    std::vector<SgEdge*> GetBridgePath(SgNode* node, size_t pri, size_t alt);
    void Repair(const std::unordered_set<const Overlap*> &ols);
    std::vector<const Overlap*> RepairIncompleteCross(const std::vector<SgEdge*>& path);
    std::vector<const Overlap*> GetAltCrossPath(Seq::EndId start, Seq::EndId end, const std::vector<SgEdge*>& path);
    bool IsAmbiguousPath(const std::vector<SgEdge*>& path);
    bool IsLinkedReversedNode(const std::vector<SgEdge*>& path, int max_depth) const;
  
    void DebugPath(const std::vector<SgEdge*> &path, const std::string& msg);

    StringGraph& graph_;   
    
    size_t max_length { 1000000 };
    size_t max_nodesize { 10 };
    bool check_reverse;
};


} // namespace fsa {