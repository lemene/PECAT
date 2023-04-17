#pragma once

#include "./simplifier.hpp"
namespace fsa {


class PTransitiveSimplifier : public Simplifier {
public:
    PTransitiveSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "ptransitive";
        desc_ = "Remove transitive edges";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);

    virtual void Running();

    std::unordered_set<const SgNode*> CollectImportantNodes(const SgEdge* e, std::unordered_set<const SgNode*>& middle) const;
    bool HasImportantInEdges(const SgNode* n) const;
    bool HasImportantOutEdges(const SgNode* n) const;
    std::list<const SgEdge*> GetPath(const SgEdge*) const;
    std::unordered_set<const SgEdge*> GetExtendPath(const SgEdge*, const std::unordered_set<SgEdge*>& rm0, const std::unordered_set<SgEdge*>& rm2) const;
    bool IsTrivialPath(const std::list<const SgEdge*> &path) const;
    bool IsSubOf(const std::unordered_set<const SgNode*>& a, const std::unordered_set<const SgNode*>& b) const;
    int MaxCoreLengthOutEdges(const SgNode* n) const;
    void RemoveSmallLoops();
protected:
    PathGraph& graph_;   
    size_t max_length_ { 500000 };
    size_t max_nodesize_ { 10 };
};


}