#pragma once

#include "./simplifier.hpp"
namespace fsa {


class ExtendSimplifier : public Simplifier {
public:
    ExtendSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "extend";
        desc_ = "extend dead notes";
    }

    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual bool PreCondition();
    virtual void Running();
    virtual void Clear() { }

protected:
    struct Path {
        Seq::EndId start;
        std::vector<std::pair<Seq::EndId, const Overlap*>> nodes;
        
        bool Empty() const { return nodes.empty(); }
    };

    void ReactiveContainedEdges();
    std::pair<Seq::EndId, const Overlap*> ExtendNode(Seq::EndId start, const std::unordered_set<Seq::Id>& cands, const std::unordered_set<Seq::Id>& visited);
    Path ExtendPath(BaseNode *node, const std::unordered_set<Seq::Id>& cands=std::unordered_set<Seq::Id>());
    std::vector<std::pair<Seq::EndId, const Overlap*>> ImproveLastEdge(const std::pair<Seq::EndId, const Overlap*>& ol);
    void ExtendLastEdge(const std::unordered_set<const Overlap*> exts, Seq::EndId s, Seq::EndId end, std::vector<std::pair<Seq::EndId, const Overlap*>>& path);
    void ModifyMainPath(const BaseNode* s, Path &main);
    Path GetSecondPath(const BaseNode* s, const Path &main);
    bool AddSecondPath(const Path& main, const Path& secondary, std::vector<Path>& paths);
    std::array<int, 2> CompareVariants(Seq::Id query, const std::unordered_set<Seq::Id>& target);
    Seq::Id GetReadId(const SgNode* n) const;
    void Repair(const std::unordered_set<const Overlap*> &ols);

    StringGraph& graph_;   

    size_t iterations { 3 };
    size_t max_number_of_nodes { 30 };
    size_t acceptable_aligned_length { 10000 };

};


}