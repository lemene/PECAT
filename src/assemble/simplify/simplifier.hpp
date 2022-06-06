#pragma once

#include "../string_graph.hpp"

namespace fsa {


class Simplifier {
public:
    Simplifier(SgGraph& sg) : ori_graph_(sg) {}
    virtual ~Simplifier();
    void Simplify(const std::vector<std::string> &params=std::vector<std::string>(), bool required=false) {

        if (ParseParameters(params)) {
            if (PreCondition()) {
                LOG(INFO)("Simplifier(%s): %s, param=(%s)", name_.c_str(), desc_.c_str(), GetParameters().c_str());
                Running();
                LOG(INFO)("Simplifier(%s) end", name_.c_str());
            } else {
                if (!required) {
                    LOG(INFO)("Preconditions of simplifier(%s) are not satisfied", name_.c_str());
                } else {
                    LOG(ERROR)("Preconditions of simplifier(%s) are not satisfied", name_.c_str());
                }
            }
        } else {
            LOG(ERROR)("An error was encountered while parsing parameters:%s", JoinStrings(params, ":").c_str());
        }
        Clear();
    }

    virtual bool ParseParameters(const std::vector<std::string> &params) { return true; }
    bool ParseParameter(const std::vector<std::string> &params);
    virtual std::string GetParameters() const { return ""; }
    virtual bool PreCondition() { return true; }
    virtual void Running() = 0;
    virtual void Clear() {}

protected:
    void Debug(const char* const format, ...);
    std::string ToString(const SgEdge* e) { return e->Id().ToString(ori_graph_.GetAsmData().GetStringPool()); }
    std::string ToString(const SgNode* n) { return n->Id().ToString(ori_graph_.GetAsmData().GetStringPool()); }
protected:
    std::string name_;
    std::string desc_;
    SgGraph& ori_graph_;
};


// TODO 未整理的而类型和函数
using EdgeScore = std::array<int, 3>;   // alignment length, number of identical SNPs, number of different SNPs




class DuplicateSimplifier : public Simplifier {
public:
    DuplicateSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "duplicate";
        desc_ = "Remove duplicate edges";
    }

    virtual void Running();
    PathGraph& graph_;   
};

EdgeScore GetEdgeScore(const BaseEdge* e, ReadVariants* rvs) ;
bool EdgeScoreSignificantlyGreater(const EdgeScore &a, const EdgeScore &b, int count=4, int rate=0.3) ;
bool IsOutEdgeInconsistent(const BaseEdge* e0, const BaseEdge* e1, size_t N, const PhaseInfoFile* pif);

}   // namespace fsa