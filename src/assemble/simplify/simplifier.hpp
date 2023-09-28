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
            LOG(ERROR)("An error was encountered while parsing parameters: %s", JoinStrings(params, ":").c_str());
        }
        Clear();
    }

    virtual bool ParseParameters(const std::vector<std::string> &params) { return true; }
    bool ParseParameter(const std::vector<std::string> &params);
    virtual std::string GetParameters() const { return ""; }
    virtual bool PreCondition() { return true; }
    virtual void Running() = 0;
    virtual void Clear() {}

    // Common Function
    static std::list<SgEdge*> GetOutLinearPath(SgEdge* e);
    static std::list<const SgEdge*> GetOutLinearPath(const SgEdge* e);
    static std::list<SgEdge*> GetInLinearPath(SgEdge* e);
    static std::list<const SgEdge*> GetInLinearPath(const SgEdge* e);
    bool IsCross(SgNode* in_0) const ;
    bool TestLength(const std::list<SgEdge*> &path) const;
    bool TestExtends(SgNode* in0, SgNode* in1, SgNode* out0, SgNode* out1) const;
    bool TestCrossExtends(const std::list<const SgEdge*> in0_2_out0, const std::list<const SgEdge*> in0_2_out1, 
                          const std::list<const SgEdge*> in1_2_out0, const std::list<const SgEdge*> in1_2_out1);
    int PathCoreLength(const std::list<const SgEdge*> &path) const;

protected:
    const AsmOptions& Options() const { return ori_graph_.GetAsmData().opts_; }
    // short name
    std::string QueryStringById(Seq::Id id) {
        return ori_graph_.GetAsmData().GetStringPool().QueryStringById(id);
    }
    
    std::pair<const Overlap*, double> ReplaceHighQualityOverlap(const Overlap* ol, double threshold);
protected:
    void Debug(const char* const format, ...) const;
    void DebugEdges(const std::string& msg, const std::unordered_set<BaseEdge*>& edges);
    std::string ToString(const SgEdge* e) const { return e->Id().ToString(ori_graph_.GetAsmData().GetStringPool()); }
    std::string ToString(const SgNode* n) const { return n->Id().ToString(ori_graph_.GetAsmData().GetStringPool()); }
    std::string ToString(Seq::Id i) const { return ori_graph_.GetAsmData().GetStringPool().QueryStringById(i); }
protected:
    std::string name_;
    std::string desc_;
    SgGraph& ori_graph_;
};


class DuplicateSimplifier : public Simplifier {
public:
    DuplicateSimplifier(PathGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "duplicate";
        desc_ = "Remove duplicate edges";
    }

    virtual void Running();
    PathGraph& graph_;   
};

bool IsOutEdgeInconsistent(const BaseEdge* e0, const BaseEdge* e1, size_t N, const PhaseInfoFile* pif);
bool IsInEdgeInconsistent(const BaseEdge* e0, const BaseEdge* e1, size_t N, const PhaseInfoFile* pif);

}   // namespace fsa