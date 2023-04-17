#pragma once

#include <vector>
#include <unordered_set>
#include "../../sequence.hpp"
#include "../read_variants.hpp"

#include "./simplifier.hpp"

namespace fsa {


class PhaseCrossSimplifier : public Simplifier {
public:
    PhaseCrossSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        name_ = "phase";
        desc_ = "Split cross points";
    }
 
    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual std::string GetParameters() const ;
    virtual void Running();
    virtual bool PreCondition() { rvs_ = graph_.GetAsmData().GetReadVariants(); return rvs_ != nullptr; }

    std::vector<std::vector<BaseNode*>> CollectCross();
    void ReplaceCross(std::vector<BaseNode*>& cand, std::vector<class PhasePath>& paths);
    StringGraph& graph_;   
    ReadVariants* rvs_ { nullptr };

    struct Options {

        int min_snp_count = 2;
        double min_snp_rate = 0.66;
        int min_read_support = 4;  // READ
        int min_valid_count = 4;
        double min_valid_rate = 0.66;
    };
    Options opts_;
    const int max_cand_size = 30;
};




class AsmDataset;
class ReadVariants;
class BaseNode;
class Overlap;

struct PhasePath {
    PhasePath(AsmDataset &ad, Seq::EndId start, ReadVariants *r, Seq::Id alt, const PhaseCrossSimplifier::Options& opts) ;
    void Extend(const std::vector<BaseNode*> cand, const std::vector<BaseNode*>& ends) ;

    void ExtendWithOtherPath(const std::vector<BaseNode*> cand, const std::vector<BaseNode*>& ends, PhasePath& other) ;

    const Overlap* FindBestExtendOverlap(Seq::EndId start, std::unordered_set<const Overlap*> &extols) ;
    bool Check(const Overlap* best,  std::unordered_set<const Overlap*> &extols) ;

    void MoveToNextTip(Seq::EndId tip, const Overlap* best,  std::unordered_set<const Overlap*> &extols);
    bool ReachEnds(const std::vector<BaseNode*> &ends) ;

    static bool IsIndependentPath(const std::vector<PhasePath> &path);

    void AddVariants(const std::vector<ReadVariants::Variants>* vars);
    void PrintVariants() const;
    std::unordered_map<int, std::unordered_map<int, int>> GetComfirmVariants() const;
    std::array<int,2> TestVariants(const std::unordered_map<int, std::unordered_map<int, int>>& avars, const std::vector<ReadVariants::Variants>* bvars) const;
    std::array<int,2> TestVariants(const std::unordered_map<int, std::unordered_map<int, int>>& avars, int name) const;

    void Debug(const char* const format, ...) const ;
    
    std::vector<Seq::EndId> tips;
    std::unordered_set<int> reads;
    std::unordered_set<Seq::EndId> visited;

    std::vector<Seq::EndId> dst;    // 
    AsmDataset &asmdata;
    ReadVariants *rvs;

    std::unordered_map<int, std::unordered_map<int, std::vector<std::array<int, 2>>>> variants;

    const PhaseCrossSimplifier::Options& opts;
};

class CrossPhaser {
public:
    CrossPhaser(PhaseCrossSimplifier& owner, const std::vector<BaseNode*> cand);
    bool Phase();

    void GetSnps(const BaseNode* node, const BaseNode* altnode);


    void Debug(const char* const format, ...) const ;
    std::vector<PhasePath> paths;
    std::vector<BaseNode*> ends;
protected:
    PhaseCrossSimplifier& owner_;
    std::vector<BaseNode*> cand_;
};

}