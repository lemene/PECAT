#pragma once

#include "./simplifier.hpp"
namespace fsa {


// TODO 未整理的而类型和函数
using EdgeScore = std::array<int, 4>;   // alignment length, number of identical SNPs, number of different SNPs

class BestOverlapsSimplifier : public Simplifier {
public:

    struct BestItem {
        BaseEdge* e;
        EdgeScore s;
    };
     
    BestOverlapsSimplifier(StringGraph& graph) : Simplifier(graph), graph_(graph) {
        
        name_ = "best";
        desc_ = "Retain the best overlaps";
    }
    
    virtual bool ParseParameters(const std::vector<std::string> &params);
    virtual std::string GetParameters() const;
    virtual bool PreCondition() { 
        rvs_ = graph_.GetAsmData().GetReadVariants();
        bad_ols_ = graph_.GetAsmData().GetInconsistentOverlaps();
        return true; 
    }
    
    virtual void Running();
    virtual void Clear() {
        cands_.clear();
        best_edges.clear();
    }


    EdgeScore GetEdgeScore(const BaseEdge *e) {
        EdgeScore score = { (int)e->Score(), 0, 0, (int)(e->Identity()*1000000) };     
        if (rvs_ != nullptr) {
            auto mm = rvs_->GetClosestSnps(*(e->ol_));
            score[1] = mm[0];
            score[2] = mm[1];
        } 
        return score;
    }

    bool EdgeScoreGreater(const EdgeScore &a, const EdgeScore &b) {
        if (a[1]*snp_weight_ - a[2] > b[1]*snp_weight_ - b[2]) {
            return true;
        } else if (a[1]*snp_weight_ - a[2] == b[1] * snp_weight_ - b[2]) {
            if (a[0] > b[0]*2) {
                return true;
            } else if (a[0]*2 < b[0]) {
                return false;
            } else {
                // TODO
                double err_a = 100 - a[3]*1.0 / 1000000;
                double err_b = 100 - b[3]*1.0 / 1000000;
                if (err_a * 2 < err_b) {
                    return true;
                } else if (err_a > err_b*2) {
                    return false;
                } else {
                    return a[0] * (1-err_a) > b[0] * (1-err_b);
                }
            }
            return a[0] > b[0];
        } else {
            return false;
        }
        // auto rate_score = [](const EdgeScore &a) {
        //     auto s = a[1] + a[2] == 0;
        //     return s > 0 ? 1.0*a[2] / s : 0.0;
        // };
        // auto arscore = rate_score(a);
        // auto brscore = rate_score(b);
        // if (a[2] < b[2]) {
        //     return true;
        // } else if (a[2] == b[2] && arscore < brscore) {
        //     return true;
        // } else if (a[2] == b[2] && arscore == brscore) {
        //     return a[0] > b[0];
        // } else {
        //     return false;
        // }
    }

    bool EdgeScoreSignificantlyGreater(const EdgeScore &a, const EdgeScore &b) {
        auto ascore = a[1] * snp_weight_ - a[2];
        auto bscore = b[1] * snp_weight_ - b[2];

        if (bscore >= 0) {
            return ascore - bscore >= std::max<double>(th_count, std::abs(bscore)*th_rate[0]);
        } else {
            return ascore - bscore >= std::max<double>(th_count, std::abs(bscore)*th_rate[1]);
        }

        // auto rate_score = [](const EdgeScore &a) {
        //     auto s = a[1] + a[2];
        //     return s > 0 ? 1.0*a[2] / s : 0.0;
        // };
        // size_t count = th_count;
        // double rate = th_rate[0];
        // auto arscore = rate_score(a);
        // auto brscore = rate_score(b);
        // return a[2] + count < b[2] && arscore * rate < brscore;
    }
    
    auto FindBestOutEdge(BaseNode* n) -> std::vector<BestItem>;

    void FindCandidateBestEdges();
    void ComfirmCandidateBestEdges();
protected:
    StringGraph& graph_;    
    int th_count { 4 };
    std::array<double,2> th_rate { {2.0, 0.66} };

    ReadVariants* rvs_ { nullptr };
    PhaseInfoFile* bad_ols_ { nullptr };
    double snp_weight_ { 0.5 };
    
    std::vector<std::array<std::vector<BestItem>,2>> cands_;
    std::unordered_set<BaseEdge*> best_edges;
};



}