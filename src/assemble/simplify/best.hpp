#pragma once

#include "./simplifier.hpp"
namespace fsa {

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
        std::array<int,3> score = { e->Score(), 0, 0 };     
        if (rvs_ != nullptr) {
            auto mm = rvs_->Test(*(e->ol_), false);
            score[1] = mm[0];
            score[2] = mm[1];
        } 
        return score;
    }

    bool EdgeScoreGreater(const EdgeScore &a, const EdgeScore &b) {
        // if ((a[1] + a[2]) < threshold_count || (b[1] + b[2]) < threshold_count) {
        //     return a[1] - a[2] > b[1] - b[2];
        // } else {
        //     double ra = (a[1] - a[2])*1.0 / (a[1] + a[2]) ;
        //     double rb = (b[1] - b[2])*1.0 / (b[1] + b[2]) ;
        //     return ra > rb;
        // }
        if (a[1]*snp_weight_ - a[2] > b[1]*snp_weight_ - b[2]) {
            return true;
        } else if (a[1]*snp_weight_ - a[2] == b[1] * snp_weight_ - b[2]) {
            return a[0] > b[0];
        } else {
            return false;
        }
    }

    bool EdgeScoreSignificantlyGreater(const EdgeScore &a, const EdgeScore &b) {

        // if (((a[1] - a[2]) - (b[1] - b[2]))  / 2>= threshold_count) {
        //     double ra = (a[1] + a[2]) > 0 ? (a[1] - a[2])*1.0 / (a[1] + a[2]) : 0.0;
        //     double rb = (b[1] + b[2]) > 0 ? (b[1] - b[2])*1.0 / (b[1] + b[2]) : 0.0;
        //     if ((ra - rb)/2 >= threshold_rate) {
        //         return true;
        //     }
        // }
        // return false;
        // auto get_score = [this](const EdgeScore &a) {
        //     return a[1]*snp_weight_ + a[2] == 0 ? 0.0 : (a[1]*snp_weight_ - a[2]) / (a[1]*snp_weight_+a[2]);
        // };
        // double sa = get_score(a);
        // double sb = get_score(b);
        // return sa - sb >= 0.33;

        auto ascore = a[1] * snp_weight_ - a[2];
        auto bscore = b[1] * snp_weight_ - b[2];

        if (bscore >= 0) {
            return ascore - bscore >= std::max<double>(th_count, std::abs(bscore)*th_rate[0]);
        } else {
            return ascore - bscore >= std::max<double>(th_count, std::abs(bscore)*th_rate[1]);
        }

    }

    auto FindBestOutEdge(BaseNode* n) -> std::vector<BestItem>;
    auto FindBestInEdge(BaseNode* n) -> std::vector<BestItem>;
    template<typename GN, typename GE>
    auto FindBestEdge(BaseNode* n, GN gn, GE ge) -> std::vector<BestItem>;

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