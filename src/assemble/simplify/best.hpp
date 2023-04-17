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
        std::array<int,3> score = { (int)e->Score(), 0, 0 };     
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
    bool EdgeScoreSignificantlyGreater2(const BestItem &a, const BestItem &b) {
 
        size_t count = th_count;
        double rate = th_rate[0];
        
        auto rate_score = [](const EdgeScore &a) {
            auto s = a[1] + a[2];
            return s > 0 ? 1.0*a[2] / s : 0.0;
        };

        if (a.e->OutNode()->InDegree() == 1 && b.e->OutNode()->InDegree() >= 2 && rvs_ != nullptr) {
            if (a.s[1] - a.s[2] >= (int)count && (a.s[1] - a.s[2]) - (b.s[1] - b.s[2]) >= (int)count) {
                auto a_in_node = a.e->InNode();
                auto arscore = rate_score(a.s);

                for (size_t i = 0; i < b.e->OutNode()->InDegree(); ++i) {
                    if (b.e->OutNode()->InEdge(i) != b.e) {
                        auto n = b.e->OutNode()->InNode<BaseNode>(i);
                        assert(rvs_ != nullptr);
                        auto mm = rvs_->Test(n->ReadId(), a_in_node->ReadId());
                        auto bs = b.s;
                        bs[1] += mm[0];
                        bs[2] += mm[1];
                        auto brscore = rate_score(bs);
                        
                        if (a.s[2] + (int)count < bs[2] && arscore * rate < brscore) {
                            return true;
                        }
                        
                    }
                }
            }
        }
        return false;
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