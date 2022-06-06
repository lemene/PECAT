#include "simplifier.hpp"


namespace fsa {

EdgeScore GetEdgeScore(const BaseEdge* e, ReadVariants* rvs) {
    std::array<int, 3> score = { e->Score(), 0, 0 };     //  length, XIANGTONG BUTONG SNP
    if (rvs != nullptr) {
        auto mm = rvs->Test(*e->ol_);
        score[1] = mm[0]; 
        score[2] = mm[1];
    } 
    return score;
};


bool EdgeScoreSignificantlyGreater(const EdgeScore &a, const EdgeScore &b, int count, int rate) {

    // if (((a[1] - a[2]) - (b[1] - b[2]))  / 2>= threshold_count) {
    //     double ra = (a[1] + a[2]) > 0 ? (a[1] - a[2])*1.0 / (a[1] + a[2]) : 0.0;
    //     double rb = (b[1] + b[2]) > 0 ? (b[1] - b[2])*1.0 / (b[1] + b[2]) : 0.0;
    //     if ((ra - rb)/2 >= threshold_rate) {
    //         return true;
    //     }
    // }
    return (a[1] - a[2]) - (b[1] - b[2]) >= std::max<int>((a[1] + a[2]) * 0.5 + (b[1]+b[2]) * 0.5, 6);
    if (((a[1] - a[2]) - (b[1] - b[2]))  / 2>= count) {
        double ra = (a[1] + a[2]) > 0 ? (a[1] - a[2])*1.0 / (a[1] + a[2]) : 0.0;
        double rb = (b[1] + b[2]) > 0 ? (b[1] - b[2])*1.0 / (b[1] + b[2]) : 0.0;
        if ((ra - rb)/2 >= rate) {
            return true;
        }
    }

    return false;
}

Simplifier::~Simplifier() {
}


bool Simplifier::ParseParameter(const std::vector<std::string> &param) {
    assert(param.size() >= 1);
    return false;
}


void Simplifier::Debug(const char* const format, ...) {
    if (DUMPER.IsWorking()) {
        va_list arglist;
        va_start(arglist, format);
        DUMPER[name_].Write(format, arglist);
        va_end(arglist);
    }
}



/** 如果两个节点的有多条长度（子节点数）<=3的的路径，则只保留第一个（排序根据节点名称） */
void DuplicateSimplifier::Running() {

    std::unordered_map<BaseEdge::ID, std::vector<PathEdge*>, BaseEdge::ID::Hash> dup_edges;
    std::vector<SgEdge*> cands = graph_.CollectEdges([&dup_edges](SgEdge* e) {
        if (e->IsType("simple")) {
            SimplePathEdge *se = static_cast<SimplePathEdge*>(e);
            if (se->path_.size() < 3) {
                auto id = BaseEdge::ID(0, e->InNode()->Id(), e->OutNode()->Id());
                auto iter = dup_edges.find(id);
                if (iter != dup_edges.end()) {
                    iter->second.push_back(se);
                }
                else {
                    dup_edges[id] = std::vector<PathEdge*>{ se };
                }
                
            }
        }
        return false;
    });

    std::unordered_set<PathEdge*> done;     // 判断对偶路径是否已经处理完成
    for (auto &i : dup_edges) {
        std::vector<PathEdge*> &dups = i.second;
        if (dups.size() > 1 && done.find(dups[0]) == done.end()) {
            done.insert(dups[0]);
            done.insert(graph_.ReverseEdge(dups[0]));
            for (size_t i = 1; i < dups.size(); ++i) {
                auto e = dups[i];
                e->Reduce("simple_dup", true);

                auto re = graph_.ReverseEdge(dups[i]);
                re->Reduce("simple_dup", true);
                
                done.insert(re);
                done.insert(e);
            }
        }
    }
}



bool IsOutEdgeInconsistent(const BaseEdge* e0, const BaseEdge* e1, size_t N, const PhaseInfoFile* pif) {
    
    auto get_nodes = [N](const BaseEdge* e) {
        std::vector<const BaseNode*> nodes;
        const BaseEdge* curr = e;
        for (size_t i = 0; i < N && curr != nullptr; ++i) {
            nodes.push_back(curr->OutNode());
            curr = nodes.back()->OutDegree() > 0 ? (const BaseEdge*)nodes.back()->OutEdge(0) : nullptr;
        }
        return nodes;
    };

    if (pif != nullptr) {
        std::vector<const BaseNode*> nodes0 = get_nodes(e0);
        std::vector<const BaseNode*> nodes1 = get_nodes(e1);

        for (auto n0 : nodes0) {
            for (auto n1 : nodes1) {
                if (pif->Contain(n0->ReadId(), n1->ReadId())) {
                    return true;
                }
            }
        }
    }

    return false;
}


} // namespace fsa