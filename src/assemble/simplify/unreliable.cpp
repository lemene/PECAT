#include "unreliable.hpp"

namespace fsa {


bool UnreliableSimplifier::ParseParameters(const std::vector<std::string> &params) {
    assert(params[0] == "unreliable");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "l") {
            min_length_ = (size_t) std::stoul(it[1]);
        } else {
            if (!ParseParameter(it)) {
                return false;
            }
        }
    }
    return true;
}

std::string UnreliableSimplifier::GetParameters() const {
    std::ostringstream oss;
    oss << "l=" << min_length_;
    return oss.str();
}


void UnreliableSimplifier::Running() {
    std::unordered_set<BaseEdge*> removed;
    std::mutex mutex;

    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() > 1;
    });

    std::atomic<size_t> index { 0 };

    auto combine_func = [&removed, &mutex](std::unordered_set<BaseEdge*> &rs) {
        std::lock_guard<std::mutex> lock(mutex);
        removed.insert(rs.begin(), rs.end());
    };

    auto work_func = [this, combine_func, &nodes, &index](size_t tid) {
        std::unordered_set<BaseEdge*> rs;

        for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
            auto n = nodes[i];
            assert(n->OutDegree() > 1);

            Debug("cand node: %s\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
           
            const auto& rinfo = graph_.GetAsmData().GetReadInfo(n->ReadId());
            auto cov = rinfo.coverage[1];
            Debug("cov: %s, %zd\n", ToString(n).c_str(), cov);


            std::vector<BaseEdge*> edges(n->GetOutEdges());
            std::sort(edges.begin(), edges.end(), [](BaseEdge*a, BaseEdge*b) { return a->Score() > b->Score(); });


            auto get_min_cov = [this](const BaseEdge* e) {
                auto& ri = graph_.GetAsmData().GetReadInfo(e->OutNode()->ReadId());
                return *std::min_element(ri.coverage.begin(), ri.coverage.end());
            };

            auto get_cliff_out_position = [](const BaseNode* n, const ReadStatInfo& rinfo) {
                return n->Id().End() == 0 ? rinfo.cliff[0] : (rinfo.cliff[1] >= 0 ? rinfo.len - rinfo.cliff[1] : rinfo.cliff[1]);
            };
            auto get_cliff_in_position = [](const BaseNode* n, const ReadStatInfo& rinfo) {
                return n->Id().End() == 0 ?  (rinfo.cliff[1] >= 0 ? rinfo.len - rinfo.cliff[1] : rinfo.cliff[1]) : rinfo.cliff[0];
            };

            int cfdlen = get_cliff_out_position(n, rinfo) ;

            //std::unordered_set<const BaseEdge*> keeps {edges[0]};
            //int accu_cov = get_min_cov(edges[0]);
            //auto best_score = edges[0]->Score();
            std::unordered_set<const BaseEdge*> keeps ;
            int accu_cov = 0;;
            size_t best_score = 0;

            for (size_t ie = 0; ie < edges.size(); ++ie) {
                auto e = edges[ie];
                Debug("check edge: %s\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
                auto icov = get_min_cov(edges[ie]);
                auto iscore = edges[ie]->Score();
                auto & ri = graph_.GetAsmData().GetReadInfo(e->OutNode()->ReadId());
                bool remove = false;

                bool consist = true;
                for (auto ke : keeps) {
                    consist = !IsOutEdgeInconsistent(ke, edges[ie], 3, graph_.GetAsmData().GetInconsistentOverlaps());
                    if (consist) break;
                }

                if (consist) {
                    if (cfdlen >= 0) {
                        remove = best_score > (size_t) cfdlen + max_sub_length_ && iscore < (size_t) cfdlen + max_sub_length_;
                        Debug("cfdlen: %d %d, %zd %zd\n", max_sub_length_, cfdlen, iscore, best_score);

                    }
                    if (!remove) {
                        int cfdlen1 = get_cliff_in_position(edges[ie]->OutNode(), ri);
                        remove = iscore < (size_t) cfdlen1 + max_sub_length_;
                        Debug("cfdlen1: %d %d, %zd\n", max_sub_length_, cfdlen1, iscore);
                    }
                     
                    if (!remove) {
                        Debug("length: %d %d, %d, %f\n", iscore, best_score, min_length_, min_length_rate_);
                        if (iscore >= min_length_ || iscore >= best_score * min_length_rate_) {
                            remove = accu_cov + icov > cov * max_cov_rate_;
                        } else {
                            remove = true;
                        }
                    }   

                }
                if (!remove) {
                    best_score = std::max<int>(iscore, best_score);
                    keeps.insert(edges[ie]);
                    accu_cov += icov;
                } else {
                    Debug("remove: %s, cov(%zd, %zd, %zd, %f)\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(), accu_cov, icov, cov, (cov * max_cov_rate_));
                }

            }
            for (auto e : edges) {
                if (keeps.find(e) == keeps.end()) {
                    rs.insert(e);
                }
            }
        }
        combine_func(rs);
    };
    

    //MultiThreadRun((size_t)graph_.Options().thread_size, work_func);
    MultiThreadRun(1, work_func);
    Debug("Remove edges: %zd\n", removed.size());
    //graph_.ReduceEdges(removed, BaseEdge::RT_UNRELIABLE);
    for (auto r : removed) {
        if (!r->IsReduce()) {
    Debug("Remove edges: %lld, %s\n", r, ToString(r).c_str());
            r->Reduce(BaseEdge::RT_UNRELIABLE);
            graph_.ReverseEdge(r)->Reduce(BaseEdge::RT_UNRELIABLE);
    Debug("Remove edges done: %s\n", ToString(r).c_str());

        }
    }


}

// void UnreliableSimplifier::Running() {
//     const size_t sub_threshold = 500;
//     std::unordered_set<BaseEdge*> removed;
//     std::mutex mutex;
//     auto nodes = graph_.CollectNodes([](BaseNode* n) {
//         return n->OutDegree() > 1;
//     });

//     std::atomic<size_t> index { 0 };

//     auto combine_func = [&removed, &mutex](std::unordered_set<BaseEdge*> &rs) {
//         std::lock_guard<std::mutex> lock(mutex);

//         removed.insert(rs.begin(), rs.end());
//     };

//     auto work_func = [this, combine_func, sub_threshold, &nodes, &index](size_t tid) {
//         std::unordered_set<BaseEdge*> rs;

//         for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
//             auto n = nodes[i];
//             assert(n->OutDegree() > 1);
           
//             auto rid = n->ReadId();
//             auto& rinfo = graph_.GetAsmData().read_infos_.find(rid)->second;
//             auto cov = rinfo.coverage[1];

//             int cfdlen = n->Id().End() == 0 ? rinfo.cliff[0] : (rinfo.cliff[1] >= 0 ? rinfo.len - rinfo.cliff[1] : rinfo.cliff[1]) ;

//             std::vector<BaseEdge*> edges(n->GetOutEdges());
//             std::sort(edges.begin(), edges.end(), [](BaseEdge*a, BaseEdge*b) { return a->Score() > b->Score(); });

//             int accu = 0;
//             auto best_score = edges[0]->Score();
//             for (size_t ie = 1; ie < edges.size(); ++ie) {
//                 auto e = edges[ie];
//                 auto i = edges[ie]->OutNode()->ReadId();
//                 auto& ri = graph_.GetAsmData().read_infos_.find(i)->second;

//                 // span cliff
//                 if (cfdlen >= 0 && best_score > cfdlen + (int)sub_threshold && e->Score() < cfdlen + (int)sub_threshold) {
//                     rs.insert(e);
//                     if (graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()) == "1050333" || 
//                         graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()) == "1050333") {
//                             printf("cfdlen %d %d %d\n", cfdlen, best_score,e->Score());
//                         }
//                     continue;
//                 }

//                 if ((e->Score() < 10000 && e->Score()  < best_score/4) || (e->Score() < best_score/2 && cov*2 < accu + ri.coverage[1])) {
//                     rs.insert(e);
//                     if (graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()) == "1050333" || 
//                         graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()) == "1050333") {
//                             printf("cfdlen2 %d %d %d, %d %d %d\n", cfdlen, best_score,e->Score(), cov, accu, ri.coverage[1]);
//                             printf("cfdlen2 %s %s\n", graph_.GetAsmData().QueryNameById(edges[0]->InNode()->ReadId()).c_str(), graph_.GetAsmData().QueryNameById(edges[0]->OutNode()->ReadId()).c_str());
//                         }
//                 } else {
//                     accu += ri.coverage[1];
//                 }
//             }
//         }
//         combine_func(rs);
//     };
    

//     MultiThreadRun((size_t)graph_.Options().thread_size, work_func);

//     graph_.ReduceEdges(removed, BaseEdge::RT_UNRELIABLE);

// }


} // namespace fsa