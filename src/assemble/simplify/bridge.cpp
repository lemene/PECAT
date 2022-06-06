#include "bridge.hpp"

namespace fsa {



bool BridgeSimplifier::ParseParameters(const std::vector<std::string> &params) { 
    //assert(params[0] == "bridge");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "length") {
            max_length = (size_t) std::stoul(it[1]);
        } else if (it[0] == "nodesize") {
            max_nodesize = (size_t) std::stoul(it[1]);
        } else {
            return false;
        }
    }
    return true;
}


// void BridgeSimplifier::Running() {
//     std::unordered_set<BaseEdge*> removed;

//     auto edges = graph_.CollectEdges([](BaseEdge* e){
//         return !e->IsReduce();
//     });

//     for (auto &e : edges) {
//         assert(!e->IsReduce());
  
//         if (e->InNode()->InDegree() == 1 && e->InNode()->OutDegree() == 2 &&
//             e->OutNode()->InDegree() == 2 && e->OutNode()->OutDegree() == 1) {
            
//             int minlen = (e->Length() + graph_.ReverseEdge(e)->Length()) / 2 * 3;      // TODO how to determine the threshold
//             int minnode = 3;

//             std::array<bool, 2> out_node_has_long_edge = {false, false};
//             for (auto ie : e->OutNode()->GetInEdges()) {
//                 if (ie == e) continue;

//                 out_node_has_long_edge[0] = ie->TestInExtend(minlen, minnode);
//                 out_node_has_long_edge[1] = e->OutNode()->GetOutEdges()[0]->TestOutExtend(minlen, minnode);
//                 break;
//             }

//             std::array<bool, 2> in_node_has_long_edge = {false, false};
//             for (auto ie : e->InNode()->GetOutEdges()) {
//                 if (ie == e) continue;

//                 in_node_has_long_edge[0] = e->InNode()->GetInEdges()[0]->TestInExtend(minlen, minnode);
//                 in_node_has_long_edge[1] = ie->TestOutExtend(minlen, minnode);
//                 break;

//             }

//             if ((out_node_has_long_edge[0] && out_node_has_long_edge[1]) || (in_node_has_long_edge[0] && in_node_has_long_edge[1])) {
//                 removed.insert(e);
//             }
            
            
//         }
//     }

//     LOG(INFO) ("Remove abnormal bridging edges: %d", removed.size()*2);
//     graph_.ReduceEdges(removed, BaseEdge::RT_BRIDGED);
// }

std::vector<SgEdge*> BridgeSimplifier::GetBridgePath(SgNode* node, size_t pri, size_t alt) {
    assert(node->InDegree() == 1 && node->OutDegree() == 2);
    assert((pri == 0 || pri == 1) && (alt == 0 || alt == 1) && pri != alt);
  
    std::vector<SgEdge*> path = graph_.GetLinearPath(node->OutEdge(pri), max_length, max_nodesize);
    if (path.size() > 0 && path.back()->OutNode()->InDegree() == 2 && path.back()->OutNode()->OutDegree() == 1) {
        auto altnode_list = graph_.GetEgoNodes(node->OutEdge(alt)->OutNode(), path.size() * 10);
        std::unordered_set<SgNode*> altnodes(altnode_list.begin(), altnode_list.end());

        bool found = false;
        for (auto e : path) {
            auto n = e->OutNode();
            //if (altnodes.find(n) != altnodes.end() || altnodes.find(graph_.ReverseNode((BaseNode*)n)) != altnodes.end()) {
            if (altnodes.find(n) != altnodes.end() ) {
                found = true;
                Debug("GetEgoNodes found:  %s\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());

                break;
            }
        }

        if (!found) {
            return path;
        }
    }

    return std::vector<SgEdge*>();
}

std::vector<const Overlap*> BridgeSimplifier::RepairIncompleteCross(const std::vector<SgEdge*>& path) {
    assert(path.front()->InNode()->OutDegree() == 2);

    auto pif = graph_.GetAsmData().GetInconsistentOverlaps();
    if (pif != nullptr) {
        auto on0 = path.front()->InNode()->OutEdge(0)->OutNode<BaseNode>()->ReadId();
        auto on1 = path.front()->InNode()->OutEdge(1)->OutNode<BaseNode>()->ReadId();

        auto in0 = path.back()->OutNode()->InEdge(0)->InNode<BaseNode>()->ReadId();
        auto in1 = path.back()->OutNode()->InEdge(1)->InNode<BaseNode>()->ReadId();
        
        bool cand_cross = pif->Contain(on0, on1) && pif->Contain(in0, in1);

        if (cand_cross) {
            // 
            BaseNode* out_node = nullptr;
            for (size_t i = 0 ; i < path.front()->InNode()->OutDegree(); ++i) {
                if (path.front()->InNode()->OutEdge(i) != path.front()) {
                    out_node = path.front()->InNode()->OutEdge(i)->OutNode<BaseNode>();
                    break;
                }
            }
        
            BaseNode* in_node = nullptr;
            for (size_t i = 0 ; i < path.back()->OutNode()->InDegree(); ++i) {
                if (path.back()->OutNode()->InEdge(i) != path.back()) {
                    in_node = path.back()->OutNode()->InEdge(i)->InNode<BaseNode>();
                    break;
                }
            }

            if (out_node->InDegree() == 1 && out_node->OutDegree() == 1 && in_node->InDegree() == 1 && in_node->OutDegree() == 1) {
                auto e = graph_.QueryEdge(BaseEdge::ID(0, in_node->Id(),  out_node->Id()));
                if (e == nullptr) {
                    return GetAltCrossPath(in_node->Id().MainNode(), out_node->Id().MainNode(), path);
                } else if (e->IsType("base") && static_cast<BaseEdge*>(e)->GetReduceType() == BaseEdge::RT_NO_BEST) {
                    return std::vector<const Overlap*>({nullptr});  // TODO 需要更换更换参数
                //} else if (pif->Contain(in_node->ReadId(), out_node->ReadId())) {
                //    return std::vector<const Overlap*>({nullptr});  // TODO 需要更换更换参数
                } else if (e->IsType("base") && static_cast<BaseEdge*>(e)->GetReduceType() == BaseEdge::RT_TRANSITIVE) {
                    return std::vector<const Overlap*>({static_cast<BaseEdge*>(e)->ol_});  // TODO 需要更换更换参数
                }
            }

        }
    }

    return std::vector<const Overlap*>();
}

void BridgeSimplifier::Repair(const std::unordered_set<const Overlap*> &ols) {

    auto get_edge = [this](const Overlap* ol) {
        if (ol->SameDirect()) {
            Seq::EndId fB = Seq::IdToEndId(ol->a_.id, 0);
            Seq::EndId gB = Seq::IdToEndId(ol->b_.id, 0);
            return graph_.GetEdge(fB, gB);
        } else {
            Seq::EndId fB = Seq::IdToEndId(ol->a_.id, 0);
	        Seq::EndId gE = Seq::IdToEndId(ol->b_.id, 1);
            return graph_.GetEdge(fB, gE);
        }
    };
    
    for (auto ol : ols) {

        auto e = get_edge(ol);
        if (e != nullptr) {
            graph_.ReactiveEdges(std::vector<BaseEdge*>({e}));

        } else {
            graph_.AddOverlap(ol);
        }

    }
}

std::vector<const Overlap*> BridgeSimplifier::GetAltCrossPath(Seq::EndId start, Seq::EndId end, const std::vector<SgEdge*>& path) {
    DebugPath(path, "cross");
    auto &asmdata = graph_.GetAsmData();

    auto get_extend = [&asmdata](Seq::EndId start) {

        Seq::Id sid = Seq::EndIdToId(start);
        int send = Seq::End(start);

        std::unordered_set<const Overlap*> ols = asmdata.GetExtendOverlaps(sid, send);
        std::unordered_map<Seq::EndId, const Overlap*> links;
        for (auto ol : ols) {
            auto& alt = ol->GetOtherRead(sid);
            auto rend = send == 0 ? (ol->SameDirect() ? 0 : 1) : (ol->SameDirect() ? 1 : 0);
            links[Seq::IdToEndId(alt.id, rend)] = ol;
        }
        return links;
    };

    auto is_contained_by_path = [&asmdata, this](Seq::Id id, const std::vector<SgEdge*>& path) {
        std::unordered_set<Seq::Id> reads ;
        reads.insert(path[0]->InNode<BaseNode>()->ReadId());
        for (size_t i = 0; i < path.size(); ++i) {
            reads.insert(path[i]->OutNode<BaseNode>()->ReadId());
        }
        
        bool contained = false;
        std::vector<Seq::Id> check { id };
        std::unordered_set<Seq::Id> done;   // avoid loop
        while (check.size() > 0) {
            auto cid = check.back();
            check.pop_back();
            done.insert(cid);

            auto rinfo = asmdata.GetReadInfo(cid);
            Debug("rinfo x: %s\n", asmdata.QueryNameById(rinfo.id).c_str());
            if (rinfo.filtered.type == RdReason::RS_CONTAINED) {
                auto oid = rinfo.filtered.sub[0];
                Debug("ccc %zd, %s %s\n", reads.size(), asmdata.QueryNameById(oid).c_str(),asmdata.QueryNameById(cid).c_str());
                if (reads.find(oid) != reads.end()) {
                Debug("ccc0\n");
                    contained = true;
                    break;
                } else {
                    if (done.find(oid) == done.end()) {
                        check.push_back(oid);
                    }
                }
            }
        }

        return contained;

    };

    auto check_overlap_quality = [this](const Overlap& ol) {
        auto qual = graph_.GetOverlapQuality(ol);
        Debug("qual: %f %f\n", qual, graph_.Options().min_identity);
        return qual >= graph_.Options().min_identity;
    };


    // 一级
    auto extend1 = get_extend(start);
    Debug("extend1 %s %zd\n", asmdata.QueryNameById(Seq::EndIdToId(start)).c_str(), extend1.size());
    for (auto ext: extend1) {
        Debug("  ext1: %s \n", asmdata.QueryNameById(Seq::EndIdToId(ext.first)).c_str());
    }
    if (extend1.find(end) == extend1.end()) {

        for (const auto &s : extend1) {
            auto rid = s.first;
            auto rinfo = asmdata.GetReadInfo(Seq::EndIdToId(rid));
            if (rinfo.filtered.type == RdReason::RS_CONTAINED) {
                Debug("rinfo: %s\n", asmdata.QueryNameById(rinfo.id).c_str());
                if (is_contained_by_path(rinfo.id, path)) {
                    auto extend2 = get_extend(rid);
                    auto e2 = extend2.find(end);
                    if (e2 != extend2.end()) {
                        Debug("mid0: %s\n", asmdata.QueryNameById(rinfo.id).c_str());
                        if (check_overlap_quality(*(s.second)) && check_overlap_quality(*(e2->second))) {
                            Debug("mid: %s\n", asmdata.QueryNameById(rinfo.id).c_str());

                            return std::vector<const Overlap*>({s.second, e2->second});
                        }

                    }
                }

            }

        }
    }


    return std::vector<const Overlap*>();
    

}

bool BridgeSimplifier::IsAmbiguousPath(const std::vector<SgEdge*>& path) {


    int minlen = graph_.PathLength(path) * 2;      // TODO how to determine the threshold
    int minnode = path.size() * 2;

    assert(path.front()->InNode()->InDegree() == 1 && path.front()->InNode()->OutDegree() == 2);
    std::array<bool, 2> in_node_has_long_edge = {false, false};
    for (size_t i = 0; i < path.front()->InNode()->OutDegree(); ++i) {
        auto e = path.front()->InNode()->OutEdge(i);
        if (e == path.front()) continue;

        in_node_has_long_edge[0] = static_cast<BaseEdge*>(path.front()->InNode()->InEdge(0))->TestInExtend(minlen, minnode,false);
        in_node_has_long_edge[1] = static_cast<BaseEdge*>(e)->TestOutExtend(minlen, minnode,false);
        break;
    }

    assert(path.back()->OutNode()->InDegree() == 2 && path.back()->OutNode()->OutDegree() == 1);
    std::array<bool, 2> out_node_has_long_edge = {false, false};
    for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
        auto e = path.back()->OutNode()->InEdge(i);
        if (e == path.back()) continue;

        out_node_has_long_edge[0] = static_cast<BaseEdge*>(path.back()->OutNode()->OutEdge(0))->TestOutExtend(minlen, minnode,false);
        out_node_has_long_edge[1] = static_cast<BaseEdge*>(e)->TestInExtend(minlen, minnode,false);
        break;
    }
    Debug("check (%d,%d): %d %d %d %d\n", minlen, minnode, out_node_has_long_edge[0], out_node_has_long_edge[1], in_node_has_long_edge[0], in_node_has_long_edge[1]);
    return out_node_has_long_edge[0] && out_node_has_long_edge[1] && in_node_has_long_edge[0] && in_node_has_long_edge[1];            
}


int BridgeSimplifier::TestConsistent(const std::vector<SgEdge*>& path) {
    int s = 0;
    
    BaseEdge* ie = static_cast<BaseEdge*>(path.front());
    
    Debug("TestConsistent IN %d", ie->subject_);
    if (ie->subject_) {
        s ++;
    } 
    
    for (size_t i = 0; i < ie->InNode()->OutDegree(); ++i) {
        auto iie = ie->InNode()->OutEdge<BaseEdge>(i);
        if (iie == ie) continue;
        Debug("TestConsistent IN other %d", iie->subject_);
        if (iie->subject_) {
            s --;
        }
    }

    BaseEdge* oe = static_cast<BaseEdge*>(path.back());
    Debug("TestConsistent Out %d", ie->subject_);
    if (oe->subject_) {
        s ++;
    } 
    for (size_t i = 0; i < oe->OutNode()->InDegree(); ++i) {
        auto ioe = oe->OutNode()->InEdge<BaseEdge>(i);
        Debug("TestConsistent Out other %d", ioe->subject_);
        if (ioe == oe) continue;
        if (ioe->subject_) {
            s --;
        }
    }
    
    return s;

  
}

bool BridgeSimplifier::IsLinkedReversedNode(const std::vector<SgEdge*>& path, int max_depth) const {
    auto rn = static_cast<SgGraph&>(graph_).ReverseNode(path.front()->InNode());
    
    auto curr = path.back()->OutNode();
    int depth = 0;
    while (depth < max_depth ) {
        if (curr == rn) return true;
        depth += 1;
        if (curr->OutDegree() == 1) {
            curr = curr->OutNode(0);
        } else {
            break;
        }
    }

    return false;
}

void BridgeSimplifier::DebugPath(const std::vector<SgEdge*> &path, const std::string& msg) {
    if (path.size() > 0) {
        Debug("%s: ", msg.c_str());

        
        Debug("%s", path[0]->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
        for (size_t i = 1; i < path.size(); ++i) {
            Debug(" -> %s", path[i]->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
        }
        Debug("\n");
    }
}

// void BridgeSimplifier::Running() {
//     graph_.SaveEdges(graph_.GetAsmData().OutputPath("test_edges.gz"));
//     std::unordered_set<BaseEdge*> removed;
//     std::unordered_set<const Overlap*> repaird;

//     auto cands = graph_.CollectNodes([](BaseNode* n) {
//         return n->InDegree() == 1 && n->OutDegree() == 2;
//     });

//     std::unordered_set<const SgEdge*> done;
//     for (auto &n : cands) {
//         assert(n->OutDegree() == 2);

//         Debug("cand node: %s\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());

//         std::vector<std::vector<SgEdge*>> paths {GetBridgePath(n, 0, 1), GetBridgePath(n, 1, 0)};

//         for (auto&path : paths) {

//             if (path.size() > 0) {
//                 DebugPath(path, "cands");
                
//                 if (done.find(path.front()) != done.end()) {
//                     continue;
//                 } else {
//                 }

//                 if (IsLinkedReversedNode(path, 20)) {
//                     Debug("found reversed: %zd\n", path.size());
//                     for (auto e : path) {
//                         removed.insert((BaseEdge*)e);
//                     }
//                     done.insert(path.front());
//                     done.insert(path.back());
//                     done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
//                     done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
//                     continue;
//                 }

//                 auto ols = RepairIncompleteCross(path);
//                 if (ols.size() > 0 && ols[0] != nullptr) {
//                     repaird.insert(ols.begin(), ols.end());
//                     Debug("found incomplete cross: %zd\n", ols.size());
//                     done.insert(path.front());
//                     done.insert(path.back());
//                     done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
//                     done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
                    
//                     continue;
//                 }

//                 if ((ols.size() == 1 && ols[0] == nullptr) || IsAmbiguousPath(path)) {
//                     Debug("found ambiguous path\n");
//                     for (auto e : path) {
//                         removed.insert((BaseEdge*)e);
//                     }
//                     done.insert(path.front());
//                     done.insert(path.back());
//                     done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
//                     done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
//                     continue;
//                 }
//             }
//         }

//         Debug("cand node ok: %s\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());

//     }
//     LOG(INFO) ("Remove abnormal bridging edges: %d", removed.size()*2);
//     graph_.ReduceEdges(removed, BaseEdge::RT_BRIDGED);
//     LOG(INFO) ("Repair cross structures: %d", repaird.size()*2);
//     Repair(repaird);
// }

void BridgeSimplifier::Running() {
    std::unordered_set<BaseEdge*> removed;
    std::unordered_set<const Overlap*> repaird;

    auto cands = graph_.CollectNodes([](BaseNode* n) {
        return n->InDegree() == 1 && n->OutDegree() == 2;
    });

    std::vector<std::vector<SgEdge*>> paths;
    for (auto &n : cands) {
        auto p0 = GetBridgePath(n, 0, 1);
        if (p0.size() > 0) {
            paths.push_back(p0);
        }
        auto p1 = GetBridgePath(n, 1, 0);
        if (p1.size() > 0) {
            paths.push_back(p1);
        }
    }

    std::sort(paths.begin(), paths.end(), [](const std::vector<SgEdge*>& a, const std::vector<SgEdge*>& b) {
        return a.size() < b.size();
    });

    std::unordered_set<const SgEdge*> done;
    for (auto &path : paths) {
        assert(path.size() > 0);
        DebugPath(path, "cands");
        
        if (done.find(path.front()) != done.end()) {
            continue;
        } else {
        }

        if (IsLinkedReversedNode(path, 20)) {
            Debug("found reversed: %zd\n", path.size());
            for (auto e : path) {
                removed.insert((BaseEdge*)e);
            }
            done.insert(path.front());
            done.insert(path.back());
            done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
            done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
            continue;
        }

        auto ols = RepairIncompleteCross(path);
        if (ols.size() > 0 && ols[0] != nullptr) {
            repaird.insert(ols.begin(), ols.end());
            Debug("found incomplete cross: %zd\n", ols.size());
            done.insert(path.front());
            done.insert(path.back());
            done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
            done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
            
            continue;
        }

        // 另外一条路径标记为no_best
        if ((ols.size() == 1 && ols[0] == nullptr) ) {
            Debug("found ambiguous path\n");
            for (auto e : path) {
                removed.insert((BaseEdge*)e);
            }
            done.insert(path.front());
            done.insert(path.back());
            done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
            done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
            continue;
        }

        if (IsAmbiguousPath(path)) {
            auto r = TestConsistent(path);
            if (r >= 2) {
                
                Debug("found Consistent path, %d\n", r);
                BaseEdge* ie = static_cast<BaseEdge*>(path.front());
                for (size_t i = 0; i < ie->InNode()->OutDegree(); ++i) {
                    auto iie = ie->InNode()->OutEdge<BaseEdge>(i);
                    if (iie != ie)  removed.insert(iie);
                }
                
                BaseEdge* oe = static_cast<BaseEdge*>(path.back());
                for (size_t i = 0; i < oe->OutNode()->InDegree(); ++i) {
                    auto ioe = oe->OutNode()->InEdge<BaseEdge>(i);
                    if (ioe != oe)  removed.insert(ioe);
                }
            } else if (r <= 0) {
                Debug("found ambiguous path\n");
                for (auto e : path) {
                    removed.insert((BaseEdge*)e);
                }
            }
            done.insert(path.front());
            done.insert(path.back());
            done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
            done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
            continue;

        }
    }

    LOG(INFO) ("Remove abnormal bridging edges: %d", removed.size()*2);
    graph_.ReduceEdges(removed, BaseEdge::RT_BRIDGED);
    LOG(INFO) ("Repair cross structures: %d", repaird.size()*2);
    Repair(repaird);
}
} // namespace fsa