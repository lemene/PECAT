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
        bool cand_cross = IsOutEdgeInconsistent(path.front()->InNode()->OutEdge<BaseEdge>(0), 
                                                path.front()->InNode()->OutEdge<BaseEdge>(1), 5, pif) &&
                          IsInEdgeInconsistent(path.back()->OutNode()->InEdge<BaseEdge>(0), 
                                               path.back()->OutNode()->InEdge<BaseEdge>(1), 5, pif);

        Debug("cand_cross: %d\n", cand_cross);
        if (cand_cross) {
            // 
            BaseNode* out_node = GetFrontAltNode(path);
            BaseNode* in_node = GetBackAltNode(path);

            if (out_node->InDegree() == 1 && out_node->OutDegree() == 1 && in_node->InDegree() == 1 && in_node->OutDegree() == 1) {
                auto e = graph_.QueryEdge(BaseEdge::ID(0, in_node->Id(),  out_node->Id()));
                if (e == nullptr) {
                    return GetAltCrossPath(in_node->Id().MainNode(), out_node->Id().MainNode(), path);
                } else if (e->IsType("base") && static_cast<BaseEdge*>(e)->GetReduceType() == BaseEdge::RT_NO_BEST) {
                    return std::vector<const Overlap*>({static_cast<BaseEdge*>(e)->ol_});  // TODO 需要更换更换参数
                    //return std::vector<const Overlap*>({nullptr});  // TODO 需要更换更换参数
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
        std::vector<std::pair<Seq::EndId, const Overlap*>> links;
        //std::unordered_map<Seq::EndId, const Overlap*> links;
        for (auto ol : ols) {
            auto& alt = ol->GetOtherRead(sid);
            auto rend = send == 0 ? (ol->SameDirect() ? 0 : 1) : (ol->SameDirect() ? 1 : 0);
            //links[Seq::IdToEndId(alt.id, rend)] = ol;
            links.push_back(std::make_pair(Seq::IdToEndId(alt.id, rend), ol));
        }
        return links;
    };

    auto check_overlap_quality = [this](const Overlap& ol) {
        auto qual = graph_.GetOverlapQuality(ol);
        Debug("qual: %f %f\n", qual, graph_.Options().min_identity);
        return qual >= graph_.Options().min_identity;
    };


    auto pif = graph_.GetAsmData().GetInconsistentOverlaps();
    assert(pif != nullptr);
    
    std::unordered_set<Seq::Id> bad;
   
    auto bs = pif->Get(Seq::EndIdToId(start));
    bad.insert(bs.begin(), bs.end());
    bs = pif->Get(Seq::EndIdToId(end));
    bad.insert(bs.begin(), bs.end());
    
    Debug("bad: %zd\n", bad.size());
    // TODO greed algorithm
    std::vector<const Overlap*> altpath;

    auto curr = start;
    for (size_t i = 0; i < 15; ++i) {
        auto extend = get_extend(curr);
        std::sort(extend.begin(), extend.end(), [](const std::pair<Seq::EndId, const Overlap*> &a,
                                                   const std::pair<Seq::EndId, const Overlap*> &b) {
            return a.second->AlignedLength() > b.second->AlignedLength();
        });
        Debug("extend %s %zd\n", asmdata.QueryNameById(Seq::EndIdToId(curr)).c_str(), extend.size());
        for (const auto &s : extend) {

            Debug("  ttt ext1: %s %zd\n", asmdata.QueryNameById(Seq::EndIdToId(s.first)).c_str(), s.second->AlignedLength());
        }
        auto iter = std::find_if(extend.begin(), extend.end(), [end](const std::pair<Seq::EndId, const Overlap*> &a){
            return a.first == end;
        });
        if (iter == extend.end()) {
            for (const auto &s : extend) {

                Debug("  ext1: %s %zd\n", asmdata.QueryNameById(Seq::EndIdToId(s.first)).c_str(), s.second->AlignedLength());
                auto rid = s.first;
                if (bad.find(Seq::EndIdToId(rid)) != bad.end()) continue;
                auto rinfo = asmdata.GetReadInfo(Seq::EndIdToId(rid));
                Debug("rinfo: %s\n", asmdata.QueryNameById(rinfo.id).c_str());
                if (rinfo.filtered.type == RdReason::RS_CONTAINED) {
                    if (check_overlap_quality(*(s.second))) {
                        altpath.push_back(s.second);
                        curr = rid;
                        break;
                    }
                }
            }
        } else {

            auto srinfo = asmdata.GetReadInfo(Seq::EndIdToId(curr));
            auto erinfo = asmdata.GetReadInfo(Seq::EndIdToId(end));
            if (erinfo.filtered.type == RdReason::RS_CONTAINED || srinfo.filtered.type == RdReason::RS_CONTAINED) {
                if (check_overlap_quality(*(iter->second)) ) {
                    altpath.push_back(iter->second);
                    curr = end;
                }

            }
            Debug("end: %s %d\n", asmdata.QueryNameById(Seq::EndIdToId(end)).c_str(), curr == end);
            break;
        }


    }
    if (curr != end) {
        altpath.clear();
    }

    return altpath;
    

}

//     b0->b1->b2
//        /
//  a0->a1->a2
// return a2
BaseNode* BridgeSimplifier::GetFrontAltNode(const std::vector<SgEdge*>& path) const {

    BaseNode* a2 = nullptr; // 参考注释
    assert(path.front()->InNode()->InDegree() == 1 && path.front()->InNode()->OutDegree() == 2);
    for (size_t i = 0; i < path.front()->InNode()->OutDegree(); ++i) {
        auto e = path.front()->InNode()->OutEdge(i);
        if (e == path.front()) continue;
        a2 = e->OutNode<BaseNode>();
        break;
    }
    assert(a2 != nullptr);
    return a2;
}

//     b0->b1->b2
//        /
//  a0->a1->a2
// return b0
BaseNode* BridgeSimplifier::GetBackAltNode(const std::vector<SgEdge*>& path) const {
        
    BaseNode* b0 = nullptr; // 参考注释
    assert(path.back()->OutNode()->InDegree() == 2 && path.back()->OutNode()->OutDegree() == 1);
    for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
        auto e = path.back()->OutNode()->InEdge(i);
        if (e == path.back()) continue;
        b0 = e->InNode<BaseNode>();
        break;
    }
    assert(b0 != nullptr);
    return b0;
}

//     b0->b1->b2
//        /
//  a0->a1->a2
// 如果a2和b0为inconsistent overlap，则删除a1->a1
bool BridgeSimplifier::IsBridgeDiploid(const std::vector<SgEdge*>& path) {

    const BaseNode* a2 = GetFrontAltNode(path); 
    const BaseNode* b0 = GetBackAltNode(path);
    
    auto pif = graph_.GetAsmData().GetInconsistentOverlaps();
    bool is_bridge_dipolid = pif != nullptr && pif->Contain(a2->ReadId(), b0->ReadId());
    Debug("TestBridgeDiplid: %s - %s: !nullptr=%d is_bridge=%d\n", ToString(a2).c_str(), ToString(b0).c_str(), pif != nullptr, is_bridge_dipolid);
    return is_bridge_dipolid;
}

bool BridgeSimplifier::IsAmbiguousPath(const std::vector<SgEdge*>& path) {

    auto path_core_length = [this](const std::vector<SgEdge*>& path) {
        int len = 0; 
        for (size_t i = 0; i < path.size(); ++i) {
            len += path[i]->Length();
        }
        len -= graph_.asmdata_.GetReadStore().GetSeqLength(static_cast<BaseEdge*>(path.back())->OutNode()->ReadId());
        return len;
    };

    int minlen = std::max(0, path_core_length(path))*2; // TODO how to determine the threshold

    //int minlen = graph_.PathLength(path) * 2;      // TODO how to determine the threshold
    int minnode = (path.size()+1) * 2;

    assert(path.front()->InNode()->InDegree() == 1 && path.front()->InNode()->OutDegree() == 2);
    std::array<bool, 2> in_node_has_long_edge = {false, false};
    for (size_t i = 0; i < path.front()->InNode()->OutDegree(); ++i) {
        auto e = path.front()->InNode()->OutEdge(i);
        if (e == path.front()) continue;

        in_node_has_long_edge[0] = static_cast<BaseEdge*>(path.front()->InNode()->InEdge(0))->TestInExtend(minlen, minnode,true);
        in_node_has_long_edge[1] = static_cast<BaseEdge*>(e)->TestOutExtend(minlen, minnode,true);
        break;
    }

    assert(path.back()->OutNode()->InDegree() == 2 && path.back()->OutNode()->OutDegree() == 1);
    std::array<bool, 2> out_node_has_long_edge = {false, false};
    for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
        auto e = path.back()->OutNode()->InEdge(i);
        if (e == path.back()) continue;

        out_node_has_long_edge[0] = static_cast<BaseEdge*>(path.back()->OutNode()->OutEdge(0))->TestOutExtend(minlen, minnode,true);
        out_node_has_long_edge[1] = static_cast<BaseEdge*>(e)->TestInExtend(minlen, minnode,true);
        break;
    }
    Debug("IsAmbiguousPath (%d,%d): %d %d %d %d\n", minlen, minnode, out_node_has_long_edge[0], out_node_has_long_edge[1], in_node_has_long_edge[0], in_node_has_long_edge[1]);
    return out_node_has_long_edge[0] && out_node_has_long_edge[1] && in_node_has_long_edge[0] && in_node_has_long_edge[1];            
}


bool BridgeSimplifier::IsLowQuality(const std::vector<SgEdge*> &path) {
    
    auto is_edge_low_quality = [this](SgEdge* e) {
        auto be = static_cast<BaseEdge*>(e);
        return graph_.GetOverlapQuality(*be->ol_) < graph_.desired_edge_quality;
    };

    bool low_qualtiy_path = false;
    for (auto p : path) {
        if (is_edge_low_quality(p)) {
            low_qualtiy_path = true;
            break;
        }
    }
    Debug("path_quality: %d\n", low_qualtiy_path);

    if (low_qualtiy_path) {
        assert(path.front()->InNode()->InDegree() == 1 && path.front()->InNode()->OutDegree() == 2);
        std::array<bool, 2> in_node_low_quality = {false, false};
        for (size_t i = 0; i < path.front()->InNode()->OutDegree(); ++i) {
            auto e = path.front()->InNode()->OutEdge(i);
            if (e == path.front()) continue;

            in_node_low_quality[0] = is_edge_low_quality(path.front()->InNode()->InEdge(0));
            in_node_low_quality[1] = is_edge_low_quality(e);
            break;
        }

        assert(path.back()->OutNode()->InDegree() == 2 && path.back()->OutNode()->OutDegree() == 1);
        std::array<bool, 2> out_node_low_quality = {false, false};
        for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
            auto e = path.back()->OutNode()->InEdge(i);
            if (e == path.back()) continue;

            out_node_low_quality[0] = is_edge_low_quality(path.back()->OutNode()->OutEdge(0));
            out_node_low_quality[1] = is_edge_low_quality(e);
            break;
        }
        Debug("check quality: %d %d %d %d\n", in_node_low_quality[0], in_node_low_quality[1], out_node_low_quality[0], out_node_low_quality[1]);
        return !in_node_low_quality[0] && !in_node_low_quality[1] && !out_node_low_quality[0] && !out_node_low_quality[1];
    } else {
        return false;
    }
   
}

int BridgeSimplifier::TestConsistent(const std::vector<SgEdge*>& path) {
    int s = 0;
    
    BaseEdge* ie = static_cast<BaseEdge*>(path.front());
    
    Debug("TestConsistent IN %d\n", ie->subject_);
    if (ie->subject_) {
        s ++;
    } 
    
    for (size_t i = 0; i < ie->InNode()->OutDegree(); ++i) {
        auto iie = ie->InNode()->OutEdge<BaseEdge>(i);
        if (iie == ie) continue;
        Debug("TestConsistent IN other %d\n", iie->subject_);
        if (iie->subject_) {
            s --;
        }
    }

    BaseEdge* oe = static_cast<BaseEdge*>(path.back());
    Debug("TestConsistent Out %d\n", ie->subject_);
    if (oe->subject_) {
        s ++;
    } 
    for (size_t i = 0; i < oe->OutNode()->InDegree(); ++i) {
        auto ioe = oe->OutNode()->InEdge<BaseEdge>(i);
        if (ioe == oe) continue;
        Debug("TestConsistent Out other %d\n", ioe->subject_);
        if (ioe->subject_) {
            s --;
        }
    }
    
    return s;
}

bool BridgeSimplifier::IsLinkedReversedNode(const std::vector<SgEdge*>& path, size_t max_depth) const {

    // collect reversed nodes
    auto curr = graph_.SgGraph::ReverseNode(path.front()->InNode());
    std::unordered_set<SgNode*> rnodes { curr };
    for (size_t i = 0; i < max_depth; ++i) {
        for (size_t i = 0; i < curr->OutDegree(); ++i) {
            rnodes.insert(curr->OutNode(i));
        }
        
        if (curr->OutDegree() == 1) {
            curr = curr->OutNode(0);
        } else {
            break;
        }
    }
    
    // check 
    curr = path.back()->OutNode();
    if (rnodes.find(curr) != rnodes.end()) {
        return true;
    }
    for (size_t i = 0; i < max_depth; ++i) {
        for (size_t i = 0; i < curr->OutDegree(); ++i) {
            if (rnodes.find(curr->OutNode(i)) != rnodes.end()) {
                return true;
            }
        }
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

        if (IsLinkedReversedNode(path, 36)) {
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

//        if (IsBridgeDiploid(path) || IsAmbiguousPath(path)) {
        if (IsAmbiguousPath(path)) { // || IsLowQuality(path)) {
            // auto r = TestConsistent(path);
            // if (r >= 2) {
                
            //     Debug("found Consistent path, %d\n", r);
            //     BaseEdge* ie = static_cast<BaseEdge*>(path.front());
            //     for (size_t i = 0; i < ie->InNode()->OutDegree(); ++i) {
            //         auto iie = ie->InNode()->OutEdge<BaseEdge>(i);
            //         if (iie != ie)  removed.insert(iie);
            //     }
                
            //     BaseEdge* oe = static_cast<BaseEdge*>(path.back());
            //     for (size_t i = 0; i < oe->OutNode()->InDegree(); ++i) {
            //         auto ioe = oe->OutNode()->InEdge<BaseEdge>(i);
            //         if (ioe != oe)  removed.insert(ioe);
            //     }
            // } else if (r <= 0) {
            //     Debug("Removed ambiguous path\n");
            //     for (auto e : path) {
            //         removed.insert((BaseEdge*)e);
            //     }
            // }
            Debug("Removed ambiguous path\n");
            for (auto e : path) {
                removed.insert((BaseEdge*)e);
            }
            done.insert(path.front());
            done.insert(path.back());
            done.insert(((SgGraph&)graph_).ReverseEdge(path.front()));
            done.insert(((SgGraph&)graph_).ReverseEdge(path.back()));
            continue;

        }
    }

    LOG(INFO) ("Remove abnormal bridging edges: %d", removed.size()*2);
    DebugEdges("rm", removed);
    graph_.ReduceEdges(removed, BaseEdge::RT_BRIDGED);
    LOG(INFO) ("Repair cross structures: %d", repaird.size()*2);
    Repair(repaird);
}
} // namespace fsa