#include "extend.hpp"

namespace fsa {

bool ExtendSimplifier::ParseParameters(const std::vector<std::string> &params) {
    assert(params[0] == "extend");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "s") {
        } else {
            return false;
        }
    }
    return true;
}

bool ExtendSimplifier::PreCondition() { 
    return !graph_.Options().read_file.empty(); 
}

void ExtendSimplifier::Running() {
    assert(PreCondition());

    ReactiveContainedEdges();
}

void ExtendSimplifier::ReactiveContainedEdges() {

    std::unordered_set<BaseNode*> done_cands;

    auto check_cand_node = [&done_cands](BaseNode* node) {
        const int N = 4;    // minimum path size

        size_t path_size = 0;
        if (node->OutDegree() == 0 && done_cands.find(node) == done_cands.end()) {
            for (SgNode* n = node; path_size < N && n->InDegree() >= 1; n = n->InNode(0)) {
                path_size ++;
            }
        }
        return path_size >= N;
    };
    

    for (size_t curr_iter = 0; curr_iter < iterations; ++curr_iter) {
        std::unordered_set<BaseEdge*> recovery;

        auto nodes = graph_.CollectNodes(check_cand_node);
        LOG(INFO)("candidate node size (%zd): %zd", curr_iter, nodes.size());
        done_cands.insert(nodes.begin(), nodes.end());

        
        std::vector<std::vector<Path>> all_paths(nodes.size());
         
        std::atomic<size_t> index { 0 };
        auto extend_func = [this, &index, &nodes, &all_paths](size_t tid) {
            for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
                auto n = nodes[i];
                auto& paths = all_paths[i];
                Debug("cand: %s\n", ToString(n).c_str());
    
                auto main = ExtendPath(n);
                if (!main.Empty() && graph_.GetNode(main.nodes.back().first)->OutDegree() > 0) {
                //if (!main.Empty()) {
                    bool main_valid = true;

                    ModifyMainPath(n, main);

                    if (!main.Empty() && graph_.GetNode(main.nodes.back().first)->OutDegree() > 0) {
                    //if (!main.Empty()) {
                        if (graph_.GetNode(main.nodes.back().first)->InDegree() > 0) {
                            auto secondary = GetSecondPath(n, main);
                            //if (!secondary.Empty() && graph_.GetNode(secondary.nodes.back().first)->OutDegree() > 0) {
                            if (!secondary.Empty()) {
                                Debug("Second path: %s: %zd\n", ToString(graph_.GetNode(secondary.nodes.back().first)).c_str(), secondary.nodes.size());
                                main_valid = AddSecondPath(main, secondary, paths);
                            }    
                        }
                        
                        if (main_valid) {
                            paths.push_back(main);
                            Debug("Add main path\n");
                        }
                    }

                } else {
                    Debug("Main path invalid: empty=%d\n", main.Empty());
                    if (!main.Empty()) Debug("Main path invalid: outdegree=%zd\n",graph_.GetNode(main.nodes.back().first)->OutDegree());
                }
            }
        };

        MultiThreadRun((size_t)graph_.Options().thread_size, extend_func);
        // for (auto n : nodes) {
        //     Debug("cand: %s\n", ToString(n).c_str());
    
        //     auto main = ExtendPath(n);
        //     if (!main.Empty() && graph_.GetNode(main.nodes.back().first)->OutDegree() > 0) {
        //         bool main_valid = true;

        //         ModifyMainPath(n, main);

        //         if (!main.Empty() && graph_.GetNode(main.nodes.back().first)->OutDegree() > 0) {
        //             auto secondary = GetSecondPath(n, main);
        //             if (!secondary.Empty() && graph_.GetNode(secondary.nodes.back().first)->OutDegree() > 0) {
        //                 Debug("Second path: %s: %zd\n", ToString(graph_.GetNode(secondary.nodes.back().first)).c_str(), secondary.nodes.size());
        //                 main_valid = AddSecondPath(main, secondary, paths);
        //             }
                    
        //             if (main_valid) {
        //                 paths.push_back(main);
        //                 Debug("Add main path\n");
        //             }
        //         }

        //     } else {
        //         Debug("Main path invalid: empty=%d\n", main.Empty());
        //         if (!main.Empty()) Debug("Main path invalid: outdegree=%zd\n",graph_.GetNode(main.nodes.back().first)->OutDegree());
        //     }
        // }
        std::vector<Path> paths;
        for (const auto& p: all_paths) {
            paths.insert(paths.end(), p.begin(), p.end());
        }
        LOG(INFO)("Recovery size2: %zd", paths.size());

        std::unordered_set<const Overlap*> repair;
        std::unordered_set<SgEdgeID, SgEdgeID::Hash> done;
        for (const auto& path : paths) {
            auto id = SgEdgeID(ID_EDGE_TYPE_MULT, SgNodeID(path.start), SgNodeID(path.nodes.back().first));
            Debug("Add path:%s, br=%zd, %d\n", id.ToString(graph_.GetAsmData().GetStringPool()).c_str(), path.nodes.size(),done.find(id) == done.end());
            if (done.find(id) == done.end()) {
                done.insert(id);
                done.insert(SgEdgeID::Reverse(id));

                for (const auto &e : path.nodes) {
                    repair.insert(e.second);
                    Debug("Add node: %s\n", graph_.GetAsmData().GetStringPool().QueryStringById(Seq::EndIdToId(e.first)).c_str());
                }
            }
            

        }
        
        LOG(INFO)("Recovery size2: %zd", repair.size()*2);
        Repair(repair);
    
    }
}


void ExtendSimplifier::ModifyMainPath(const BaseNode* s, Path &main) {
    auto end = graph_.GetNode(main.nodes.back().first);
    assert(end != nullptr);
    Debug("The first path: %s->%s(Indegree=%zd,OutDegree=%zd)\n", 
        ToString(s).c_str(), ToString(end).c_str(), end->InDegree(), end->OutDegree());

    if (end->InDegree() == 1 && graph_.GetAsmData().GetReadVariants() != nullptr) {
        std::unordered_set<Seq::Id> altnodes;
        auto curr = end->InNode<BaseNode>(0);    // skip the first element
        for (size_t i = 0; i < main.nodes.size(); ++i) {
            altnodes.insert(curr->ReadId());
            if (curr->InDegree() == 1) {
                curr = curr->InNode<BaseNode>(0);
            } else {
                break;
            }
        }

        Debug("The first altnodes nodesize: %zd\n", altnodes.size());
        main = ExtendPath(const_cast<BaseNode*>(s), altnodes);
        Debug("The first path(1): pathsize=%zd\n", main.nodes.size());
    }
}

auto  ExtendSimplifier::GetSecondPath(const BaseNode* s, const Path &main) -> Path {
    std::unordered_set<Seq::Id> altnodes;
    for (size_t i = 1; i < main.nodes.size(); ++i) {
        altnodes.insert(Seq::EndIdToId(main.nodes[i].first));
        Debug("  add alt: %s\n", ToString(Seq::EndIdToId(main.nodes[i].first)).c_str());
    }
    Debug("altnodes: %zd\n", altnodes.size());
    return ExtendPath(const_cast<BaseNode*>(s), altnodes); 
}

std::pair<Seq::EndId, const Overlap*> ExtendSimplifier::ExtendNode(Seq::EndId start, const std::unordered_set<Seq::Id>& altnodes, const std::unordered_set<Seq::Id>& visited) {

    double threshold = graph_.actual_min_identity;
    auto &asmdata = graph_.GetAsmData();

    Seq::Id sid = Seq::EndIdToId(start);
    int send = Seq::End(start);

    struct OlInfo {
        const Overlap* ol;
        size_t aligned;
        bool in_graph { false };
        Seq::Id eid { Seq::NID };
        int eend {0};
        BaseNode* node { nullptr };
        bool linked { false };
        std::array<int, 2> var_diff;
    };

    auto ols= asmdata.GetExtendOverlapsEx(sid, send);
    std::vector<OlInfo> infos(ols.size());
    std::transform(ols.begin(), ols.end(), infos.begin(), [sid, send, &altnodes, this](const Overlap* ol) {
        OlInfo its;
        its.ol = ol;
        its.aligned = ol->AlignedLength();
        its.eid = ol->GetOtherRead(sid).id;
        its.var_diff = CompareVariants(its.eid, altnodes);
        its.eend = send == 0 ? (ol->SameDirect() ? 0 : 1) : (ol->SameDirect() ? 1 : 0);
        its.node = graph_.GetNode(Seq::IdToEndId(its.eid, its.eend));
        its.linked = its.node != nullptr && its.node->OutDegree() > 0;
        return its;
    });


    std::sort(infos.begin(), infos.end(), [this](const OlInfo& a, const OlInfo& b) {

        bool sa = a.var_diff[1] >= a.var_diff[0];
        bool sb = b.var_diff[1] >= b.var_diff[0];

        if (sa && !sb) {
            return true;
        } else if (!sa && sb) {
            return false;
        } else {
            bool sa = a.linked && a.aligned >= acceptable_aligned_length;
            bool sb = b.linked && b.aligned >= acceptable_aligned_length;
            if (sa && !sb) {
                return true;
            } else if (!sa && sb) {
                return false;
            } else {
                return a.aligned > b.aligned;
            }            
        }
    });

    for (const auto& i : infos) {
        Debug("  check ExtendNode(%s):linked=(%d,%d) var=(%d, %d)\n", 
            asmdata.GetStringPool().QueryStringById(i.eid).c_str(), i.linked, i.node != nullptr,  i.var_diff[0], i.var_diff[1]  );
    }

    for (auto& i : infos) {
        auto &ol = i.ol; 
        auto it = asmdata.read_infos_.find(i.eid);
        assert(it != asmdata.read_infos_.end());

        if (visited.find(Seq::IdToEndId(i.eid, i.eend)) != visited.end()) {
            Debug(" Node(%s, %d) is visited\n", ToString(i.eid).c_str(), i.eend);
            continue;
        }

        auto qual = graph_.GetOverlapQuality(*ol);  
        Debug("  ExtendNode(%s): linked=(%d,%d), qual=%f, var=(%d,%d)\n", 
            asmdata.GetStringPool().QueryStringById(i.eid).c_str(), i.linked, i.node == nullptr, 
            qual, i.var_diff[0], i.var_diff[1] );
        
        if (qual < threshold) {
            auto nol = ReplaceHighQualityOverlap(ol, threshold);
            if (nol.first != nullptr) {
                i.ol = nol.first;
                qual = nol.second;
                ol = nol.first;
            }

        }

        if (qual >= threshold) {
            if (i.node != nullptr || it->second.filtered.type == RdReason::RS_CONTAINED || 
                (it->second.filtered.type == RdReason::RS_OK && i.node == nullptr)) {
                return {Seq::IdToEndId(i.eid, i.eend), ol};
            }
        }
    
    }

    return {0, nullptr};

}

auto ExtendSimplifier::ExtendPath(BaseNode* n, const std::unordered_set<Seq::Id>& cands) -> Path {
    assert(n->OutDegree() == 0);
    Debug("reactive contained cand:%d %s\n", n->Id().Value(0), ToString(n).c_str());
        
    Seq::EndId curr = n->Id().Value(0);
    
    Path path = {curr};

    std::unordered_set<Seq::Id> visited;
    auto prev_nodes = graph_.GetEgoNodes(graph_.ReverseNode(n), 30);
    for (auto pn : prev_nodes) {
        visited.insert(graph_.ReverseNode(static_cast<BaseNode*>(pn))->Id().Value(0));

    }
    Debug("visited size: %zd\n", visited.size());


    for (size_t i = 0; i < max_number_of_nodes; ++i) {
        Debug("ExtendNode curr:%d %s\n", curr, QueryStringById(Seq::EndIdToId(curr)).c_str());

        auto ext = ExtendNode(curr, cands, visited);
        
        if (ext.second != nullptr) {
            if (graph_.GetNode(ext.first) != nullptr && graph_.GetNode(ext.first)->OutDegree() > 0) {
                // improve overlap
                auto ps = ImproveLastEdge(ext);
                if (ps.size() > 0 && ext.first == ps.back().first) {
                    path.nodes.insert(path.nodes.end(), ps.begin(), ps.end());
                } else {
                    path.nodes.push_back(ext);
                }

                break;
            } else {
                path.nodes.push_back(ext);
                curr = ext.first;
                visited.insert(curr);
            }
        } else {
            break;
        }
    }

    if (path.nodes.size() > 0 && graph_.GetNode(path.nodes.back().first) == nullptr) {
        Debug("Invalid path: clear nodes\n");
        path.nodes.clear();
    }
    return path;
}


std::vector<std::pair<Seq::EndId, const Overlap*>> ExtendSimplifier::ImproveLastEdge(const std::pair<Seq::EndId, const Overlap*>& ol) {
    std::vector<std::pair<Seq::EndId, const Overlap*>> path;

    Seq::Id eid = Seq::EndIdToId(ol.first);
    int eend = Seq::End(ol.first);

    auto get_other_read_end = [](const Overlap* ol, int end) {
        return end == 0 ? (ol->SameDirect() ? 0 : 1) : (ol->SameDirect() ? 1 : 0);
    };

    Seq::Id sid = ol.second->GetOtherRead(eid).id;
    int send = get_other_read_end(ol.second, eend);
    Debug("Improve start: pathsize=%zd\n", path.size());


    auto ext_set = graph_.GetAsmData().GetExtendOverlaps(sid, send);
    ExtendLastEdge(ext_set, Seq::IdToEndId(sid, send), ol.first, path);
    Debug("Improve end: pathsize=%zd\n", path.size());
    //assert(path.size() > 0 && path.back().first == ol.first);

    return path;
}

void ExtendSimplifier::ExtendLastEdge(const std::unordered_set<const Overlap*> ext_set, Seq::EndId s, Seq::EndId end, std::vector<std::pair<Seq::EndId, const Overlap*>>& path) {
    Debug("get_extend: %zd, %zd\n", ext_set.size(), path.size());
    std::vector<const Overlap*> exts(ext_set.begin(), ext_set.end());

    std::sort(exts.begin(), exts.end(), [](const Overlap* a, const Overlap* b) {
        return a->AlignedLength() > b->AlignedLength();
    });

    auto &asmdata = graph_.GetAsmData();
    double threshold = graph_.actual_min_identity;
    Seq::Id curr_id = Seq::EndIdToId(s);
    int curr_end = Seq::End(s);

    auto get_other_read_end = [](const Overlap* ol, int end) {
        return end == 0 ? (ol->SameDirect() ? 0 : 1) : (ol->SameDirect() ? 1 : 0);
    };

    for (size_t i = 0; i < exts.size(); ++i) {
        auto next_id = exts[i]->GetOtherRead(curr_id).id;
        auto next_end = get_other_read_end(exts[i], curr_end);

        auto it = asmdata.read_infos_.find(next_id);
        assert(it != asmdata.read_infos_.end());

        // 避免循环
        if (std::find_if(path.begin(), path.end(), [next_id](const std::pair<Seq::EndId, const Overlap*> &a) {
            return Seq::EndIdToId(a.first) == next_id;
        }) != path.end()) {
            continue;
        };

        auto qual = graph_.GetOverlapQuality(*exts[i]);  
        Debug("  Improve(%s):  %f\n", QueryStringById(next_id).c_str(), qual);
        Debug("  Improve(%s->%s->%s):%zd  %f > %f\n", QueryStringById(curr_id).c_str(),QueryStringById(next_id).c_str(),QueryStringById(Seq::EndIdToId(end)).c_str(), i, qual,threshold);
        
        if (qual >= threshold) {
            if (end == Seq::IdToEndId(next_id, next_end)) {
                path.push_back({Seq::IdToEndId(next_id, next_end), exts[i]});
                break;
            } else if (it->second.filtered.type == RdReason::RS_CONTAINED) {
                auto o = asmdata.QueryOverlap(next_id, Seq::EndIdToId(end));
                if (o != nullptr) {
                    auto nnext = asmdata.GetExtendOverlaps(next_id, next_end);
                    if (nnext.find(o) != nnext.end()) {
                        path.push_back({Seq::IdToEndId(next_id, next_end), exts[i]});
                        Debug("nnext: %zd, %zd", nnext.size(), path.size());
                        ExtendLastEdge(nnext, path.back().first, end, path);
                        break;
                    }
                }
            }
        }
    }
    
}

Seq::Id ExtendSimplifier::GetReadId(const SgNode* n) const {
    assert(n != nullptr && n->IsType(BaseNode::TypeName()));
    auto bn = static_cast<const BaseNode*>(n);
    return bn->ReadId();
}

bool ExtendSimplifier::AddSecondPath(const Path& a, const Path& b, std::vector<Path>& paths) {
    assert(a.nodes.size() > 0 && b.nodes.size() > 0);

    auto enode_a = graph_.GetNode(a.nodes.back().first);
    auto enode_b = graph_.GetNode(b.nodes.back().first);
    assert(enode_a != nullptr && enode_b != nullptr);

    const int DEPTH = 30;
    Debug("PathBetter(end): %s(in=%zd,out=%zd) <-> %s(in=%zd,out(%zd)\n", 
        ToString(enode_a).c_str(), enode_a->InDegree(), enode_a->OutDegree(),
        ToString(enode_b).c_str(), enode_b->InDegree(), enode_b->OutDegree());

    if (enode_a != enode_b) {
        if (enode_a->InDegree() != 0 && enode_b->InDegree() != 0) {
            auto neighbor_a = graph_.GetNeighborNodes(enode_a, DEPTH);
            Debug("Find Neighbor nodes(depth=%d): %zd\n", DEPTH, neighbor_a.size());    
            if (std::find(neighbor_a.begin(), neighbor_a.end(), enode_b) != neighbor_a.end()) {
                Debug("b is one of a's neighbor\n");    
                
                auto next_a = graph_.GetEgoNodes(enode_a, DEPTH); 
                Debug("Find next nodes(depth=%d): %zd\n", DEPTH, next_a.size());   
                if (std::find(next_a.begin(), next_a.end(), enode_b) == next_a.end()) {
                    Debug("b is not one of a's subsequent nodes\n");    
                    paths.push_back(b);
                    Debug("Add second path\n");    
                    auto next_b = graph_.GetEgoNodes(enode_b, DEPTH);  
                    return std::find(next_b.begin(), next_b.end(), enode_a) == next_b.end();
                }
            }

        } else {
            if (enode_b->InDegree() == 0) {
                paths.push_back(b);
                Debug("Add second path\n"); 
            }

            return enode_a->InDegree() == 0;
        }
    }
    return true;
}

std::array<int,2> ExtendSimplifier::CompareVariants(Seq::Id q, const std::unordered_set<Seq::Id>& target) {
    std::array<int,2> result = {0, 0};
    auto rvs = graph_.GetAsmData().GetReadVariants();
    if (rvs != nullptr) {
        for (auto t : target) {
            auto r = rvs->Test(q, t);
            result[0] += r[0];
            result[1] += r[1];
            Debug("TestDiff: (%s, %s) = (%d, %d)\n", QueryStringById(q).c_str(), QueryStringById(t).c_str(), r[0], r[1]);

        }
    }
    return result;
}

void ExtendSimplifier::Repair(const std::unordered_set<const Overlap*> &ols) {

    auto get_edge = [this](const Overlap* ol) {
        Seq::EndId fB = Seq::IdToEndId(ol->a_.id, 0);
	    Seq::EndId fE = Seq::IdToEndId(ol->a_.id, 1);
	    Seq::EndId gB = Seq::IdToEndId(ol->b_.id, 0);
	    Seq::EndId gE = Seq::IdToEndId(ol->b_.id, 1);
        if (ol->SameDirect()) {
            if (ol->a_.start > 0) {
                return graph_.GetEdge(fE, gE);
            } else {
                return graph_.GetEdge(fB, gB);
            }
        } else {
            if (ol->a_.start > 0) {
                return graph_.GetEdge(fE, gB);
            } else {
                return graph_.GetEdge(fB, gE);
            }
        }
    };
    
    for (auto ol : ols) {

        auto e = get_edge(ol);
        if (e != nullptr) {
            graph_.ReactiveEdges(std::vector<BaseEdge*>({e}));

        } else {
            Debug("add_overlap %s, %s: %s\n", ToString(ol->a_.id).c_str(), ToString(ol->b_.id).c_str(), ol->ToM4Line().c_str());
            graph_.AddOverlap(ol);
        }

    }
}


} // namespace fsa
