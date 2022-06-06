#include "bubble.hpp"

namespace fsa {

bool BubbleSimplifier::ParseParameters(const std::vector<std::string> &params) {
    //assert(params[0] == "bubble");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "s") {
            simple = std::stoi(it[1]) == 1;
        } else {
            return false;
        }
    }
    return true;
}

// void BubbleSimplifier::Running() {
//     LOG(INFO)("BubbleSimplifier %d", simple);
//     graph_.FindBubbles(simple);
// }

void BubbleSimplifier::Running() {

    LOG(INFO)("BubbleSimplifier %d", simple);
    size_t thread_size = graph_.Options().thread_size;

    auto cands = graph_.CollectNodes([this](const SgNode* n) { return n->OutDegree() > 1; });
    std::vector<BubbleEdge*> bubbles(cands.size());
    MultiThreadMap(thread_size, cands, bubbles, [this](SgNode* n) {
        auto rr = Detect(static_cast<PathNode*>(n), simple);
        //auto rr = graph_.FindBubble(static_cast<PathNode*>(n), simple);
        if (rr != nullptr) {
            if (!rr->Validate(graph_.GetAsmData().GetInconsistentOverlaps(), 
                              graph_.GetReadStore(), 
                              graph_.string_graph_, 
                              graph_.Options().max_bubble_length)) {
                delete rr;
                rr =  nullptr;
            }

        }
        return rr;
    });

    std::vector<BubbleEdge*>  compound_path0;
    std::for_each(bubbles.begin(), bubbles.end(), [&compound_path0](BubbleEdge* e) { 
        if (e != nullptr) compound_path0.push_back(e); 
    });
    LOG(INFO)("End FindBubble %zd", compound_path0.size());

    std::sort(compound_path0.begin(), compound_path0.end(), [](const BubbleEdge* a, const BubbleEdge *b) { return a->simple_paths_.size() > b->simple_paths_.size(); });


    std::unordered_set<PathEdge*> edge_to_cpath;
    std::vector<BubbleEdge*>  compound_path1;
    std::unordered_set<PathEdge::ID, PathEdge::ID::Hash> path1_id;
    for (size_t i = 0; i < compound_path0.size(); ++i) {
        auto path = compound_path0[i];

        bool overlapped = false;
        std::list<PathEdge*> reverse_simple_paths;
        for (auto e : path->simple_paths_) {
            auto re = graph_.ReverseEdge(e);
            assert(re != nullptr);
            if (edge_to_cpath.find(e) == edge_to_cpath.end()) {
                edge_to_cpath.insert(e);
                edge_to_cpath.insert(re);
                reverse_simple_paths.push_back(re);
            } else {
                overlapped = true;
                break;
            }
        }

        if (!overlapped) {
            compound_path1.push_back(path);
            BubbleEdge* reverse_path = new BubbleEdge(graph_.ReverseNode(path->OutNode()), 
                                                      graph_.ReverseNode(path->InNode()), 
                                                      reverse_simple_paths, 
                                                      path->length_, 
                                                      path->width_, 
                                                      path->score_);
            reverse_path->IdentifySimplePaths(*graph_.string_graph_);
            assert(path1_id.find(reverse_path->Id() ) == path1_id.end());
            compound_path1.push_back(reverse_path);
            path1_id.insert(path->Id());
            path1_id.insert(reverse_path->Id());
        } else {
            delete path;
            compound_path0[i] = nullptr;    // TODO 使用智能指针
        }
    }
    
    for (auto path : compound_path1) {
        graph_.InsertEdge(path);
    }

    graph_.TestGraph();
}



bool IsClearBubble(const SgNode* start, const SgNode* end, const std::list<PathEdge*> &edges, const std::unordered_set<SgNode*> valids) {
    std::unordered_set<SgNode*> nodes;
    for (auto e : edges) {
        if (e->InNode() != start) {
            nodes.insert(e->InNode());
        }
        if (e->OutNode() != end) {
            nodes.insert(e->OutNode());
        }
    }

    for (auto n : nodes) {
        for (size_t i = 0; i < n->InDegree(); ++i) {
            auto ie = n->InEdge(i);
            auto in = ie->InNode();
            if (valids.find(in) == valids.end()) {
                return false;
            }
        }
    }
    return true;
}

BubbleEdge* BubbleSimplifier::Detect(PathNode* start_node, bool check, int depth_cutoff, int width_cutoff) {
    int length_cutoff = graph_.Options().max_bubble_length;

    SgNode* end_node = nullptr;

    std::list<PathEdge*> bundle_edges;     
    std::unordered_map<SgNode*, std::pair<int, int>> visited; // length, score

    std::list<SgNode*> local_node_list = graph_.GetEgoNodes(start_node, depth_cutoff);
    std::unordered_set<SgNode*> local_nodes(local_node_list.begin(), local_node_list.end());
    std::unordered_set<SgNode*> tips;

    int depth = 0;
    double width = 1.0;
    size_t length = 0;

    bool loop_detect = false;
    bool meet_error = false;
    bool spur = false;

    visited[start_node] = std::make_pair(0, 0);
    for (size_t i = 0; i < start_node->OutDegree(); ++i) {
        auto e = start_node->OutEdge(i);
        tips.insert(e->OutNode());
        bundle_edges.push_back(static_cast<PathEdge*>(e));
        length = std::max(length, e->Length());
        //LOG(INFO)("detect-bubble -- : %s, %zd, %s", start_node->Id().ToString(graph_.asmdata_.GetStringPool()).c_str(), e->Length(), e->Id().ToString(graph_.asmdata_.GetStringPool()).c_str());
    }



    do {
        std::unordered_map<SgNode*, SgEdge*> new_visited;     // 最新被访问节点，延后加入visited
        std::unordered_set<SgNode*> newtips, oldtips;     // 新产生的末梢节点和未处理的末梢节点
        
        for (auto n : tips) {
            //if (n->out_edges_.size() == 0) continue;        // dead end

            SgEdge* best_in_edge = nullptr;
            for (size_t i = 0; i < n->InDegree(); ++i) {
                auto e = n->InEdge(i);
                // 检查入节点，分成三类：不在局部集合中、已经访问、没有访问
                // 如果所有入节点都已经访问，则找出分数最高的边。并且可以扩展它的出节点
                // 否则改节点延后处理

                if (local_nodes.find(e->InNode()) != local_nodes.end()) {

                    if (visited.find(e->InNode()) != visited.end()) {
                        if (best_in_edge == nullptr || best_in_edge->Score() < e->Score()) {
                            best_in_edge = e;
                        }
                    }
                    else {
                        best_in_edge = nullptr;     // 
                        break;
                    }
                }
                else {
                    // 忽略这个入节点
                }
            }

            if (best_in_edge != nullptr) {

                assert(n == best_in_edge->OutNode());
                new_visited[n] = best_in_edge;

                // 如果气泡没有收敛，继续添加新的末梢节点
                if (tips.size() > 1) {
                    for (size_t i = 0; i < n->OutDegree(); ++i) {
                        auto e = n->OutEdge(i);
                        
                        if (visited.find(e->OutNode()) != visited.end() ||
                            new_visited.find(e->OutNode()) != new_visited.end()) {
                            loop_detect = true;
                            break;
                        }

                        SgNode *revese_node = graph_.ReverseNode(static_cast<PathNode*>(e->OutNode()));
                        if (local_nodes.find(e->OutNode()) != local_nodes.end() && 
                            visited.find(revese_node) == visited.end() &&
                            new_visited.find(revese_node) == new_visited.end() ) {

                            if (tips.find(e->OutNode()) == tips.end()) {
                                newtips.insert(e->OutNode());
                            }
                            bundle_edges.push_back(static_cast<PathEdge*>(e));
                        }
                        else {
                            meet_error = true;
                            break;
                        }
                    }

                    if (n->OutDegree() == 0) {
                        spur = true;
                        break;
                    }
                }
                else {
                    //end_node = n; // tips[0]
                    if (visited.find(n) != visited.end()) {
                        loop_detect = true;
                        end_node = nullptr;
                    } else {
                        end_node = n; // tips[0]
                    }
                }
            } else {
                if (tips.size() > 1) {
                    oldtips.insert(n);
                } else {
                    if (visited.find(n) != visited.end()) {
                        loop_detect = true;
                        end_node = nullptr;
                    } else {
                        end_node = n; // tips[0]
                    }
                }
            }

        }

        for (auto &i : new_visited) {
            visited[i.first] = std::make_pair(
                visited[i.second->InNode()].first + i.second->Length(),
                visited[i.second->InNode()].second + i.second->Score());

            // 更新当前长度
        //LOG(INFO)("detect-bubble up: %s, %zd, %zd", start_node->Id().ToString(graph_.asmdata_.GetStringPool()).c_str(), length, visited[i.first].first);
            if (length < visited[i.first].first) {
                length = visited[i.first].first;
            }

        }

       // LOG(INFO)("detect-bubble: %s, %zd", start_node->Id().ToString(graph_.asmdata_.GetStringPool()).c_str(), length);

        depth += 1;
        width = 1.0 * bundle_edges.size() / depth;


        if (tips.size() >= 1) {
            tips.clear();
            tips.insert(newtips.begin(), newtips.end());
            tips.insert(oldtips.begin(), oldtips.end());
        }

    } while (tips.size() >= 1 && tips.size() < 6 && !loop_detect && !meet_error && !spur && depth <= depth_cutoff && length <= length_cutoff && (depth <= 10 || width <= width_cutoff));

    if (end_node != nullptr && !loop_detect && !meet_error && !spur && depth <= depth_cutoff && length <= length_cutoff && (depth <= 10 || width <= width_cutoff) && (!check || check && IsClearBubble(start_node, end_node,bundle_edges, local_nodes))) {
        
        //LOG(INFO)("detect-bubble new: %s, %zd, %zd", start_node->Id().ToString(graph_.asmdata_.GetStringPool()).c_str(), visited[end_node].first, length);
        return new BubbleEdge(start_node, end_node, bundle_edges, visited[end_node].first, width, visited[end_node].second);
    }
    else {
        return nullptr;
    }
}

} // namespace fsa