#include "pspur.hpp"

namespace fsa {


bool Spur2Simplifier::ParseParameters(const std::vector<std::string> &params) { 
    //assert(params[0] == "bridge");

    CreateDebugFile();
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


// void Spur2Simplifier::Running() {

//     size_t depth_threshold = graph_.Options().max_spur_nodesize;
//     size_t length_threshold = graph_.Options().max_spur_length;

//     auto cands = graph_.CollectNodes([](SgNode* n) {
//         return n->InDegree() == 0;
//     });
    
//     for (auto n : cands) {
//         assert(n->InDegree() == 0);
        
//         std::vector<PathNode*> path;
//         PathNode *curr = static_cast<PathNode*>(n);
//         while (curr->InDegree() <= 1 && curr->OutDegree() == 1) {
//             path.push_back(curr);
//             curr = static_cast<PathNode*>(curr->OutNode(0));            
//             assert(curr != nullptr);
//         }
//         path.push_back(curr);   // 

//         if (path.size() >= 2 && path.back()->InDegree() >= 2 && path.back()->OutDegree() >= 1) {
//             size_t path_len = 0;
//             size_t path_node = 0;
//             for (size_t i = 0; i < path.size() - 1; ++i) {
//                 PathEdge* e = static_cast<PathEdge*>(path[i]->OutEdge(0));
//                 path_len += (graph_.ReverseEdge(e)->Length() + e->Length()) / 2;
//                 path_node += e->NodeSize();
//             }
            
//             if (path_len < length_threshold && path_node < depth_threshold) {

//                 bool is_spur = true;
//                 for (size_t i = 0; i < path.back()->InDegree(); ++i) {
//                     auto e = path.back()->InEdge(i);
                    
//                     if (e->InNode() != path[path.size()-2]) {
//                         if (!TestInExtend(e, path_len*2, path_node*2)) {
//                             is_spur = false;
//                             break;
//                         }
//                     }
//                 }

//                 if (is_spur) {
                    
//                     for (size_t i = 0; i + 1 < path.size(); ++i) {
//                         auto e = static_cast<PathEdge*>(path[i]->OutEdge(0));
//                         if (!e->IsReduced()) {
//                             e->Reduce("spur:2", true);
//                             graph_.ReverseEdge(e)->Reduce("spur:2", true);
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     LOG(INFO)("End Find spurs in the graph");
// }
void Spur2Simplifier::Running() {
    auto cands = graph_.CollectEdges([](SgEdge* e) {
        return e->InNode()->InDegree() == 0;
    });
    
    std::unordered_set<SgEdge*> removed;

    for (auto e : cands) {
        assert(e->InNode()->InDegree() == 0);
        
        std::vector<SgEdge*> path { e };
        size_t length = (e->Length() + graph_.SgGraph::ReverseEdge(e)->Length()) / 2;
        size_t nodesize = e->NodeSize();

        while (path.back()->OutNode()->InDegree() == 1 && path.back()->OutNode()->OutDegree() == 1) {
            path.push_back(path.back()->OutNode()->OutEdge(0));
            length += (path.back()->Length() + graph_.SgGraph::ReverseEdge(path.back())->Length()) / 2;
            nodesize += path.back()->NodeSize();
        }

        if (path.back()->OutNode()->InDegree() >= 2 && path.back()->OutNode()->OutDegree() >= 1 &&
            length < max_length && nodesize < max_nodesize) {


            bool is_spur = true;
            for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
                auto e = path.back()->OutNode()->InEdge(i);
                
                if (e != path.back()) {
                    if (!TestInExtend(e, length*2, nodesize)) {
                        is_spur = false;
                        break;
                    }
                }
            }

            if (is_spur) {
                for (auto e : path) {
                    removed.insert(e);

                }
            }
        }
    }

    for (auto e : removed) {
        auto pe = static_cast<PathEdge*>(e);
        if (!pe->IsReduced()) {
            pe->Reduce("spur:2", true);
            graph_.ReverseEdge(pe)->Reduce("spur:2", true);
        }
    }
}


bool Spur2Simplifier::TestInExtend(SgEdge* e, int minlen, int minnode) {
    int len = e->Length();
    int node = e->NodeSize();
    const SgNode* curr = e->InNode();

    while ((len < minlen || node < minnode) && curr->InDegree() >= 1) {//} && curr->OutDegree() == 1) {
        len += curr->InEdge(0)->Length();
        node += curr->InEdge(0)->NodeSize(); 
        curr = curr->InEdge(0)->InNode();
    }

    return len >= minlen && node >= minnode;
}


} // namespace fsa