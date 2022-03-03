#include "loop.hpp"

namespace fsa {


void LoopSimplifier::Running() {
    FindLoops();
}



void LoopSimplifier::FindLoops() {

    std::vector<SgNode*> cands = graph_.CollectNodes([](SgNode *n) { return n->InDegree() == 2; });
    std::vector<NodeOrEdge> loops(cands.size());

    MultiThreadMap(graph_.Options().thread_size, cands, loops, [&](SgNode* n) {
        return DetectLoop(n);
    });

    LOG(INFO)("Found loop structures: (%zd, %zd)", 
        std::count_if(loops.begin(), loops.end(), [](const NodeOrEdge& i) {
            return i.node != nullptr;
        }),
        std::count_if(loops.begin(), loops.end(), [](const NodeOrEdge& i) {
            return i.edge != nullptr;
        })
    );

    for (auto i : loops) {
        if (i.edge != nullptr) {
            assert(i.node == nullptr);

            auto e = i.edge;
            
            if (graph_.QueryEdge(e->Id()) == nullptr) {
                assert(graph_.QueryEdge(e->ReverseId())  == nullptr);
                graph_.InsertEdge(e);
                graph_.InsertEdge(e->Reverse(graph_));
            } else {
                delete e;
            }
        } else if (i.node != nullptr) {
            // TODO
            auto n = i.node;
            LOG(WARNING)("TODO find loopnode");
            if (graph_.QueryNode(n->Id()) == nullptr) {
                assert(graph_.QueryNode(n->RId()) == nullptr);

                auto rn = n->Reverse(graph_);
                assert(rn != nullptr);
                graph_.InsertLoopNode(n);
                graph_.InsertLoopNode(rn);

            } else {
                delete n;
            }
        }

   
    }
}


// void LoopSimplifier::FindLoops2() {

//     std::vector<SgNode*> cands = graph_.CollectNodes([](SgNode *n) { return n->InDegree() == 2; });
//     std::vector<LoopNode*> loops(cands.size());

//     MultiThreadMap(graph_.Options().thread_size, cands, loops, [&](SgNode* n) {
//         return DetectLoop(n);
//     });
 
//     LOG(INFO)("Found loop structures: %zd", std::count_if(loops.begin(), loops.end(), [](const LoopNode* n) {
//         return n != nullptr;
//     }));

//     for (auto n : loops) {
//         if (n == nullptr) continue; 

//         if (graph_.QueryNode(n->Id()) == nullptr) {
//             assert(graph_.QueryNode(n->RId()) == nullptr);

//             auto rn = n->Reverse(graph_);
//             assert(rn != nullptr);
//             graph_.InsertLoopNode(n);
//             graph_.InsertLoopNode(rn);

//         } else {
//             delete n;
//         }
//     }
// }


auto LoopSimplifier::DetectLoop(SgNode* start) -> NodeOrEdge {
    assert(start->InDegree() == 2);
    if (start->OutDegree() == 1) {

        SgEdge* forward = start->OutEdge(0);
        auto end = forward->OutNode();
        if (end->OutDegree() == 2) {
            //       <- - <- 
            //       \      /
            //   -> - -> - -> -> 
            // start_node's outdegree == 1 意味着它只有唯一的延长方式
            // end_node's outdegree >= 2 找到到一个出边 能够返回 start_node，

            SgEdge* backward = nullptr;
            size_t bestlen = 0;
            for (size_t i = 0; i < end->OutDegree(); ++i) {
                auto e = end->OutEdge(i);
                if (e->OutNode() == start) {
                    size_t len = (e->Length() + graph_.ReverseEdge(static_cast<PathEdge*>(e))->Length()) / 2;
                    if (backward == nullptr || len < bestlen) {
                        backward = e;
                        bestlen = len;
                    }
                }
            }

            if (backward != nullptr){
                return  new LoopEdge(start, end, {forward}, {backward}) ;
            }

        } else if (end == start) {    // 
            //       <- <-
            //       \   /
            //   -> - -> 
            //  loop处于末端
            assert(start->OutDegree() == 1 && start->InDegree() == 2); // 前面条件暗含的条件
            
            return new LoopNode({forward}) ; // backward == forward
        } 

    } else if (start->OutDegree() == 2 ) { // out degree > 2 无法确定的 loop展开到哪个分支
        //       <- <-
        //       \   /
        //   -> - -> - -> 
        //  气泡处于中间，但start_node和end_node相同

        SgEdge* backward = nullptr;
        size_t count = 0;
        for (size_t i = 0; i < start->OutDegree(); ++i) {
            auto e = start->OutEdge(i);
            if (e->OutNode() == start) {
                backward = e;
                count ++;
            }
        }

        if (count == 1) {
            return new LoopNode({backward}) ;
        } 
    }

    return NodeOrEdge();  
}


//  寻找string graph中的loop，为了简化，loop个边的长度为 1 
//   
//   
LoopEdge* LoopSimplifier::FindLoop(SgNode* start_node, int max_loop_length) {
    assert(start_node->InDegree() == 2);
    if (start_node->OutDegree() == 1) {

        SgEdge* start_to_end = start_node->OutEdge(0);
        auto end_node = start_to_end->OutNode();
        if (end_node->OutDegree() >= 2) {
            //       <- - <- 
            //       \      /
            //   -> - -> - -> -> 
            // start_node's outdegree == 1 意味着它只有唯一的延长方式
            // end_node's outdegree >= 2 找到到一个出边 能够返回 start_node，

            SgEdge* end_to_start = nullptr;
            for (size_t i = 0; i < end_node->OutDegree(); ++i) {
                auto e = end_node->OutEdge(i);
                if (e->OutNode() == start_node) {
                    end_to_start = e;
                    break;
                }

            } 

            if (end_to_start != nullptr) {
                return new LoopEdge(start_node, end_node, {start_to_end}, {end_to_start});
            }
        } else if (end_node == start_node) {    // 
            //       <- <-
            //       \   /
            //   -> - -> 
            //  loop处于末端
            assert(end_node->OutDegree() == 1 && end_node->InDegree() == 2); // 前面条件暗含的条件

            return nullptr;
            // return new LoopEdge(start_node, end_node, {{start_node->out_edges_[0]}});
        } else {
            return nullptr;
        }
    } else if (start_node->OutDegree() == 2 ) { // out degree > 2 无法确定的 loop展开到哪个分支
        //       <- <-
        //       \   /
        //   -> - -> - -> 
        //  气泡处于中间，但start_node和end_node相同

        auto end_node = start_node;

        PathEdge* end_to_start = nullptr;
        for (size_t i = 0; i < end_node->OutDegree(); ++i) {
            auto e = end_node->OutEdge(i);
            if (e->OutNode() == start_node) {
                end_to_start = (PathEdge*)e;
                break;
            }
        }
    
        if (end_to_start != nullptr) {
            return nullptr;
        }
    }

    return nullptr;

}


} // namespace fsa