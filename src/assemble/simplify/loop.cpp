#include "loop.hpp"

namespace fsa {


void LoopSimplifier::Running() {
    FindLoops();
}


void LoopSimplifier::FindLoops() {

    std::vector<SgNode*> cands = graph_.CollectNodes([](SgNode *n) { return n->InDegree() == 2; });
    LOG(INFO)("Found cand loop structures: %zd", cands.size()); 

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


auto LoopSimplifier::DetectLoop(SgNode* start) -> NodeOrEdge {
    assert(start->InDegree() == 2);
    Debug("cand: %s %zd %zd\n", ToString(start).c_str(), start->InDegree(), start->OutDegree());

    if (start->OutDegree() == 1) {

        SgEdge* forward = start->OutEdge(0);
        auto end = forward->OutNode();
        Debug("cand end: %s %zd %zd\n", ToString(end).c_str(), end->InDegree(), end->OutDegree());
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


} // namespace fsa