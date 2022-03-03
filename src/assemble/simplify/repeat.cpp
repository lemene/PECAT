#include "repeat.hpp"

namespace fsa {



void RepeatSimplifier::Running() {
    std::unordered_set<BaseEdge*> edges_to_reduce;

    auto nodes_to_test = graph_.CollectNodes([](BaseNode* n) {
        return n->InDegree() == 1 && n->OutDegree() == 1;
    });
    
    for (auto n : nodes_to_test) {
        auto in_node = n->GetInEdges()[0]->InNode();
        auto out_node = n->GetOutEdges()[0]->OutNode();

        for (auto e : in_node->GetOutEdges()) {
            //auto vv = e->in_node_;
            auto ww = e->OutNode();

            auto ww_out_nodes = ww->GetAllOutNodes();
            auto v_out_nodes = n->GetAllOutNodes();
            bool overlap = HasCommon(ww_out_nodes, v_out_nodes);

            int ww_in_count = ww->InDegree();

            if (ww != n && !e->IsReduce() && ww_in_count > 1 && !overlap) {
                edges_to_reduce.insert(e);
            }

        }

        for (auto e : out_node->GetInEdges()) {
            auto vv = e->InNode();
            //auto ww = e->out_node_;

            auto vv_in_nodes = vv->GetAllInNodes();
            auto v_in_nodes = n->GetAllInNodes();
            bool overlap = HasCommon(vv_in_nodes, v_in_nodes);

            int vv_out_count = vv->OutDegree();

            if (vv != n && !e->IsReduce() && vv_out_count > 1 && !overlap) {
                edges_to_reduce.insert(e);
            }

        }
    }
    graph_.ReduceEdges(edges_to_reduce, BaseEdge::RT_REPEAT);
}


} // namespace fsa