#include "transitive.hpp"

namespace fsa {


bool TransitiveSimplifier::ParseParameters(const std::vector<std::string> &params) { 
    assert(params[0] == "transitive");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "fuzz") {
            fuzz_ = (size_t) std::stoul(it[1]);
        } else {
            return false;
        }
    }
    return true;
}

void TransitiveSimplifier::Running() {
    std::unordered_set<BaseEdge*> reduced;

    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() > 0;
    });

    for (auto n : nodes) { n->mark_ = 'V'; }

	for (auto n : nodes) {
        assert(n->OutDegree() > 0);

		std::vector<BaseEdge*> &out_edges = n->GetOutEdges();

		std::sort(out_edges.begin(), out_edges.end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < b->Length(); });

		for (auto e : out_edges) {
			e->OutNode()->mark_ = 'I';
		}

		int max_len = out_edges.back()->Length()  + fuzz_;

		for (auto e : out_edges) {
			BaseNode* w = e->OutNode();
			if (w->mark_ == 'I') {
				std::sort(w->GetOutEdges().begin(), w->GetOutEdges().end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < a->Length(); });
				for (auto e2 : w->GetOutEdges()) {
					if (e2->Length() + e->Length() < (size_t)max_len) {
						if (e2->OutNode()->mark_ == 'I') {
							e2->OutNode()->mark_ = 'E';
						}
					}
				}
			}
		}
		for (auto e : out_edges) {
			BaseNode* w = e->OutNode();
			std::sort(w->GetOutEdges().begin(), w->GetOutEdges().end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < b->Length(); });
			if (w->OutDegree() > 0) {
				//if (w->out_edges_.front()->OutNode()->mark_ == 'I') {
				//	w->out_edges_.front()->OutNode()->mark_ = 'E';
                //}
			}				
			for (auto e2 : w->GetOutEdges()) {
				if (e2->Length() < fuzz_) {
					if (e2->OutNode()->mark_ == 'I') {
						e2->OutNode()->mark_ = 'E';
					}
				}
			}
		}

		for (auto e : out_edges) {
			if (e->OutNode()->mark_ == 'E') {
                reduced.insert(e);
                //reduced.insert(StringGraph::ReverseEdge(e));
			}
			e->OutNode()->mark_ = 'V';
		}
	}

    graph_.ReduceEdges(reduced, BaseEdge::RT_TRANSITIVE);

}


} // namespace fsa