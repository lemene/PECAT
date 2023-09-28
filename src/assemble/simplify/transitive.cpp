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
        return n->OutDegree() >= 2;
    });

	for (auto n : nodes) {
        assert(n->OutDegree() > 0);

		std::vector<BaseEdge*> &out_edges = n->GetOutEdges();

		std::sort(out_edges.begin(), out_edges.end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < b->Length(); });

		std::unordered_map<BaseNode*, BaseEdge*> node_infos;
		std::unordered_set<BaseEdge*> rmd;
		for (auto e : out_edges) {
			node_infos[e->OutNode()] = e;
		}

		for (auto e : out_edges) {
			Debug("cand: %s\n", ToString(e).c_str());
			
			BaseNode* w = e->OutNode();
			//if (rmd.find(e) == rmd.end()) {
				std::sort(w->GetOutEdges().begin(), w->GetOutEdges().end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < a->Length(); });
				for (auto e2 : w->GetOutEdges()) {
					auto e0 = node_infos.find(e2->OutNode());
					if (e0 != node_infos.end()) {
						Debug("check: %zd, %zd, %zd, %d, %d, %d\n", e0->second->Length(), e->Length(), e2->Length(), 
							graph_.GetAsmData().HasDup(e2->ol_->a_.id, e2->ol_->b_.id),
							graph_.GetAsmData().HasDup(e->ol_->a_.id, e->ol_->b_.id),
							graph_.GetAsmData().HasDup(e0->second->ol_->a_.id, e0->second->ol_->b_.id));

						bool consist = e0->second->Length() + fuzz_ >= e2->Length() + e->Length() && e0->second->Length() <= e2->Length() + e->Length() + fuzz_;
						if (!consist) {
							auto dup2 = graph_.GetAsmData().GetDup(e2->ol_->a_.id, e2->ol_->b_.id);
							auto dup = graph_.GetAsmData().GetDup(e->ol_->a_.id, e->ol_->b_.id);
							//auto dup0 = graph_.GetAsmData().GetDup(e0->second->ol_->a_.id, e0->second->ol_->b_.id);
							if (dup2.size() > 0 && dup.size() > 0) { // && dup0.size() > 0) {
								for (auto ol2 : dup2) {
									for (auto ol : dup) {
										//for (auto ol0 : dup0) {
											//if (Overlap::IsConsistent(*ol2, *ol, *ol0, fuzz_)) {
											if (Overlap::IsConsistent(*ol2, *ol, *(e0->second->ol_), fuzz_)) {
												consist = true;
												goto label_consist;
											}
										//}
									}
								}
								label_consist:;
							}
						}
						if (consist) {
						 	rmd.insert(e0->second);
						 	Debug("removed: %s, %s -> %s -> %s\n", ToString(e0->second).c_str(), ToString(e0->second).c_str(), ToString(e).c_str(), ToString(e2).c_str());

						}
						// if (e0->second->Length() + fuzz_ >= e2->Length() + e->Length() && e0->second->Length() <= e2->Length() + e->Length() + fuzz_  ||
						// 	graph_.GetAsmData().HasDup(e2->ol_->a_.id, e2->ol_->b_.id) ||
						// 	graph_.GetAsmData().HasDup(e->ol_->a_.id, e->ol_->b_.id) ||
						// 	graph_.GetAsmData().HasDup(e0->second->ol_->a_.id, e0->second->ol_->b_.id) ) {
						// 	rmd.insert(e0->second);
						// 	Debug("removed: %s, %s -> %s -> %s\n", ToString(e0->second).c_str(), ToString(e0->second).c_str(), ToString(e).c_str(), ToString(e2).c_str());
						// }
					}
				}
			//}
		}
		Debug("rmd_size %zd %zd\n", rmd.size(), out_edges.size());
		reduced.insert(rmd.begin(), rmd.end());

	}

    graph_.ReduceEdges(reduced, BaseEdge::RT_TRANSITIVE);
}


} // namespace fsa