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

    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() >= 2;
    });

    std::unordered_set<BaseEdge*> reduced;
    std::mutex mutex;
    auto combine_func = [this, &mutex, &reduced](const std::unordered_set<BaseEdge*> &edges) {
        std::lock_guard<std::mutex> lock(mutex);
		reduced.insert(edges.begin(), edges.end());
    };


    std::atomic<size_t> index {0};
    auto work_func = [this, &index, &nodes, combine_func](size_t tid) {
		std::unordered_set<BaseEdge*> removed;
        for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
			auto n = nodes[i];
			assert(n->OutDegree() > 0);

			auto rmd = SimplifyNode(n);
			removed.insert(rmd.begin(), rmd.end());
        }
        combine_func(removed);
    };

	auto work_func_sort = [&index, &nodes](size_t tid) {
        for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
			auto n = nodes[i];
			std::vector<BaseEdge*> &out_edges = n->GetOutEdges();
			std::sort(out_edges.begin(), out_edges.end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < b->Length(); });
        }
	};

	index.store(0);
	MultiThreadRun(graph_.Options().thread_size, work_func_sort);
	
	index.store(0);
	MultiThreadRun(graph_.Options().thread_size, work_func);
	LOG(INFO)("Transtive reduce %zd", reduced.size());
    graph_.ReduceEdges(reduced, BaseEdge::RT_TRANSITIVE);
}

std::unordered_set<BaseEdge*> TransitiveSimplifier::SimplifyNode(BaseNode* n) {
		std::vector<BaseEdge*> &out_edges = n->GetOutEdges();

		//std::sort(out_edges.begin(), out_edges.end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < b->Length(); });

		std::unordered_map<BaseNode*, BaseEdge*> node_infos;
		std::unordered_set<BaseEdge*> rmd;
		for (auto e : out_edges) {
			node_infos[e->OutNode()] = e;
		}

		for (auto e : out_edges) {
			//Debug("cand: %s\n", ToString(e).c_str());
			
			BaseNode* w = e->OutNode();
			//if (rmd.find(e) == rmd.end()) {
				//std::sort(w->GetOutEdges().begin(), w->GetOutEdges().end(), [](BaseEdge* a, BaseEdge *b) { return a->Length() < a->Length(); });
				for (auto e2 : w->GetOutEdges()) {
					auto e0 = node_infos.find(e2->OutNode());
					if (e0 != node_infos.end()) {
						// Debug("check: %zd, %zd, %zd, %d, %d, %d\n", e0->second->Length(), e->Length(), e2->Length(), 
						// 	graph_.GetAsmData().HasDup(e2->ol_->a_.id, e2->ol_->b_.id),
						// 	graph_.GetAsmData().HasDup(e->ol_->a_.id, e->ol_->b_.id),
						// 	graph_.GetAsmData().HasDup(e0->second->ol_->a_.id, e0->second->ol_->b_.id));

						bool consist = e0->second->Length() + fuzz_ >= e2->Length() + e->Length() && e0->second->Length() <= e2->Length() + e->Length() + fuzz_;
						if (!consist) {
							auto dup2 = graph_.GetAsmData().GetDup(e2->ol_->a_.id, e2->ol_->b_.id);
							auto dup = graph_.GetAsmData().GetDup(e->ol_->a_.id, e->ol_->b_.id);
							if (dup2.size() > 0 && dup.size() > 0) {
								for (auto ol2 : dup2) {
									for (auto ol : dup) {
										if (Overlap::IsConsistent(*ol2, *ol, *(e0->second->ol_), fuzz_)) {
											consist = true;
											goto label_consist;
										}
									}
								}
								label_consist:;
							}
						}
						if (consist) {
						 	rmd.insert(e0->second);
						 	//Debug("removed: %s, %s -> %s -> %s\n", ToString(e0->second).c_str(), ToString(e0->second).c_str(), ToString(e).c_str(), ToString(e2).c_str());

						}
					}
				}
			//}
		}
		//Debug("rmd_size %zd %zd\n", rmd.size(), out_edges.size());
		return rmd;
}

} // namespace fsa