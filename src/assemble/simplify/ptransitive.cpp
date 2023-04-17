#include "ptransitive.hpp"

namespace fsa {


bool PTransitiveSimplifier::ParseParameters(const std::vector<std::string> &params) { 
    assert(params[0] == "transitive");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "fuzz") {
        } else {
            return false;
        }
    }
    return true;
}

void PTransitiveSimplifier::Running() {
	for (size_t _ = 0; _ < 5; _++) {
		RemoveSmallLoops();
	}

	for (size_t _ = 0; _ < 5; _++) {
		Debug("iter(%zd)\n", _);
    std::unordered_set<SgEdge*> removed;
    auto nodes = graph_.CollectNodes([this](SgNode* n) {
        return n->OutDegree() > 1 && HasImportantInEdges(n);
    });
	LOG(INFO)("Candidate size: %zd", nodes.size());

	struct Item {
		SgNode* node;
		int max_len_out_edge;
	};
	std::vector<Item> node_items(nodes.size());
	std::transform(nodes.begin(), nodes.end(), node_items.begin(), [this](SgNode* n) ->Item {
		return {n, MaxCoreLengthOutEdges(n)};
	});

	std::sort(node_items.begin(), node_items.end(), [](const Item& a, const Item &b) {
		return a.max_len_out_edge < b.max_len_out_edge;
	});

	std::unordered_set<const SgNode*> involved;
	for (auto it : node_items) {
		auto n = it.node;
		bool is_cross = IsCross(n);
		bool is_involved = involved.find(n) != involved.end();

		Debug("cand: %s (IsCross=%d, Involved=%d)\n", ToString(n).c_str(), is_cross, is_involved);
		if (is_cross || is_involved) continue;

		std::vector<std::unordered_set<const SgNode*>> ends;
		std::unordered_set<const SgNode*> middles;
		bool failure = false;

		for (size_t i = 0; i < n->OutDegree(); ++i) {
			std::unordered_set<const SgNode*> m;
			ends.push_back(CollectImportantNodes(n->OutEdge(i), m));
			Debug("end(%zd): %zd\n", i, ends.back().size() );
			for (auto e : ends.back()) {
				Debug("end(-): %s\n", ToString(e).c_str());
			}
			middles.insert(m.begin(), m.end());
		}


		std::unordered_set<size_t> contained;
		for (size_t i = 0; i < ends.size(); ++i) {
			if (ends[i].size() == 0) {
				failure = true;
				break;
			}
			for (auto e : ends[i]) {
				Debug("End node(%zd, %s): %d\n", i, ToString(e).c_str(), involved.find(e) != involved.end());
				if (involved.find(e) != involved.end()) {
					failure = true;
					Debug("End node(%s) is involved\n", ToString(e).c_str());
					break;
				}
				Debug("End node\n");
			}
			if (failure) break;

			if (contained.find(i) != contained.end()) continue;
			for (size_t j = i+1; j < ends.size(); ++j) {
				if (contained.find(j) != contained.end()) continue;
				if (IsSubOf(ends[j], ends[i])) {
					contained.insert(j);
					Debug("End node include: %zd in %zd\n", j, i);
				} else if (IsSubOf(ends[i], ends[j])) {
					contained.insert(i);
					Debug("End node include: %zd in %zd\n", i, j);
				}
			}
		}

		if (!failure) {
			
			for (auto i : contained) {
				removed.insert(n->OutEdge(i));
				Debug("Add removed: %s\n", ToString(n->OutEdge(i)).c_str());
			}
			if (contained.size() > 0) {
				involved.insert(n);
				for (size_t i = 0; i < ends.size(); ++i) {
					for (auto n : ends[i]) {
						involved.insert(n);
						involved.insert(graph_.SgGraph::ReverseNode(n));
						
						Debug("Add end: %s\n", ToString(n).c_str());

					}
				}
				

				for (auto m : middles) {
					involved.insert(m);
					involved.insert(graph_.SgGraph::ReverseNode(m));
					
					Debug("Add middle: %s\n", ToString(m).c_str());

				}
				
			}
		}
		
	}

    std::unordered_set<SgEdge*> removed2;
	for (auto r : removed) {
		Debug("remove: %s\n", r->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
		auto path = GetExtendPath(r, removed, removed2);
		for (auto p : path) {
			removed2.insert(const_cast<SgEdge*>(p));
		} 
	}

	for (auto r : removed2) {
		Debug("remove2: %s\n", r->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());

		auto pe = static_cast<PathEdge*>(const_cast<SgEdge*>(r));
		if (!pe->IsReduced()) {
			pe->Reduce("ptransitive", true);
			graph_.ReverseEdge(pe)->Reduce("ptransitive", true);
		}
	}

	}
}

bool PTransitiveSimplifier::HasImportantInEdges(const SgNode* n) const {
	int max_len_in_edge = 0;
	for (size_t i = 0; i < n->InDegree(); ++i) {
		auto path = GetPath(n->InEdge(i));
		max_len_in_edge = std::max(max_len_in_edge, PathCoreLength(path));
	}

	int max_len_out_edge = MaxCoreLengthOutEdges(n);

	return max_len_in_edge > max_len_out_edge;
}

std::list<const SgEdge*> PTransitiveSimplifier::GetPath(const SgEdge* e) const {
	std::list<const SgEdge*> path {e};

	assert(e != nullptr);
	Debug("GetPath in %s, %zd, %zd, %zd, %zd\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(),
		e->InNode()->InDegree(), e->InNode()->OutDegree(), e->OutNode()->InDegree(), e->OutNode()->OutDegree() );
	const SgEdge *curr = e;
	while (curr->OutNode()->InDegree() == 1 && curr->OutNode()->OutDegree() == 1) {
		path.push_back(curr->OutNode()->OutEdge(0));
		curr = path.back();
		assert(curr != nullptr);
		Debug("GetPath in %s\n", curr->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
	}
	
	curr = e;
	Debug ("GetPath -- %zd %zd\n", curr->InNode()->InDegree(), curr->InNode()->OutDegree());
	while (curr->InNode()->InDegree() == 1 && curr->InNode()->OutDegree() == 1) {
		path.push_front(curr->InNode()->InEdge(0));
		curr = path.front();
		Debug("GetPath out %s\n", curr->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
	}
	Debug ("GetPath out %zd\n", path.size());

	return path;
}


std::unordered_set<const SgEdge*> PTransitiveSimplifier::GetExtendPath(const SgEdge* e, const std::unordered_set<SgEdge*> &rm0, const std::unordered_set<SgEdge*>& rm1) const {
	std::unordered_set<const SgEdge*> edges;

	auto is_dead_end = [&rm0, &rm1, &edges](const SgNode* node) {	
		for (size_t i = 0; i < node->InDegree(); ++i) {
			auto e = const_cast<SgEdge*>(node->InEdge(i));
			if (rm0.find(e) == rm0.end() && rm1.find(e) == rm1.end() && edges.find(e) == edges.end()) {
				return false;
			}
		}
		return true;
	};

	assert(e != nullptr);
	std::list<const SgEdge*> stack { e };
	
	while (stack.size() > 0 && stack.size() < 15) {
		auto curr = stack.front();
		stack.pop_front();
		auto path = GetPath(curr);
		if (IsTrivialPath(path)) {
			edges.insert(path.begin(), path.end());
			//if (path.back()->OutNode()->InDegree() == 1) {
			if (is_dead_end(path.back()->OutNode())) {
				for (size_t i = 0; i < path.back()->OutNode()->OutDegree(); ++i) {
					auto out_e = path.back()->OutNode()->OutEdge(i);
					if (edges.find(out_e) == edges.end()) {
						Debug("find extend path: found loop %s\n", ToString(out_e).c_str());
						stack.push_back(out_e);
					}
				}
			}
		}
	}
	if (stack.size() > 0) {
		Debug("find extend path: too deep %zd\n", stack.size());
		edges.clear();
	}

	return edges;
}

bool PTransitiveSimplifier::IsTrivialPath(const std::list<const SgEdge*> &path) const {
	size_t nodesize = SgGraph::PathNodeSize(path);
	int length = PathCoreLength(path);
	Debug("IsTrivialPath: length(%d < %zd) && nodesize(%zd < %zd) \n", length, max_length_, nodesize, max_nodesize_);
	return length < (int)max_length_ && nodesize < max_nodesize_;	// TODO
}


std::unordered_set<const SgNode*> PTransitiveSimplifier::CollectImportantNodes(const SgEdge* e, std::unordered_set<const SgNode*>& middle) const {
	std::unordered_set<const SgNode*> imp;


	std::unordered_set<const SgNode*> visited;

	std::list<const SgNode*> stack;
	if (IsTrivialPath(GetPath(e))) { 
		stack.push_back(e->OutNode()); 
		middle.insert(e->OutNode());
	}
	size_t count = 0;

	while (stack.size() > 0) {
		auto curr = stack.front();
		visited.insert(curr);
		stack.pop_front();
		for (size_t i = 0; i < curr->OutDegree(); ++i) {
			auto e = curr->OutEdge(i);
			if (IsTrivialPath(GetPath(e)) && !IsCross(const_cast<SgNode*>(curr))) {
				if (visited.find(e->OutNode()) == visited.end()) {
					stack.push_back(e->OutNode());
					middle.insert(e->OutNode());
				}
			} else {
				imp.insert(curr);
				if (IsTrivialPath(GetPath(e)) && IsCross(const_cast<SgNode*>(curr))) {
					if (visited.find(e->OutNode()) == visited.end()) {
						stack.push_back(e->OutNode());
						middle.insert(e->OutNode());
					}
				}
			}
		}
		Debug("stack: %zd\n", stack.size());
		count++;
		if (count > 20) break;
	}

	if (stack.size() > 0) {
		imp.clear();
		middle.clear();
	}

	return imp;
}


bool PTransitiveSimplifier::IsSubOf(const std::unordered_set<const SgNode*>& a, const std::unordered_set<const SgNode*>& b) const {
	for (auto i : a) {
		if (b.find(i) == b.end()) return false;
	}
	return true;
}

int PTransitiveSimplifier::MaxCoreLengthOutEdges(const SgNode* n) const {
	int max_len_out_edge = 0;
	for (size_t i = 0; i < n->OutDegree(); ++i) {
		auto path = GetPath(n->OutEdge(i));
		max_len_out_edge = std::max(max_len_out_edge, PathCoreLength(path));
	}
	return max_len_out_edge;
}

void PTransitiveSimplifier::RemoveSmallLoops() {
	
    std::unordered_set<SgEdge*> removed;
    auto nodes = graph_.CollectNodes([this](SgNode* n) {
        return n->InDegree() > 1;
    });

	for (auto n : nodes) {
		if (n->InDegree() <= 1) continue;	// be modified

		auto ext_list = graph_.GetEgoNodes(n, max_nodesize_, max_length_);
		
		std::unordered_set<SgEdge*> rmd;
		size_t count = 0;
		for (size_t i = 0; i < n->InDegree(); ++i) {
			auto ie = n->InEdge(i);
			if (std::find(ext_list.begin(), ext_list.end(), ie->InNode()) != ext_list.end()) {
				
				std::unordered_set<SgEdge*> rs;
				rs.insert(ie);
				auto curr = ie->InNode();
				while (curr->InDegree() == 1 && curr->OutDegree() == 1) {
					rs.insert(curr->InEdge(0));
					curr = curr->InNode(0);
				}
				if (curr->OutDegree() > 1) {
					rmd.insert(rs.begin(), rs.end());
					count ++;
				}
			}
		}
		if (count < n->InDegree()) {
			//removed.insert(rmd.begin(), rmd.end());
			for (auto r : rmd) {
				Debug("remove_small_loops: %s\n", r->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());

				auto pe = static_cast<PathEdge*>(const_cast<SgEdge*>(r));
				if (!pe->IsReduced()) {
					pe->Reduce("small_loops", true);
					graph_.ReverseEdge(pe)->Reduce("small_loops", true);
				}
	}
		}
	}

	// for (auto r : removed) {
	// 	Debug("remove_small_loops: %s\n", r->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());

	// 	auto pe = static_cast<PathEdge*>(const_cast<SgEdge*>(r));
	// 	if (!pe->IsReduced()) {
	// 		pe->Reduce("small_loops", true);
	// 		graph_.ReverseEdge(pe)->Reduce("small_loops", true);
	// 	}
	// }
}

} // namespace fsa