#include "string_graph.hpp"

#include <algorithm>
#include <cstdio>
#include <cassert>
#include <list>
#include <unordered_set>
#include <iterator> 
#include <iomanip>

#include "asm_dataset.hpp"
#include "../overlap.hpp"
#include "../overlap_store.hpp"
#include "../utils/logger.hpp"
#include "../file_io.hpp"
#include "edlib.h"
#include "simplify/simplifier.hpp"
#include "simplify/spur.hpp"
#include "simplify/pspur.hpp"
#include "simplify/bridge.hpp"
#include "simplify/pbridge.hpp"
#include "simplify/transitive.hpp"
#include "simplify/cross.hpp"
#include "simplify/low_quality.hpp"
#include "simplify/unreliable.hpp"
#include "simplify/loop.hpp"
#include "simplify/bubble.hpp"
#include "simplify/semi_bubble.hpp"
#include "simplify/best.hpp"
#include "simplify/repeat.hpp"
#include "simplify/phase.hpp"

namespace fsa {
SgGraph::~SgGraph() {
    for (auto& i : org_nodes_) {
        delete i.second;
    }

    for (auto& i : org_edges_) {
        delete i.second;
    }
}
    
std::vector<SgEdge*> SgGraph::GetLinearPath(SgEdge* start, size_t max_length, size_t max_nodesize) {
    std::vector<SgEdge*>  path { start };

    size_t vlength = ReverseEdge(start)->Length();
    size_t length = start->Length();
    size_t nodesize = start->NodeSize();
    SgNode* curr = start->OutNode();
    assert(nodesize > 0);

    while (vlength < max_length && length < max_length && nodesize < max_nodesize && curr->InDegree() == 1 && curr->OutDegree() == 1) {
        SgEdge* next = curr->OutEdge(0);
        path.push_back(next);
        vlength += ReverseEdge(next)->Length();
        length += next->Length();
        nodesize += next->NodeSize(); 
        curr = next->OutNode();
    }

    if (vlength < max_length && length < max_length && nodesize < max_nodesize) {
        return path;
    } else {
        return std::vector<SgEdge*>();
    }
}

    
size_t SgGraph::PathLength(const std::vector<SgEdge*> &path) {
    size_t len = 0;
    for (auto e : path) {
        len += ReverseEdge(e)->Length();
        len += e->Length();
    }
    return len / 2;
}

std::list<SgNode*> SgGraph::GetEgoNodes(SgNode* n, int max_depth) {
    std::list<SgNode*> nodes{ n };
    std::unordered_set<SgNode*> nodes_set {n};

    int depth = 0;
    auto curr = nodes.begin();
    auto level_end = nodes.end();
    level_end--;

    while (depth < max_depth && curr != nodes.end()) {

        for (size_t i = 0; i < (*curr)->OutDegree(); ++i) {
            auto e = (*curr)->OutEdge(i);

            if (nodes_set.find(e->OutNode()) == nodes_set.end()) {
                nodes.push_back(e->OutNode());
                nodes_set.insert(e->OutNode());
            }
        }
        
        if (curr == level_end) {
            depth++;
            level_end = nodes.end();
            level_end--;
        }
        curr++;
    }


    return std::list<SgNode*>(nodes.begin(), curr); 
}

std::list<SgNode*> SgGraph::GetEgoNodes(SgNode* n, int max_depth, int max_len) {
    std::list<SgNode*> nodes{ n }; 
    std::list<int> lens{ 0 };

    int depth = 0;
    auto curr = nodes.begin();
    auto curr_len = lens.begin();

    auto level_end = nodes.end();
    level_end--;

    while (depth < max_depth && curr != nodes.end()) {

        for (size_t i = 0; i < (*curr)->OutDegree(); ++i) {
            auto e = (*curr)->OutEdge(i);

            if (std::find(nodes.begin(), nodes.end(), e->OutNode()) == nodes.end() && 
                *curr_len + e->Length() < max_len) {
                nodes.push_back(e->OutNode());
                lens.push_back(*curr_len + e->Length());
            }
        }

        if (curr == level_end) {
            depth++;
            level_end = nodes.end();
            level_end--;
        }
        curr++;
        curr_len++;
    }


    return std::list<SgNode*>(nodes.begin(), curr);
}


std::list<SgNode*> SgGraph::GetEgoNodes(SgNode* n, int max_depth, int max_len, int max_nodesize) {
    std::list<SgNode*> nodes{ n }; 
    std::list<int> lens{ 0 };

    int nodesize = 0;
    int depth = 0;
    auto curr = nodes.begin();
    auto curr_len = lens.begin();

    auto level_end = nodes.end();
    level_end--;

    while (depth < max_depth && curr != nodes.end()) {

        for (size_t i = 0; i < (*curr)->OutDegree(); ++i) {
            auto e = (*curr)->OutEdge(i);

            if (std::find(nodes.begin(), nodes.end(), e->OutNode()) == nodes.end() && 
                *curr_len + e->Length() < max_len) {
                nodes.push_back(e->OutNode());
                lens.push_back(*curr_len + e->Length());
                nodesize += e->NodeSize();
            }
        }

        if (curr == level_end) {
            depth++;
            level_end = nodes.end();
            level_end--;

            if (nodesize > max_nodesize) break;
        }
        curr++;
        curr_len++;
    }


    return std::list<SgNode*>(nodes.begin(), curr);
}

std::vector<std::vector<SgEdge*>> SgGraph::MaximumFlow(const SgNode* src, const SgNode *dst, const std::unordered_set<const SgEdge*> &_edges) {

    std::unordered_map<const SgEdge*, int> edges;
    for (auto e : _edges) {
        edges[e] = 1;
    }
    
    struct Work {
        SgNode* n;
        size_t i;
    };

    auto find_path = [&]() {
        std::vector<const SgEdge*> path;
        std::unordered_set<const SgNode*> done;

        std::vector<Work> stack;
        stack.push_back({const_cast<SgNode*>(src), 0});
        done.insert(src);

        while (stack.size() > 0 && stack.back().n != dst) {
            auto& curr = stack.back();
            bool found = false;
            if (curr.i < curr.n->InDegree()) {
                for (; curr.i < curr.n->InDegree(); ++curr.i) {
                    auto ie = curr.n->InEdge(curr.i);
                    if (done.find(ie->InNode()) != done.end()) continue;
                    
                    auto iter = edges.find(ie);
                    if (iter != edges.end() && iter->second == 0) {
                        stack.push_back({ie->InNode(), 0});
                        done.insert(ie->InNode());
                        found = true;
                        break;
                    }
                    
                }
            }

            if (!found && curr.i >= curr.n->InDegree()) {
                for (; curr.i < curr.n->InDegree() + curr.n->OutDegree(); ++curr.i) {
                    auto oe = curr.n->OutEdge(curr.i-curr.n->InDegree());
                    if (done.find(oe->OutNode()) != done.end()) continue;
                    
                    auto iter = edges.find(oe);
                    if (iter != edges.end() && iter->second == 1) {
                       
                        stack.push_back({oe->OutNode(), 0});
                        done.insert(oe->OutNode());
                        found = true;
                        break;
                    }
                    
                }

            }

            if (!found) {
                done.erase(curr.n);
                stack.pop_back();
                if (stack.size() > 0) stack.back().i++;
            }
        }

        if (stack.size() > 0) {
            std::vector<const SgEdge*> path;
            for (size_t i=0; i<stack.size()-1; ++i){
                if (stack[i].i < stack[i].n->InDegree()) {
                    path.push_back(stack[i].n->InEdge(stack[i].i));
                } else {
                    path.push_back(stack[i].n->OutEdge(stack[i].i-stack[i].n->InDegree()));
                }
            }
            return path;
        } else {
            return std::vector<const SgEdge*>();
        }
    };

    std::vector<const SgEdge*> path = find_path();
    while (path.size() > 0) {
        const SgNode* curr = src;
        for (auto p : path) {
            if (p->InNode() == curr) {
                assert(edges[p] == 1);
                edges[p] = 0;
                curr = p->OutNode();
            } else {
                assert(p->OutNode() == curr);
                assert(edges[p] == 0);
                edges[p] = 1;
                curr = p->InNode(); 
            }
        }
        assert(curr == dst);
        path = find_path();
    }


    std::vector<std::vector<SgEdge*>> paths;
    std::unordered_set<const SgEdge*> done;

    std::vector<Work> stack;
    stack.push_back({const_cast<SgNode*>(src), 0});

    while (stack.size() > 0 && stack.back().n != dst) {
        auto& curr = stack.back();
        //printf("stack %zd, %zd, %d %zd\n", stack.size(), edges.size(), curr.n->Id(), curr.i); fflush(stdout);

        for (; curr.i < curr.n->OutDegree(); ++curr.i) {
            auto oe = curr.n->OutEdge(curr.i);
            
            auto iter = edges.find(oe);
            if (iter != edges.end() && iter->second == 0 && done.find(oe) == done.end()) {
                stack.push_back({oe->OutNode(), 0});
                break;
            }
        }



        if (stack.size() > 0 && stack.back().n == dst) {
            std::vector<SgEdge*> path;
            for (size_t i=0; i<stack.size()-1; ++i){
                assert(stack[i].i < stack[i].n->OutDegree());
                path.push_back(stack[i].n->OutEdge(stack[i].i));
                done.insert(path.back());
            }
            paths.push_back(path);
        } 

        assert(stack.size() > 0);
        if (stack.back().n == dst || stack.back().i >= stack.back().n->OutDegree()) {
            stack.pop_back();
            if (stack.size() > 0) stack.back().i++;
        }
    }

    return paths;
}


void SgGraph::Simplify(const std::string &strategy, const std::string &reducers_str) {

    std::unordered_map<std::string, std::vector<std::string>> params;
    auto reducers = SplitStringByChar(reducers_str, '|');
    for (const auto &i : reducers) {
        auto items = SplitStringByChar(i, ':');
        assert(items.size() >= 1);
        params[items[0]] = items;

    }

    std::unordered_map<std::string, size_t> index;
    auto actions = SplitStringByChar(strategy, '|');

    for (const auto &i : actions) {
        auto items = SplitStringByChar(i, ':');
        assert(items.size() >= 1);
        auto simplifier = simplifiers_.find(items[0]);
        if (simplifier != simplifiers_.end()) {
            auto i = index[items[0]];
            std::string n = items[0] + "-" + std::to_string(i);
            auto p = params.find(n);
            if (p == params.end()) {
                p = params.find(items[0]);
            }

            if (p != params.end()) {
                items.insert(items.end(), p->second.begin()+1, p->second.end());
            }

            simplifier->second->Simplify(items);

            index[items[0]]++;
        } else {
            LOG(ERROR)("Failed to recognize string graph simplification strategy: %s", i.c_str());
        }
    }
}

StringGraph::StringGraph(AsmDataset &asmdata) 
     : SgGraph(asmdata)
     , nodes_(reinterpret_cast<std::unordered_map<BaseNode::ID, BaseNode*, BaseNode::ID::Hash> &>(org_nodes_))
     , edges_(reinterpret_cast<std::unordered_map<BaseEdge::ID, BaseEdge*, BaseEdge::ID::Hash> &>(org_edges_)) {
//{
    simplifiers_["transitive"].reset(new TransitiveSimplifier(*this));
    simplifiers_["spur"].reset(new SpurSimplifier(*this));
    simplifiers_["best"].reset(new BestOverlapsSimplifier(*this));
    simplifiers_["bridge"].reset(new BridgeSimplifier(*this));
    simplifiers_["phase"].reset(new PhaseCrossSimplifier(*this));
    simplifiers_["quality"].reset(new LowQualitySimplifier(*this));
    simplifiers_["repeat"].reset(new RepeatSimplifier(*this));
    simplifiers_["unreliable"].reset(new UnreliableSimplifier(*this));

}

StringGraph::~StringGraph() {
}

std::string StringGraph::IdToString(BaseNode::ID id) const {
    if (id.MainNode() > 0) {
        return asmdata_.GetReadStore().QueryNameById(id.MainNode()-1) + ":B";
    }
    else if(id.MainNode() < 0){
        return asmdata_.GetReadStore().QueryNameById(-id.MainNode()-1) + ":E";
    }
    else {
        return "NA";
    }

}

void StringGraph::AddOverlap(const Overlap* overlap) {

	Seq::EndId fB = Seq::IdToEndId(overlap->a_.id, 0);
	Seq::EndId fE = Seq::IdToEndId(overlap->a_.id, 1);
	Seq::EndId gB = Seq::IdToEndId(overlap->b_.id, 0);
	Seq::EndId gE = Seq::IdToEndId(overlap->b_.id, 1);
    
	if (overlap->a_.start > 0) {
        assert(overlap->a_.end == overlap->a_.len);
		if (overlap->SameDirect()) {
			//  f.B         f.E
			//	f----------->
			//	g       ------------->
			//	      g.B           g.E
            assert(overlap->b_.start == 0 && overlap->b_.end < overlap->b_.len);
           
			AddEdge(overlap, fE, gE, overlap->b_.id);
			AddEdge(overlap, gB, fB, overlap->a_.id);
		} else {
			// f.B         f.E
			//	f----------->
			//	g         <------------ -
			//	          g.E           g.B
            assert(overlap->b_.start > 0 && overlap->b_.end == overlap->b_.len);

			AddEdge(overlap, fE, gB, overlap->b_.id);
			AddEdge(overlap, gE, fB, overlap->a_.id);
			
		}
	} else {
        assert(overlap->a_.end < overlap->a_.len);
		if (overlap->SameDirect()) {
			
			//       f.B         f.E
			//  f     ----------->
			//	g------------->
			//	g.B           g.E
            assert(overlap->b_.start > 0 && overlap->b_.end == overlap->b_.len);
			AddEdge(overlap, fB, gB, overlap->b_.id);
			AddEdge(overlap, gE, fE, overlap->a_.id);
		}
		else {
			//        f.B         f.E
			// f       ----------->
			//	g <------------ 
			//	g.E           g.B
            assert(overlap->b_.start == 0 && overlap->b_.end < overlap->b_.len);
			AddEdge(overlap, fB, gE, overlap->b_.id);
			AddEdge(overlap, gB, fE, overlap->a_.id);
		}
	}

}

void StringGraph::AddOverlaps(const OverlapStore &ol_store) {

    std::unordered_set<std::array<Seq::Id, 2>, ArrayHash<int,2>, ArrayEqual<int,2>> done;
    for (size_t i=0; i<ol_store.Size(); ++i) {
        const auto& o = ol_store.Get(i);
        auto ids = std::minmax(o.a_.id, o.b_.id);
        std::array<Seq::Id, 2> id_pair = {ids.first, ids.second};

        if (done.find(id_pair) == done.end()) {

            if (o.attached == 0) {
                AddOverlap(&(o));
                done.insert(id_pair);

            }
        }
    }
    LOG(INFO)("Done = %zd", done.size());
}

void StringGraph::AddOverlaps(const OverlapStore &ol_store, int min_length, int min_aligned_lenght, float min_identity) {

    std::unordered_set<BaseNode::ID, BaseNode::ID::Hash> contained;

    for (size_t i=0; i<ol_store.Size(); ++i) {
        const auto& ol = ol_store.Get(i);
        if (ol.IsContaining(0)) {
            contained.insert(ol.b_.id);
        } else if (ol.IsContained(0)) {
            contained.insert(ol.a_.id);
        }
    }

    LOG(INFO)("contained = %zd", contained.size());
    std::unordered_set<std::array<Seq::Id, 2>, ArrayHash<int,2>, ArrayEqual<int,2>> done;
    for (size_t i=0; i<ol_store.Size(); ++i) {
        const auto& o = ol_store.Get(i);
        auto ids = std::minmax(o.a_.id, o.b_.id);
        std::array<Seq::Id,2> id_pair = {ids.first, ids.second};

        if (done.find(id_pair) == done.end()) {

            if (o.attached == 0 && !FilterOverlap(o, contained, min_length, min_aligned_lenght, min_identity)) {
                AddOverlap(&(o));
                done.insert(id_pair);

            }
        }
    }
    LOG(INFO)("Done = %zd", done.size());
}


bool StringGraph::FilterOverlap(const Overlap &ovlp, const std::unordered_set<BaseNode::ID, BaseNode::ID::Hash> &contained, int min_length, int min_aligned_length, float min_identity) {
    if (contained.find(ovlp.a_.id) != contained.end() || contained.find(ovlp.b_.id) != contained.end())
        return true;

    if (ovlp.a_.id == ovlp.b_.id) 
        return true;

    if (!ovlp.IsProper(0))
        return true;
    if (ovlp.identity_ < min_identity)
        return true;
    if (ovlp.a_.len < min_length || ovlp.b_.len < min_length)
        return true;

    if (ovlp.AlignedLength() < (size_t)min_aligned_length)
        return true;
    return false;
}

void StringGraph::AddEdge(const Overlap* ol, int in_node, int out_node, int read) {
	auto in = nodes_.find(in_node);
	if (in == nodes_.end()) {
		auto r = nodes_.insert(std::make_pair(in_node, new BaseNode(in_node)));
		assert(r.second);
		in = r.first;
	}

	auto out = nodes_.find(out_node);
	if (out == nodes_.end()) {
		auto r = nodes_.insert(std::make_pair(out_node, new BaseNode(out_node)));
		assert(r.second);
		out = r.first;
	}
	
	BaseEdge *e = new BaseEdge(in->second, out->second);
	edges_[BaseEdge::ID(0, in_node, out_node)] = e;
    e->read_ = read;
    e->ol_ = ol;

    in->second->AddOutEdge(e);
    out->second->AddInEdge(e);
}

std::unordered_set<BaseNode*> StringGraph::BfsNodes(BaseNode* n, BaseNode *exclude, int depth) {
	std::unordered_set<BaseNode*> result;
	result.insert(n);

	std::list<BaseNode*> cand;
	cand.push_back(n);
    BaseNode* depth_node = n;
	int dp = 1;
	while (dp < depth && cand.size() > 0) {
		BaseNode* v = cand.front();
		cand.pop_front();

		for (auto e : v->GetOutEdges()) {
			if (e->OutNode() != exclude) {
                if (result.find(e->OutNode()) == result.end()) {
                    result.insert(e->OutNode());
                    if (e->OutNode()->OutDegree() > 0) {
                        cand.push_back(e->OutNode());
                    }
                }
			}
		}

        for (auto e : v->GetReducedOutEdge()) {
            if (e->OutNode() != exclude) {
                if (result.find(e->OutNode()) == result.end()) {
                    result.insert(e->OutNode());
                    if (e->OutNode()->OutDegree() > 0) {
                        cand.push_back(e->OutNode());
                    }
                }
            }
        }
        if (v == depth_node) {

		    dp++;
            depth_node = cand.size() > 0 ? cand.back() : nullptr;
        }

	}
	return result;
}


StringGraph::LinearPath StringGraph::FindBridgePath(BaseEdge* start, int max_length, int max_nodesize) {
    
    std::vector<BaseEdge*> edge = {start};
    int length = (edge.back()->Length() + ReverseEdge(edge.back())->Length()) / 2;

    while (length < max_length  && (int)edge.size() < max_nodesize) {
        if (edge.back()->OutNode()->InDegree() >= 2 && edge.back()->OutNode()->OutDegree() == 1) {
            
            return { edge, length };

        } else if (edge.back()->OutNode()->InDegree() == 1 && edge.back()->OutNode()->OutDegree() == 1) {
            edge.push_back(edge.back()->OutNode()->GetOutEdges()[0]);
            length += (edge.back()->Length() +  ReverseEdge(edge.back())->Length() ) / 2;
        } else {
            break;
        }
    } 
    return LinearPath();
}





bool StringGraph::IsDiffType(const BaseNode* n0, const BaseNode* n1) {
    
    auto rvs = asmdata_.GetReadVariants();
    auto rr = rvs->Test(n0->ReadId(), n1->ReadId());
    return rr[1] > rr[0];
}


void StringGraph::IdentifySimplePaths() {
    std::unordered_set<BaseEdge*> visited;

    for (auto &i : edges_) {
        BaseEdge* e = i.second; // short name
        if (!e->IsReduce() && visited.find(e) == visited.end()) {
            paths_.push_back(ExtendSimplePath(e, visited));
            auto vpath = Reverse(paths_.back());
            for (auto e : vpath) visited.insert(e);
            paths_.push_back(vpath);
        }
    }


//    assert(Assert_PathDual(paths_));  TOO SLOW
}

std::list<BaseEdge*> StringGraph::ExtendSimplePath(BaseEdge *e, std::unordered_set<BaseEdge*> &visited) {
    assert(visited.find(e) == visited.end());

    std::list<BaseEdge*> path;

    std::unordered_set<BaseNode*> rnodes;

    visited.insert(e);
    path.push_back(e);
    rnodes.insert(ReverseNode(e->InNode()));
    rnodes.insert(ReverseNode(e->OutNode()));
    rnodes.insert(e->InNode());
    rnodes.insert(e->OutNode());

    BaseEdge* curr = path.back();
    while (curr->OutNode()->InDegree() == 1 && curr->OutNode()->OutDegree() == 1 && 
        visited.find(curr->OutNode()->GetOutEdges().front()) == visited.end() &&
        rnodes.find(curr->OutNode()->GetOutEdges().front()->OutNode()) == rnodes.end()) {

        path.push_back(curr->OutNode()->GetOutEdges().front());
        visited.insert(curr->OutNode()->GetOutEdges().front());
        rnodes.insert(ReverseNode(curr->OutNode()->GetOutEdges().front()->OutNode()));
        rnodes.insert(curr->OutNode()->GetOutEdges().front()->OutNode());
        curr = curr->OutNode()->GetOutEdges().front();
    }

    curr = path.front();
    while (curr->InNode()->InDegree() == 1 && curr->InNode()->OutDegree() == 1 &&
        visited.find(curr->InNode()->GetInEdges().front()) == visited.end() &&
        rnodes.find(curr->InNode()->GetInEdges().front()->InNode()) == rnodes.end()) {

        path.push_front(curr->InNode()->GetInEdges().front());
        visited.insert(curr->InNode()->GetInEdges().front());
        rnodes.insert(ReverseNode(curr->InNode()->GetInEdges().front()->InNode()));
        rnodes.insert(curr->InNode()->GetInEdges().front()->InNode());

        curr = curr->InNode()->GetInEdges().front();
    }

    return path;
}


template<typename T, typename U>
bool BelongTo(const T& src, const U& dst) {
    for (auto i : src) {
        if (dst.find(i) == dst.end()) return false;
    }
    return true;
}


std::vector<BaseEdge*> StringGraph::ShortestPath(const BaseNode* src, const BaseNode *dst, std::unordered_set<BaseEdge*> doable, int(*score)(BaseEdge*)) {

    typedef std::tuple<const BaseNode*, BaseEdge*, int> WorkType;
    std::vector<WorkType> nodes;
    std::unordered_map<const BaseNode*, WorkType> done;
    WorkType srcdst{ dst, nullptr, 0 };       

    auto cmp = [](const WorkType &a, const WorkType &b) { return std::get<2>(a) < std::get<2>(b); };

    std::make_heap(nodes.begin(), nodes.end(), cmp);

    nodes.push_back(std::make_tuple(src, (BaseEdge*)nullptr, 0));
    std::push_heap(nodes.begin(), nodes.end());

    while (nodes.size() > 0) {
        WorkType i = nodes.front();
        std::pop_heap(nodes.begin(), nodes.end());
        nodes.pop_back();

        if (std::get<0>(i) == src && src == dst) {
            if ((std::get<1>(srcdst) == nullptr && std::get<1>(i) != nullptr) ||
                (std::get<1>(srcdst) != nullptr && std::get<1>(i) != nullptr && std::get<2>(i) < std::get<2>(srcdst)))
                
                srcdst = i;
        }

        if (done.find(std::get<0>(i)) == done.end()) {
            done[std::get<0>(i)] = i;
            
            for (auto e : std::get<0>(i)->GetOutEdges()) {
                if (doable.find(e) != doable.end()) {
                    nodes.push_back(std::make_tuple(e->OutNode(), e, std::get<2>(i) + score(e)));
                    std::push_heap(nodes.begin(), nodes.end());
                }
            }
        }
    }

    std::vector<BaseEdge*> path;

    auto r = done.end();
    if (src == dst && std::get<1>(srcdst) != nullptr) {
        path.push_back(std::get<1>(srcdst));
        r = done.find(std::get<1>(srcdst)->InNode());
    }
    else {
        r = done.find(dst);
    }

    if (r != done.end() && r->first != src) {
        while (r->first != src) {
            path.push_back(std::get<1>(r->second));
            r = done.find(std::get<1>(r->second)->InNode());
        }
    }

    std::reverse(path.begin(), path.end());

    return path;
   
}



std::list<BaseEdge*> StringGraph::Reverse(const std::list<BaseEdge*>& path) {
    std::list<BaseEdge*> vpath;
    for (auto &e : path) {
        vpath.push_front(ReverseEdge(e));
    }
    return vpath;
    
}

bool StringGraph::Assert_PathDual(const std::list<std::list<BaseEdge*>> paths) {
    auto reverse_path = [&](const std::list<BaseEdge*> &p) {
        std::list<BaseEdge*> vp;
        for (auto &e : p) {
            vp.push_front(ReverseEdge(e));
        }
        return vp;
    };

    auto is_equal = [](const std::list<BaseEdge*> &a, const std::list<BaseEdge*> &b) {
        if (a.size() != b.size()) return false;

        auto ia = a.begin();
        auto ib = b.begin();
        for (; ia != a.end(); ++ia, ++ib) {
            if (*ia != *ib) return false;
        }
        return true;
    };
    
    std::unordered_set<const std::list<BaseEdge*>*> done;

    for (auto &p1 : paths) {
        if (done.find(&p1) == done.end()) {
            done.insert(&p1);

            auto vp1 = reverse_path(p1);
            bool found = false;
            for (auto &p2 : paths) {
                if (done.find(&p2) == done.end()) {
                    if (is_equal(p2, vp1)) {
                        done.insert(&p2);
                        found = true;
                        break;
                    }
                }
            }
            if (!found) return false;
        }
    }
    return true;
}


void StringGraph::SaveEdges(const std::string &fname) {
    auto flush_oss = [](GzFileWriter &writer, std::ostringstream &oss) {
        writer << oss.str();
        oss.str("");
    };

    GzFileWriter writer(fname);
    if (writer.Valid()) {
        std::ostringstream oss;
        oss << std::setprecision(3);
        for (auto i : edges_) {
            BaseEdge* e = i.second;
            auto tile = e->GetTile();
            oss << std::fixed 
                << std::setw(14) << IdToString(e->InNode()->Id()) << " "
                << std::setw(14) << IdToString(e->OutNode()->Id()) << " "
                << tile.start << " " << tile.end << " "
                << e->Score() << " " << e->Identity() << " "
                << e->GetReduceTypeName() << "\n";
            if (oss.tellp() >= 100000000) {
                flush_oss(writer, oss);
            }
        }
        flush_oss(writer, oss);
 


    } else {
        LOG(WARNING)("Failed to open file: %s", fname.c_str());
    }

}


double StringGraph::GetOverlapQuality(const Overlap &ol) {
    const auto &rd_store = asmdata_.GetReadStore();

    const auto & query = rd_store.GetSeq(ol.a_.id);
    const auto & target = rd_store.GetSeq(ol.b_.id);

    auto tseq = target.ToUInt8(ol.b_.start, ol.b_.end);
    auto qseq = query.ToUInt8(ol.a_.start, ol.a_.end, !ol.SameDirect());

    auto r = edlibAlign((const char*)&qseq[0], qseq.size(), (const char*)&tseq[0], tseq.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    if (r.status == EDLIB_STATUS_OK) {
        return 1.0 - r.editDistance * 1.0 / ol.AlignedLength();
    } else {
        return 0.0;
    }
}


void StringGraph::ReduceOtherEdges(const std::unordered_set<BaseEdge*> reserved, BaseEdge::ReduceType type) {
    for (auto &i: edges_) { 
        auto e = i.second;
        if (!e->IsReduce()) {
            auto re = ReverseEdge(e);
            if (reserved.find(e) == reserved.end() && reserved.find(re) == reserved.end()) {
                e->Reduce(type);
                re->Reduce(type);
            }
        }
    }
}



PathGraph::PathGraph(AsmDataset &asmdata) : SgGraph(asmdata) 
     , nodes_(reinterpret_cast<std::unordered_map<SgNode::ID, PathNode*, SgNode::ID::Hash> &>(org_nodes_))
     , edges_(reinterpret_cast<std::unordered_map<PathEdge::ID, PathEdge*, PathEdge::ID::Hash> &>(org_edges_)) {

    simplifiers_["spur"].reset(new Spur2Simplifier(*this));
    simplifiers_["duplicate"].reset(new DuplicateSimplifier(*this));
    simplifiers_["bubble"].reset(new BubbleSimplifier(*this));
    simplifiers_["bridge"].reset(new PathBridgeSimplifier(*this));
    simplifiers_["semi"].reset(new SemiBubbleSimplifier(*this));
    simplifiers_["loop"].reset(new LoopSimplifier(*this));
    simplifiers_["cross"].reset(new CrossSimplifier(*this));

}

PathGraph::~PathGraph() {
}


void PathGraph::BuildFrom(StringGraph& sg) {
    string_graph_ = &sg;
    for (auto &path : sg.GetPaths()) {
        AddEdge(path);
    }
}

void PathGraph::AddEdge(std::list<BaseEdge*>& path) {
    assert(path.size() > 0);

    auto in = nodes_.find(PathNode::CreateId(path.front()->InNode()->Id()));
    if (in == nodes_.end()) {
        auto node = new PathNode(path.front()->InNode());
        auto r = nodes_.insert(std::make_pair(node->Id(), node));
        // TODO
        assert(r.second);
        in = r.first;
    }

    auto out = nodes_.find(PathNode::CreateId(path.back()->OutNode()->Id()));
    if (out == nodes_.end()) {
        auto node = new PathNode(path.back()->OutNode());
        auto r = nodes_.insert(std::make_pair(node->Id(), node));
        // TODO
        assert(r.second);
        out = r.first;
    }

    InsertEdge(new SimplePathEdge(in->second, out->second, path));
}

bool TestOutExtend(PathEdge* e, int minlen, int minnode) {
    int len = e->Length();
    int node = e->NodeSize();
    const PathNode* curr = e->OutNode();

    while ((len < minlen || node < minnode) && curr->InDegree() == 1 && curr->OutDegree() == 1) {
        len += curr->GetOutEdge(0)->Length();
        node += curr->GetOutEdge(0)->NodeSize(); 
        curr = curr->GetOutEdge(0)->OutNode();
        
    }
    return len >= minlen && node >= minnode;
}


std::vector<PathEdge*> PathGraph::ShortestPath(const PathNode* src, const PathNode *dst,
    std::unordered_set<PathNode*> candnodes, int(*score)(PathEdge*)) {

    typedef std::tuple<const PathNode*, PathEdge*, int> WorkType;
    std::vector<WorkType> nodes;
    std::unordered_map<const SgNode*, WorkType> done;

    auto cmp = [](const WorkType &a, const WorkType &b) { return std::get<2>(a) < std::get<2>(b); };
    std::make_heap(nodes.begin(), nodes.end(), cmp);

    nodes.push_back(std::make_tuple(src, (PathEdge*)nullptr, 0));
    std::push_heap(nodes.begin(), nodes.end());

    while (nodes.size() > 0) {
        WorkType i = nodes.front();
        std::pop_heap(nodes.begin(), nodes.end());
        nodes.pop_back();

        if (done.find(std::get<0>(i)) == done.end()) {
            done[std::get<0>(i)] = i;
        
            for (auto e : std::get<0>(i)->out_edges_) {
                if (candnodes.find(e->OutNode()) != candnodes.end()) {
                    nodes.push_back(std::make_tuple(e->OutNode(), e, std::get<2>(i) + score(e)));
                    std::push_heap(nodes.begin(), nodes.end());
                }
            }
        }
    }

    std::vector<PathEdge*> path;


    auto r = done.find(dst);
    if (r != done.end() && r->first !=src) {
        while (r->first != src) {
            path.push_back(std::get<1>(r->second));
            r = done.find(std::get<1>(r->second)->InNode());
        }
    }

    std::reverse(path.begin(), path.end());

    return path;

}

PathGraph::LinearPath PathGraph::FindBridgePath(PathEdge* start, int max_length, int max_nodesize) {
    
    std::vector<PathEdge*> edge = {start};
    int length = (edge.back()->Length() + ReverseEdge(edge.back())->Length()) / 2;
    int nodesize = edge.back()->NodeSize();
        

    while (length < max_length  && nodesize < max_nodesize) {
        if (edge.back()->OutNode()->InDegree() >= 2 && edge.back()->OutNode()->OutDegree() == 1) {
  
            return { edge, length, nodesize };

        } else if (edge.back()->OutNode()->InDegree() && edge.back()->OutNode()->OutDegree() == 1) {
            edge.push_back(edge.back()->OutNode()->OutEdge<PathEdge>(0));
            length += (edge.back()->Length() +  ReverseEdge(edge.back())->Length() ) / 2;
            nodesize += edge.back()->NodeSize();

        } else {
            break;
        }
    } 
    return LinearPath();
}

PathGraph::LinearPath PathGraph::FindBridgePath1(PathEdge* start, int max_length, int max_nodesize) {
    
    std::vector<PathEdge*> edge = {start};
    int length = (edge.back()->Length() + ReverseEdge(edge.back())->Length()) / 2;
    int nodesize = edge.back()->NodeSize();
        

    while (length < max_length  && nodesize < max_nodesize) {
        if (edge.back()->OutNode()->InDegree() == 1 && edge.back()->OutNode()->OutDegree() == 1) {
            edge.push_back(edge.back()->OutNode()->OutEdge<PathEdge>(0));
            length += (edge.back()->Length() +  ReverseEdge(edge.back())->Length() ) / 2;
            nodesize += edge.back()->NodeSize();

        } else {
            return { edge, length, nodesize };
        }
    } 
    return LinearPath();
}

bool PathGraph::HasBridgeJunction(const LinearPath& path, int max_depth) {
    std::unordered_set<SgNode*> locals;
    auto start = path.path[0];
    for (auto ie : start->InNode()->out_edges_) {
        if (ie == start) continue;
        std::list<SgNode*> nodes = GetEgoNodes(ie->OutNode(), max_depth, path.length*3, path.nodesize*3);
        locals.insert(nodes.begin(), nodes.end());
    }

    return locals.find(path.path.back()->OutNode()) == locals.end();
}


void PathGraph::IdentifyPaths(const std::string &method) {

    std::vector<std::list<PathEdge*>> paths;
    std::unordered_set<PathEdge*> visited;

    for (auto &i : edges_) {
        std::deque<PathEdge*> stack;
        stack.push_back(i.second);

        size_t start = paths.size();

        
        while (!stack.empty()) {
            PathEdge* e = stack.front();
            stack.pop_front();

            if (!e->IsReduced() && visited.find(e) == visited.end()) {

                paths.push_back(ExtendPath(e, visited, method));

                for (auto ie : paths.back().back()->OutNode()->out_edges_) {
                    if (!ie->IsReduced() && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
                for (auto ie : paths.back().back()->OutNode()->in_edges_) {
                    if (!ie->IsReduced() && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
                for (auto ie : paths.back().front()->InNode()->out_edges_) {
                    if (!ie->IsReduced() && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
                for (auto ie : paths.back().front()->InNode()->in_edges_) {
                    if (!ie->IsReduced() && visited.find(ie) == visited.end()) {
                        stack.push_back(ie);
                    }
                }
            }

        }

        if (start < paths.size()) {
            clusters_.push_back(Cluster());
            for (size_t i = start; i < paths.size(); ++i) {
                for (auto e : paths[i]) {
                    clusters_.back().edges.insert(e);
                    edge2cluster_[e] = clusters_.size() - 1;
                }
            }
        }
    }

    LOG(INFO)("Contig cluster size: %zd, %zd", clusters_.size(), paths.size());

    for (auto & c : clusters_) c.FindLongest();

    LOG(INFO)("Test whether the paths overlap");

    visited.clear();
    for (auto & path : paths) {
        bool overlapped = false;
        for (auto p : path) {
            if (p->IsReduced() || visited.find(p) != visited.end() || 
                ReverseEdge(p)->IsReduced() || visited.find(ReverseEdge(p)) != visited.end()) {
                overlapped = true;
                break;
            }
        }

        assert(!overlapped);
        if (!overlapped) {
            paths_.push_back(path);
            std::list<PathEdge*> rpath;
            for (auto p : path) {
                visited.insert(p);
                PathEdge* rp = ReverseEdge(p);
                assert(rp != nullptr);
                visited.insert(rp);
                rpath.push_back(rp);
            }
            std::reverse(rpath.begin(), rpath.end());
            paths_.push_back(rpath);
        }
    }

    SortPaths();
}

template<typename TI, typename TO>
std::list<PathEdge*> PathGraph::ExtendPathWithMethod(PathEdge* e, std::unordered_set<PathEdge*> &visited, TI get_in_edge, TO get_out_edge) {
    assert(visited.find(e) == visited.end());

    std::list<PathEdge*> path;
    std::unordered_set<PathNode*> rnodes;

    visited.insert(e);
    visited.insert(ReverseEdge(e));
    path.push_back(e);
    rnodes.insert(ReverseNode(e->InNode()));
    rnodes.insert(ReverseNode(e->OutNode()));

    PathEdge* next = get_out_edge(path.back(), visited);
    while (next != nullptr && rnodes.find(next->OutNode()) == rnodes.end()) {
        assert(visited.find(next) == visited.end());
        path.push_back(next);
        visited.insert(next);
        visited.insert(ReverseEdge(next));
        rnodes.insert(ReverseNode(next->OutNode()));
        next = get_out_edge(next, visited);
    }


    PathEdge* prev = get_in_edge(path.front(), visited);
    while (prev != nullptr && rnodes.find(prev->InNode()) == rnodes.end()) {
        assert(visited.find(prev) == visited.end());

        path.push_front(prev);
        visited.insert(prev);
        visited.insert(ReverseEdge(prev));
        rnodes.insert(ReverseNode(prev->InNode()));
        prev = get_in_edge(prev, visited);
    }
 
    return path;
}


std::list<PathEdge*> PathGraph::ExtendPath(PathEdge *e, std::unordered_set<PathEdge*> &visited, const std::string &method) {
    assert(visited.find(e) == visited.end());

    if (method == "no") {

        auto get_in_edge = [](const PathEdge *e, const std::unordered_set<PathEdge*> &visited) {
            if (e->InNode()->InDegree() == 1 && e->InNode()->OutDegree() == 1) {
                if (visited.find(e->InNode()->in_edges_.front()) == visited.end()) {
                    return e->InNode()->in_edges_.front();
                }
            }
            return (PathEdge*)nullptr;
        };

        auto get_out_edge = [](const PathEdge *e, const std::unordered_set<PathEdge*> &visited) {
            if (e->OutNode()->InDegree() == 1 && e->OutNode()->OutDegree() == 1) {
                if (visited.find(e->OutNode()->out_edges_.front()) == visited.end()) {
                    return e->OutNode()->out_edges_.front();
                }
            }
            return (PathEdge*)nullptr;
        };

        return ExtendPathWithMethod(e, visited, get_in_edge, get_out_edge);
    }
    else if (method == "best") {

        // auto get_in_edge = [](const PathEdge* e, const std::unordered_set<PathEdge*> &visited) {
        //     if (e->InNode()->InDegree() == 1) {
        //         PathEdge* best_out = e->InNode()->GetBestOutEdge();
        //         if (best_out == e) {
        //             if (visited.find(e->InNode()->in_edges_.front()) == visited.end())
        //                 return e->InNode()->in_edges_.front();
        //         }
        //     }
        //     return (PathEdge*)nullptr;
        // };
        // auto get_out_edge = [](const PathEdge* e, const std::unordered_set<PathEdge*> &visited) -> PathEdge* {
        //     if (e->OutNode()->OutDegree() == 1) {
        //         PathEdge* best_in = e->OutNode()->GetBestInEdge();
        //         if (best_in == e) {
        //             if (visited.find(e->OutNode()->out_edges_.front()) == visited.end())
        //                 return e->OutNode()->out_edges_.front();
        //         }
        //     }
        //     return (PathEdge*)nullptr;
        // };
        
        auto get_in_edge = [](const PathEdge* e, const std::unordered_set<PathEdge*> &visited) -> PathEdge* {
            PathEdge* best_out = e->InNode()->GetBestOutEdge();
            if (best_out != e) {
                return nullptr;
            }

            PathEdge* best_in = e->InNode()->GetBestInEdge();
            if (visited.find(best_in) == visited.end()) {
                return best_in;
            } else {
                return nullptr;
            }
        };
        auto get_out_edge = [](const PathEdge* e, const std::unordered_set<PathEdge*> &visited) -> PathEdge* {
            PathEdge* best_in = e->OutNode()->GetBestInEdge();
            if (best_in != e) {
                return nullptr;
            }

            PathEdge* best_out = e->OutNode()->GetBestOutEdge();
            if (visited.find(best_out) == visited.end()) {
                return best_out;
            } else {
                return nullptr;
            }

        };

        return ExtendPathWithMethod(e, visited, get_in_edge, get_out_edge);
    }
    else {
        LOG(ERROR)("Unknow --select-branch = %s", method.c_str());
        return std::list<PathEdge*>{e};
    }


}


void PathGraph::SaveEdges(const std::string &fname) {
    GzFileWriter writer(fname);
    if (writer.Valid()) {
        std::ostringstream oss;

        for (auto &i : edges_) {
            auto e = i.second;
            oss << e->Id().ToString(asmdata_.GetStringPool()) << " " 
                << e->type_ << " " << e->Length() << " " 
                << e->Score() << " " << e->ToString(asmdata_.GetStringPool())  << "\n";

            if (oss.tellp() > 1000000) {
                writer.Flush(oss);
            }
        }
        writer.Flush(oss);

    } else {
        LOG(WARNING)("Failed to open file: %s", fname.c_str());
    }
}

void PathGraph::SortPaths() {
    assert(paths_.size() % 2 == 0);
    std::vector<std::array<size_t,2>> weights(paths_.size() / 2);
    for (size_t i = 0; i < paths_.size(); i += 2) {
        weights[i/2] = { i, std::accumulate(paths_[i].begin(), paths_[i].end(), (size_t)0, [](size_t a, PathEdge* b) { return a + b->Length(); }) };
    }
    std::sort(weights.begin(), weights.end(), [](const std::array<size_t,2>& a, const std::array<size_t,2> &b) {
        return a[1] > b[1];
    });
    
    std::vector<std::list<PathEdge*>> newpaths (paths_.size());
    for (size_t i = 0; i < weights.size(); ++i) {
        std::swap(paths_[weights[i][0]], newpaths[i*2]);
        std::swap(paths_[weights[i][0]+1], newpaths[i*2+1]);
    }
    std::swap(newpaths, paths_);
}

void PathGraph::SavePaths(const std::string &fname) {
    GzFileWriter writer(fname);
    if (writer.Valid()) {
        std::ostringstream oss;

        assert(paths_.size() % 2 == 0);

        for (size_t i = 0; i < paths_.size(); ++i) {
            for (auto e : paths_[i]) {
                oss << "ctg" << i/2 << " " <<  ((i % 2) == 0 ? "F ":"R ")
                    << e->InNode()->Id().ToString(asmdata_.GetStringPool()) 
                    << " -> " 
                    << e->OutNode()->Id().ToString(asmdata_.GetStringPool()) << "\n";
            }

            if (oss.tellp() > 1000000) {
                writer.Flush(oss);
            }
        }

        writer.Flush(oss);

    } else {
        LOG(WARNING)("Failed to open file: %s", fname.c_str());
    }
}

size_t PathGraph::PathLength(const std::list<PathEdge*> &path) {
    size_t len =  std::accumulate(path.begin(), path.end(), 0, [](int a, PathEdge* b) { return a + b->score_; });
    return len;
}


void PathGraph::TestGraph() {
    for (auto &i : nodes_) {
        auto n = i.second;
        for (auto e : n->org_out_edges_) {
            assert(e->InNode() == n);
        }
        for (auto e : n->org_in_edges_) {
            assert(e->OutNode() == n);
        }
    }
}

void PathGraph::InsertLoopNode(LoopNode* n) {

    std::vector<AAAEdge*> new_edges;

    for (size_t i = 0; i < n->OriginInDegree(); ++i) {
        auto e = QueryEdge(n->OriginInEdge(i)->Id());
        
        AAAEdge* ne = nullptr;
        if (n->IsConsistOf(e->InNode())) {
            LOG(INFO)("XXXXXXXX");
            ne = new AAAEdge(static_cast<PathEdge*>(e), n, n);
        } else {
            ne = new AAAEdge(static_cast<PathEdge*>(e), nullptr, n);
        }
        new_edges.push_back(ne);
    }
    
    for (size_t i = 0; i < n->OriginOutDegree(); ++i) {
        auto e = QueryEdge(n->OriginOutEdge(i)->Id());
        assert(e != nullptr);
        
        AAAEdge* ne = nullptr;
        if (n->IsConsistOf(e->OutNode())) {
            //ne = new AAAEdge(static_cast<PathEdge*>(e), n, n);
        } else {
            ne = new AAAEdge(static_cast<PathEdge*>(e), n, nullptr);
            new_edges.push_back(ne);
        }
    }
    

    for (auto e : new_edges) {
        InsertEdge(e);
    }
    
    InsertNode(n);
    n->ReduceSub();
}


void PathGraph::InsertCrossNode(CrossNode* n) {


    std::vector<AAAEdge*> new_edges;
    for (size_t i = 0; i < n->OriginInDegree(); ++i) {
        auto e = QueryEdge(n->OriginInEdge(i)->Id());
        assert(e != nullptr);

        AAAEdge* ne = nullptr;
        if (n->IsConsistOf(e->InNode())) {
            ne = new AAAEdge(static_cast<PathEdge*>(e), n, n);
        } else {
            ne = new AAAEdge(static_cast<PathEdge*>(e), nullptr, n);
        }
        new_edges.push_back(ne);
    }
    
    for (size_t i = 0; i < n->OriginOutDegree(); ++i) {
        auto e = QueryEdge(n->OriginOutEdge(i)->Id());
        assert(e != nullptr);

        AAAEdge* ne = nullptr;
        if (n->IsConsistOf(e->OutNode())) {
            //ne = new AAAEdge(static_cast<PathEdge*>(e), n, n);
        } else {
            ne = new AAAEdge(static_cast<PathEdge*>(e), n, nullptr);
            new_edges.push_back(ne);
        }
    }

    for (auto e : new_edges) {
        InsertEdge(e);
    }
    
    InsertNode(n);
    n->ReduceSub();
}
bool IsClearBubble(const PathNode* start, const PathNode* end, const std::list<PathEdge*> &edges, const std::unordered_set<SgNode*> valids) {
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

BubbleEdge* PathGraph::FindBubble(PathNode* start_node, bool check, int depth_cutoff, int width_cutoff) {
    int length_cutoff = Options().max_bubble_length;

    PathNode* end_node = nullptr;

    std::list<PathEdge*> bundle_edges;      // 瀛樻斁姘旀场鐨勮竟
    std::unordered_map<SgNode*, std::pair<int, int>> visited; // length, score

    std::list<SgNode*> local_node_list = GetEgoNodes(start_node, depth_cutoff);
    std::unordered_set<SgNode*> local_nodes(local_node_list.begin(), local_node_list.end());
    std::unordered_set<PathNode*> tips;

    visited[start_node] = std::make_pair(0, 0);
    for (auto e : start_node->out_edges_) {
        tips.insert(e->OutNode());
        bundle_edges.push_back(e);
    }

    int depth = 0;
    double width = 1.0;
    int length = 0;

    bool loop_detect = false;
    bool meet_error = false;
    bool spur = false;

    do {
        std::unordered_map<PathNode*, PathEdge*> new_visited;     // 鏈€鏂拌￠?璁块棶鑺傜偣锛屽欢鍚庡姞鍏￥visited
        std::unordered_set<PathNode*> newtips, oldtips;     // 鏂颁骇鐢熺殑鏈?姊㈣妭鐐瑰拰鏈a澶勭悊鐨勬湯姊㈣妭鐐
        
        for (auto n : tips) {
            //if (n->out_edges_.size() == 0) continue;        // dead end

            PathEdge *best_in_edge = nullptr;
            for (auto e : n->in_edges_) {
                // 妫€鏌ュ叆鑺傜偣锛屽垎鎴愪笁绫伙細涓嶅湪灞€閮ㄩ泦鍚堜腑銆佸凡缁忚?块棶銆佹病鏈夎?块棶
                // 濡傛灉鎵€鏈夊叆鑺傜偣閮藉凡缁忚?块棶锛屽垯鎵惧嚭鍒嗘暟鏈€楂樼殑杈广€傚苟涓斿彲浠ユ墿灞曞畠鐨勫嚭鑺傜偣
                // 鍚﹀垯鏀硅妭鐐瑰欢鍚庡¤勭悊

                if (local_nodes.find(e->InNode()) != local_nodes.end()) {

                    if (visited.find(e->InNode()) != visited.end()) {
                        if (best_in_edge == nullptr || best_in_edge->score_ < e->score_) {
                            best_in_edge = e;
                        }
                    }
                    else {
                        best_in_edge = nullptr;     // 
                        break;
                    }
                }
                else {
                    // 蹇界暐杩欎釜鍏ヨ妭鐐
                }
            }

            if (best_in_edge != nullptr) {

                assert(n == best_in_edge->OutNode());
                new_visited[n] = best_in_edge;

                // 濡傛灉姘旀场娌℃湁鏀舵暃锛岀户缁-娣诲姞鏂扮殑鏈?姊㈣妭鐐
                if (tips.size() > 1) {
                    for (auto e : n->out_edges_) {
                        if (visited.find(e->OutNode()) != visited.end() ||
                            new_visited.find(e->OutNode()) != new_visited.end()) {
                            loop_detect = true;
                            break;
                        }

                        PathNode *revese_node = ReverseNode(e->OutNode());
                        if (local_nodes.find(e->OutNode()) != local_nodes.end() && 
                            visited.find(revese_node) == visited.end() &&
                            new_visited.find(revese_node) == new_visited.end() ) {

                            if (tips.find(e->OutNode()) == tips.end()) {
                                newtips.insert(e->OutNode());
                            }
                            bundle_edges.push_back(e);
                        }
                        else {
                            meet_error = true;
                            break;
                        }
                    }

                    if (n->out_edges_.size() == 0) {
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
                visited[i.second->InNode()].first + i.second->length_,
                visited[i.second->InNode()].second + i.second->score_);

            // 鏇存柊褰撳墠闀垮害
            if (length < visited[i.first].first) {
                length = visited[i.first].first;
            }

        }

        depth += 1;
        width = 1.0 * bundle_edges.size() / depth;


        if (tips.size() >= 1) {
            tips.clear();
            tips.insert(newtips.begin(), newtips.end());
            tips.insert(oldtips.begin(), oldtips.end());
        }

    } while (tips.size() >= 1 && tips.size() < 6 && !loop_detect && !meet_error && !spur && depth <= depth_cutoff && length <= length_cutoff && (depth <= 10 || width <= width_cutoff));

    if (end_node != nullptr && !loop_detect && !meet_error && !spur && depth <= depth_cutoff && length <= length_cutoff && (depth <= 10 || width <= width_cutoff) && (!check || check && IsClearBubble(start_node, end_node,bundle_edges, local_nodes))) {
        
        return new BubbleEdge(start_node, end_node, bundle_edges, visited[end_node].first, width, visited[end_node].second);
    }
    else {
        return nullptr;
    }
}

std::vector<PathEdge*> PathGraph::Cluster::GetStarts() const {

    std::vector<PathEdge*> ends;
    for (auto e : edges) {        
        bool is_end = true;
        for (size_t i = 0; i < e->InNode()->InDegree(); ++i) {
            if (edges.find(e->InNode()->InEdge<PathEdge>(i)) != edges.end()) {
                is_end = false;
                break;
            }
        }

        if (is_end) {
            ends.push_back(e);
        } 
    }

    if (ends.size() == 0) {
        ends.push_back(*edges.begin());
    }

    return ends;
}

std::vector<PathEdge*> PathGraph::Cluster::GetEnds() const {

    std::vector<PathEdge*> ends;
    for (auto e : edges) {        
        bool is_end = true;
        for (size_t i = 0; i < e->InNode()->InDegree(); ++i) {
            if (edges.find(e->InNode()->InEdge<PathEdge>(i)) != edges.end()) {
                is_end = false;
                break;
            }
        }

        if (is_end) {
            ends.push_back(e);
        } 
    }

    if (ends.size() == 0) {
        ends.push_back(*edges.begin());
    }

    return ends;
}

void PathGraph::Cluster::FindLongest() {

    size_t maxlen = 0;

    std::vector<PathEdge*> starts = GetStarts();
    DUMPER["asm"]("FFF: Starts: %zd, %zd", starts.size(), edges.size());
    for (auto s : starts) {
        std::unordered_set<PathEdge*> visited;
        std::unordered_map<PathEdge*, LPath> longests;

        FindLongest0(s, visited, longests);
        auto lg = longests.find(s);
        assert(lg != longests.end());

        DUMPER["asm"]("FFFe: %zd", lg->second.second.size());
        for (auto e : lg->second.second) {
            nontrivial.insert(e);
        }
        maxlen = std::max(maxlen, lg->second.first);
    }

    std::vector<PathEdge*> ends = GetEnds();
    DUMPER["asm"]("ends: %zd, %zd", ends.size(), edges.size());
    for (auto s : ends) {
        std::unordered_set<PathEdge*> visited;
        std::unordered_map<PathEdge*, LPath> longests;
        FindLongest1(s, visited, longests);
        auto lg = longests.find(s);
        assert(lg != longests.end());

        DUMPER["asm"]("FFFe: %zd", lg->second.second.size());
        for (auto e : lg->second.second) {
            nontrivial.insert(e);
        }
        maxlen = std::max(maxlen, lg->second.first);
    }

    length = maxlen;
}

void PathGraph::Cluster::FindLongest0(PathEdge* e, std::unordered_set<PathEdge*> &visited, std::unordered_map<PathEdge*, LPath>& longests) {

    assert(visited.find(e) == visited.end());
    visited.insert(e);
    DUMPER["asm"]("FFF0 visited %zd", visited.size());
    
    std::unordered_map<PathEdge*, LPath>::const_iterator best = longests.end();
    for (size_t i = 0; i < e->OutNode()->OutDegree(); ++i) {
        auto oe = e->OutNode()->OutEdge<PathEdge>(i);
        if (visited.find(oe) == visited.end()) {
            auto iter = longests.find(oe);
            if (iter == longests.end()) {
                FindLongest0(oe, visited, longests);
                iter = longests.find(oe);
                assert(iter != longests.end());
            }
            if (best == longests.end() || best->second.first < iter->second.first) {
                best = iter;
            }
        }
    }

    if (best == longests.end()) {
        longests[e] = std::make_pair(e->Length(), std::vector<PathEdge*>({e}));
    } else {
        LPath& lp = longests[e];
        lp.first = e->Length() + best->second.first;
        lp.second = best->second.second;
        lp.second.push_back(e);
    }
    visited.erase(e);
}

void PathGraph::Cluster::FindLongest1(PathEdge* e, std::unordered_set<PathEdge*> &visited, std::unordered_map<PathEdge*, LPath>& longests) {

    assert(visited.find(e) == visited.end());

    visited.insert(e);
    DUMPER["asm"]("FFF1 visited %zd", visited.size());
    std::unordered_map<PathEdge*, LPath>::const_iterator best = longests.end();
    for (size_t i = 0; i < e->InNode()->InDegree(); ++i) {
        auto ie = e->InNode()->InEdge<PathEdge>(i);
        if (visited.find(ie) == visited.end()) {
            auto iter = longests.find(ie);
            if (iter == longests.end()) {
                FindLongest1(ie, visited, longests);
                iter = longests.find(ie);
                assert(iter != longests.end());
            }
            if (best == longests.end() || best->second.first < iter->second.first) {
                best = iter;
            }
        }
    }

    if (best == longests.end()) {
        longests[e] = std::make_pair(e->Length(), std::vector<PathEdge*>({e}));
    } else {
        LPath& lp = longests[e];
        lp.first = e->Length() + best->second.first;
        lp.second = best->second.second;
        lp.second.push_back(e);
    }
    visited.erase(e);
}

} // namespace fsa {
    
