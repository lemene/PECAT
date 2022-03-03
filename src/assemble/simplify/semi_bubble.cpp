#include "semi_bubble.hpp"

namespace fsa {


bool SemiBubbleSimplifier::PreCondition() {
    return graph_.GetAsmData().GetInconsistentOverlaps() != nullptr;
}


bool SemiBubbleSimplifier::ParseParameters(const std::vector<std::string> &params) {
    assert(params[0] == "semi");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (false) {
   
        } else {
            if (!ParseParameter(it)) {
                return false;
            }
        }
    }
    return true;
}


void SemiBubbleSimplifier::Running() {
    assert(graph_.GetAsmData().GetInconsistentOverlaps() != nullptr);

    auto cands = graph_.CollectNodes([](const SgNode* n) {
        return n->OutDegree() > 1;
    });

    std::vector<SemiBubbleEdge*> semi(cands.size(), nullptr);
    MultiThreadMap(1, cands, semi, [this](SgNode* n) -> SemiBubbleEdge* {
        Debug("Find semi node id0: %s\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
        auto e = Detect(static_cast<PathNode*>(n));
        if (e != nullptr) {
            Debug("Find semi id0: %s\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
            if (!e->Validate(graph_.GetAsmData().GetInconsistentOverlaps(), graph_.GetReadStore(), graph_.string_graph_, graph_.Options().max_bubble_length)) { 
                delete e;
                e = nullptr;
            }
        }
        return e;
    });
 
    LOG(INFO)("Found semi-bubble structures %zd", std::count_if(semi.begin(), semi.end(), [](const SgEdge* n) {
        return n != nullptr;
    }));

    for (auto e : semi) {
        if (e == nullptr) continue;
            
        assert(e->Validate(graph_.GetAsmData().GetInconsistentOverlaps(), graph_.GetReadStore(), graph_.string_graph_, graph_.Options().max_bubble_length));

        Debug("Find semi id: %s\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
        auto oe = graph_.QueryEdge(e->Id());
        if (oe == nullptr) {
            auto re = e->Reverse(graph_);
            assert(graph_.QueryEdge(re->Id()) == nullptr);
            graph_.InsertEdge(e);
            graph_.InsertEdge(re);

        } else {
            auto re = e->Reverse(graph_);
            Debug("Find semi type: %s\n", oe->Type().c_str());
            assert(oe->IsType("semi-bubble"));
            static_cast<SemiBubbleEdge*>(oe)->Merge(*e);

            auto roe = graph_.QueryEdge(re->Id());
            Debug("Find semi type r: %s\n", oe->Type().c_str());
            assert(roe != nullptr && roe->IsType("semi-bubble"));
            static_cast<SemiBubbleEdge*>(roe)->Merge(*re);

            delete re;
            delete e;
        }
    }

    graph_.TestGraph();
}

SemiBubbleEdge* SemiBubbleSimplifier::Detect(PathNode* start_node, int depth_cutoff, int width_cutoff, int length_cutoff) {

    assert (start_node->OutDegree() >= 2);
    std::vector< std::vector<PathEdge*>> paths;
    bool err = false;

    auto is_loop_end = [](PathEdge* curr) {
        return curr->InNode()->InDegree() == 2 && curr->OutNode()->OutDegree() == 1 && 
                curr->OutNode()->OutEdge(0)->OutNode() == curr->InNode();
    };

    auto is_edge_valid = [](PathEdge* e) {
        return e->IsType("simple") || e->IsType("alt") || (e->IsType("bubble") && e->MinPathSize() < 9);
    };

    for (size_t i = 0; i < start_node->OutDegree(); ++i) {
        auto e = (PathEdge*)start_node->OutEdge(i);
        Debug("Find semi contine 0: %s, %d\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(), is_edge_valid(e));
        if (!is_edge_valid(e)) continue;

        std::vector<PathEdge*>  path;
        Debug("Find semi contine 0: %s, - %d\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(), e->OutNode()->OutDegree());
        if (e->OutNode()->OutDegree() == 1) {
            path.push_back(e);
            PathEdge* curr = e->OutNode()->OutEdge<PathEdge>(0);
            while (is_edge_valid(curr) && curr->InNode()->InDegree() == 1 && curr->OutNode()->OutDegree() == 1) {

                path.push_back(curr);
                curr = curr->OutNode()->OutEdge<PathEdge>(0);
            }
        Debug("Find semi contine e: %s\n", curr->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
            
            if (is_edge_valid(curr) && curr->InNode()->InDegree() == 1 && curr->OutNode()->OutDegree() == 0 ) {
                path.push_back(curr);
            } else if (is_edge_valid(curr) && curr->InNode()->InDegree() == 2 && curr->OutNode()->OutDegree() == 1 && 
                curr->OutNode()->OutEdge(0)->OutNode() == curr->InNode()) {                      // loop
                path.push_back(curr);
            } 
            paths.push_back(path);
        } else {
            path.push_back(e);
            paths.push_back(path);
            //err = true;
        }
    }
    Debug("Find semi contine x %d: %zd\n", err, paths.size());

    if (!err && paths.size() >= 2) {

        for (size_t i=1; i<paths.size(); ++i) {
            if (paths[i].back()->OutNode() == paths[0].back()->OutNode()) {
                return nullptr;
            }
        }

        std::sort(paths.begin(), paths.end(), [is_loop_end](const std::vector<PathEdge*>& a, const std::vector<PathEdge*> &b) {
            int aend = a.back()->OutNode()->OutDegree() > 0 && !is_loop_end(a.back()) ? 1 : 0;
            int bend = b.back()->OutNode()->OutDegree() > 0 && !is_loop_end(b.back()) ? 1 : 0;
            return aend > bend || (aend == bend && a.size() > b.size());
        });

        Debug("Find semi contine xxxx  %d: %zd, %zd\n", paths[1].back()->OutNode()->OutDegree() == 0 || is_loop_end(paths[1].back()), paths[0].size(), paths[1].size());
        if (paths[1].back()->OutNode()->OutDegree() == 0 || is_loop_end(paths[1].back())) {

            if (paths[0].size() > paths[1].size() + 30) {
                paths[0].erase(paths[0].begin()+paths[1].size(), paths[0].end());
            }
            SemiBubbleEdge *edge = new SemiBubbleEdge(paths);
            return edge;
        }
        return nullptr;
            
    } else {
        return nullptr;
    }
}




} // namespace fsa