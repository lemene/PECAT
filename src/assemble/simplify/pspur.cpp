#include "pspur.hpp"

namespace fsa {


bool Spur2Simplifier::ParseParameters(const std::vector<std::string> &params) { 
    //assert(params[0] == "bridge");

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

void Spur2Simplifier::Running() {
    auto cands = graph_.CollectEdges([](SgEdge* e) {
        return !static_cast<PathEdge*>(e)->IsReduced() && e->InNode()->InDegree() == 0;
    });
    
    std::unordered_set<SgEdge*> removed;

    for (auto e : cands) {
        assert(e->InNode()->InDegree() == 0);
        Debug("conds: %s\n", e->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str());
        
        std::vector<SgEdge*> path { e };
        size_t length = (e->Length() + graph_.SgGraph::ReverseEdge(e)->Length()) / 2;
        size_t nodesize = e->NodeSize();
    
        Debug("conds size0: %zd, %zd\n", length, nodesize);

        while (path.back()->OutNode()->InDegree() == 1 && path.back()->OutNode()->OutDegree() == 1 && 
               length < max_length && nodesize < max_nodesize) {
            path.push_back(path.back()->OutNode()->OutEdge(0));
            auto ee = path.back();
            auto re = graph_.SgGraph::ReverseEdge(path.back());
            length += (path.back()->Length() + graph_.SgGraph::ReverseEdge(path.back())->Length()) / 2;
            nodesize += path.back()->NodeSize();
        }

        Debug("conds size: %zd, %zd\n", length, nodesize);
        if (path.back()->OutNode()->InDegree() >= 2 && path.back()->OutNode()->OutDegree() >= 1 &&
            length < max_length && nodesize < max_nodesize) {


            bool is_spur = true;
            for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
                auto e = path.back()->OutNode()->InEdge(i);
                
                if (e != path.back()) {
                    if (!TestInExtend(e, length*1.001, nodesize)) {
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
        
        Debug("end: %s\n", e->ToString(graph_.GetAsmData().GetStringPool()).c_str());
    }

    Debug("remove: %zd\n", removed.size());
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