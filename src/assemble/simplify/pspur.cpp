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
            length += (path.back()->Length() + graph_.SgGraph::ReverseEdge(path.back())->Length()) / 2;
            nodesize += path.back()->NodeSize();
        }

        Debug("conds size: %zd, %zd, %zd, %zd\n", length, max_length, nodesize, max_nodesize);
        if (path.back()->OutNode()->InDegree() >= 2 && path.back()->OutNode()->OutDegree() >= 1 &&
            length < max_length && nodesize < max_nodesize) {


            bool is_spur = true;
            for (size_t i = 0; i < path.back()->OutNode()->InDegree(); ++i) {
                auto e = path.back()->OutNode()->InEdge(i);
                
                if (e != path.back()) {
                    if (!TestInExtend(e, length*2, nodesize*2)) {
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
    //RemoveDeadEnds();
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


void Spur2Simplifier::RemoveDeadEnds() {
    auto cands = graph_.CollectNodes([](SgNode* n) {
        return n->InDegree() == 1 && n->OutDegree() > 1;
    });

    auto is_dead_end = [this](SgNode* n, const std::list<SgNode*> &nodes) {
        if (n->InDegree() > 1) return false;

        if (nodes.size() >= max_nodesize) return false;

        std::unordered_set<SgNode*> nodeset(nodes.begin(), nodes.end());
        nodeset.insert(n->InNode(0));

        size_t node_size = 0;
        for (auto inode : nodes) {

            for (size_t i = 0; i < inode->InDegree(); i++) {
                if (nodeset.find(inode->InNode(i)) == nodeset.end()) return false;
            }
            for (size_t i = 0; i < inode->OutDegree(); i++) {
                if (nodeset.find(inode->OutNode(i)) == nodeset.end()) return false;
            }

            for (size_t i = 0; i < inode->OutDegree(); i++) {
                node_size += inode->OutEdge(i)->NodeSize();
            }
        }

        if (node_size < max_nodesize) return false;

        return true;

    };

    std::unordered_set<SgEdge*> removed;

    for (auto node : cands) {
        Debug("cand dead end: %s\n", ToString(node).c_str());
        assert( node->InDegree() == 1 && node->OutDegree() > 1);

        auto in_nodes = graph_.GetEgoNodes(graph_.SgGraph::ReverseNode(node->InNode(0)), max_nodesize, max_length);

        if (is_dead_end(graph_.SgGraph::ReverseNode(node->InNode(0)), in_nodes)) continue;

        std::vector<std::list<SgNode*>> deads;
        for (size_t i = 0; i < node->OutDegree(); ++i) {
            auto out_nodes = graph_.GetEgoNodes(node->OutNode(i), max_nodesize, max_length);
            if (is_dead_end(node->OutNode(i), out_nodes)) {
                deads.push_back(out_nodes);
            }
        }

        if (deads.size() < node->OutDegree()) {
            for (auto nlist : deads) {
                assert(nlist.front()->InDegree() == 1);
                Debug("remove dead end: %s\n", ToString(nlist.front()->InEdge(0)).c_str());
                removed.insert(nlist.front()->InEdge(0));
                for (auto n : nlist) {
                    for (size_t i = 0; i < n->OutDegree(); ++i) {
                        removed.insert(n->OutEdge(i));
                    }

                }
            }
             
        }
    }
    
    Debug("remove: %zd\n", removed.size());
    for (auto e : removed) {
        auto pe = static_cast<PathEdge*>(e);
        if (!pe->IsReduced()) {
            pe->Reduce("spur:3", true);
            graph_.ReverseEdge(pe)->Reduce("spur:3", true);
        }
    }

}

} // namespace fsa