#include "cross_node.hpp"
#include "../string_node.hpp"

#include "../string_graph.hpp"

namespace fsa {

CrossNode::CrossNode(const std::vector<SgNode*>& in_nodes, const std::vector<SgNode*>& out_nodes) 
 : PathNode(nullptr, 0), in_nodes_(in_nodes), out_nodes_(out_nodes) {

    assert(in_nodes.size() == 2 && out_nodes.size() == 2);
    assert(in_nodes[0]->OutDegree() == 2 && in_nodes[1]->OutDegree() == 2);
    assert(out_nodes[0]->InDegree() == 2 && out_nodes[1]->InDegree() == 2);

    for (auto n : in_nodes) {
        for (size_t i = 0; i < n->InDegree(); ++i) {

            origin_in_edges_.push_back(n->InEdge(i));
        }
    }
    for (auto n : out_nodes) {
        for (size_t i = 0; i < n->OutDegree(); ++i) {
            origin_out_edges_.push_back(n->OutEdge(i));
        }
    }

    for (auto n : in_nodes) {
        for (size_t i = 0; i < n->OutDegree(); ++i) {
            LinearPath path;
            path.edges = GetOutLinearPath(n->OutEdge(i));
            DUMPER["cross"]("OutEdges: %s, %s\n", path.edges.front()->Id().ToString().c_str()
                                                , path.edges.back()->Id().ToString().c_str());
            origin_edges_.push_back(path);
        }
    }
    DUMPER["cross"]("OutNodes %s / %s, %s / %s\n", 
        origin_edges_[0].OutNode()->Id().ToString().c_str(), 
        origin_edges_[1].OutNode()->Id().ToString().c_str(), 
        origin_edges_[2].OutNode()->Id().ToString().c_str(), 
        origin_edges_[3].OutNode()->Id().ToString().c_str());
    assert(origin_edges_.size() == 4);
    assert(origin_edges_[0].OutNode() == origin_edges_[2].OutNode() || origin_edges_[0].OutNode() == origin_edges_[3].OutNode());
    assert(origin_edges_[1].OutNode() == origin_edges_[2].OutNode() || origin_edges_[1].OutNode() == origin_edges_[3].OutNode());

    id_ = CreateID();
    AttachType("cross");    

}

CrossNode* CrossNode::Reverse(PathGraph& sg) const {
    
    std::vector<SgNode*> ins(out_nodes_.size());
    std::transform(out_nodes_.begin(), out_nodes_.end(), ins.begin(), [&sg](SgNode*n) {
        return sg.ReverseNode((PathNode*)n);
    });
        
    std::vector<SgNode*> outs(in_nodes_.size());
    std::transform(in_nodes_.begin(), in_nodes_.end(), outs.begin(), [&sg](SgNode*n) {
        return sg.ReverseNode((PathNode*)n);
    });
        DUMPER["cross"]("Reverse CrossNode 0: %zd %zd \n", ins.size(), outs.size());
        DUMPER["cross"]("Reverse CrossNode 1: %zd %zd, %zd %zd\n", ins[0]->OutDegree(), ins[1]->OutDegree(), outs[0]->InDegree(), outs[1]->InDegree());
        DUMPER["cross"]("Reverse CrossNode 3: %s %s, %s %s\n", 
            ins[0]->Id().ToString(sg.GetAsmData().GetStringPool()).c_str(), 
            ins[1]->Id().ToString(sg.GetAsmData().GetStringPool()).c_str(), 
            outs[0]->Id().ToString(sg.GetAsmData().GetStringPool()).c_str(), 
            outs[1]->Id().ToString(sg.GetAsmData().GetStringPool()).c_str());
    auto r = new CrossNode(ins, outs);
    DUMPER["cross"]("REVERS %s, %s, %s\n", Id().ToString(sg.GetAsmData().GetStringPool()).c_str(), r->Id().ToString(sg.GetAsmData().GetStringPool()).c_str(), 
        SgNodeID::Reverse(r->Id()).ToString().c_str());
    
    assert(r->Id() == SgNodeID::Reverse(Id()));
    return r;
}

bool CrossNode::IsRaw() const {
    for (auto &p : origin_edges_) {
        for (auto e : p.edges) {
            if (static_cast<PathEdge*>(e)->IsReduced()) {
                return false;
            }
        }
    }
    return true;
}

SgNodeID CrossNode::CreateID() const {
    std::vector<int> values;
    
    values.insert(values.end(), in_nodes_[0]->Id().Values().begin(), in_nodes_[0]->Id().Values().end());
    values.insert(values.end(), in_nodes_[1]->Id().Values().begin(), in_nodes_[1]->Id().Values().end());
    values.insert(values.end(), out_nodes_[0]->Id().Values().begin(), out_nodes_[0]->Id().Values().end());
    values.insert(values.end(), out_nodes_[1]->Id().Values().begin(), out_nodes_[1]->Id().Values().end());

    std::sort(values.begin(), values.end());

    return SgNodeID(ID_NODE_TYPE_CROSS, values);
}


void CrossNode::ReduceSub() {
    for (auto &p : origin_edges_) {
        for (auto e : p.edges) {
            static_cast<PathEdge*>(e)->Reduce("Contain");
        }
    }
}

SgNode* CrossNode::LinearPath::OutNode() const { return edges.back()->OutNode(); }
SgNode* CrossNode::LinearPath::InNode()  const { return edges.front()->InNode(); }


std::list<SgEdge*> CrossNode::GetOutLinearPath(SgEdge* e) {
    std::list<SgEdge*> edges {e};
    std::unordered_set<SgEdge*> done {e};

    while (edges.back()->OutNode()->OutDegree() == 1 && edges.back()->OutNode()->InDegree() == 1) {
        auto oe = edges.back()->OutNode()->OutEdge(0);
        if (done.find(oe) == done.end()) {
            edges.push_back(oe);
            done.insert(oe);
        } else {
            break;
        }
    }
    return edges;
}

std::list<SgEdge*> CrossNode::GetInLinearPath(SgEdge* e) {
    std::list<SgEdge*> edges {e};
    std::unordered_set<SgEdge*> done {e};

    while (edges.front()->InNode()->InDegree() == 1 && edges.front()->InNode()->OutDegree() == 1) {
        auto ie = edges.front()->InNode()->InEdge(0);
        if (done.find(ie) == done.end()) {
            edges.push_front(ie);
            done.insert(ie);
        } else {
            break;
        }
    }
    return edges;
}

}