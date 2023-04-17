#include "string_node.hpp"

#include "string_graph.hpp"

namespace fsa {

std::string SgNodeID::ToString(const StringPool &sp) const {
    std::ostringstream oss;
    oss << t ;
    for (auto iv : v) {
        oss << "_" << sp.QueryStringById(Seq::EndIdToId(iv));
    }
    return oss.str();
}

std::string SgNodeID::ToString() const {
    std::ostringstream oss;
    oss <<  (int)t;
    for (auto iv : v) {
        oss << "_" << iv;
    }
    return oss.str();
}

SgNode* SgNode::InNode(size_t i) {
    return InEdge(i)->InNode();
}

const SgNode* SgNode::InNode(size_t i) const {
    return InEdge(i)->InNode();
}

SgNode* SgNode::OutNode(size_t i) { 
    return OutEdge(i)->OutNode(); 
}

const SgNode* SgNode::OutNode(size_t i) const { 
    return OutEdge(i)->OutNode(); 
}


std::vector<BaseNode*> BaseNode::GetAllOutNodes() const {
    std::vector<BaseNode*> nodes;
    for (auto e : out_edges_) {
        nodes.push_back(e->OutNode());
    }
    for (auto e : reduced_out_edges_) {
        nodes.push_back(e->OutNode());
    }
    return std::move(nodes);
}

std::vector<BaseNode*> BaseNode::GetAllInNodes() const {
    std::vector<BaseNode*> nodes;
    for (auto e : in_edges_) {
        nodes.push_back(e->InNode());
    }
    for (auto e : reduced_in_edges_) {
        nodes.push_back(e->InNode());
    }
    return std::move(nodes);

}


BaseEdge* BaseNode::GetBestInEdge() const {
    return in_edges_.size() > 0 ? *std::max_element(in_edges_.begin(), in_edges_.end(), BaseEdge::Less) : nullptr;  
}

BaseEdge* BaseNode::GetBestOutEdge() const {
    return out_edges_.size() > 0 ? *std::max_element(out_edges_.begin(), out_edges_.end(), BaseEdge::Less) : nullptr; 
}

PathEdge* PathNode::GetBestInEdge() const {
    for (auto e : in_edges_) {
        e->Contain(string_node_->GetBestInEdge());
        return e;
    }
    return nullptr;
}

PathEdge* PathNode::GetBestOutEdge() const {
    for (auto e : out_edges_) {
        e->Contain(string_node_->GetBestOutEdge());
        return e;
    }
    return nullptr;
}

LoopNode::LoopNode(const std::vector<SgEdge*> &backward) 
 : PathNode(nullptr, 0), backward_(backward) {
    assert(backward_.front()->InNode() == backward.back()->OutNode());
    id_ = CreateID();
    AttachType("loop"); 
}

LoopNode::LoopNode(const std::vector<SgEdge*> &backward, const std::vector<SgEdge*> &forward) 
 : PathNode(nullptr, 0), backward_(backward), forward_(forward) {
     
    if (forward.size() > 0) {
        assert(backward_.front()->InNode() == forward_.back()->OutNode());
        assert(backward_.back()->OutNode() == forward_.front()->InNode());
    } else {
        assert(backward_.front()->InNode() == backward.back()->OutNode());
    }
    id_ = CreateID();
    AttachType("loop"); 


}

LoopNode* LoopNode::Reverse(PathGraph& sg) const {
    LOG(INFO)("ReversePath B %zd", backward_.size());
    auto rbackward = ((SgGraph&)sg).ReversePath(backward_);
    LOG(INFO)("ReversePath F %zd %zd", forward_.size(), rbackward.size());

    auto rforward = ((SgGraph&)sg).ReversePath(forward_);
    LOG(INFO)("NEW LOOPNODE %zd", rforward.size());
    return new LoopNode(rbackward, rforward);
}

size_t LoopNode::OriginOutDegree() const { 
    return backward_.front()->InNode()->OutDegree() - 1; 
}
size_t LoopNode::OriginInDegree() const { 
    return backward_.back()->OutNode()->InDegree() - 1; 
}
const SgEdge* LoopNode::OriginOutEdge(size_t i) const {
    SgNode* n = backward_.front()->InNode();

    size_t true_ie = 0;
    for (size_t ie = 0; ie < n->OutDegree(); ++ie) {
        SgEdge* e =  n->OutEdge(ie);
        if (e != backward_.front()) {
            if (true_ie == i) {
                return e;
            }
            true_ie ++;
        }
    }
    assert(!"never come here");
    return nullptr;
}

const SgEdge* LoopNode::OriginInEdge(size_t i) const {
    SgNode* n = backward_.back()->OutNode();

    size_t true_ie = 0;
    for (size_t ie = 0; ie < n->InDegree(); ++ie) {
        SgEdge* e =  n->InEdge(ie);
        if (e != backward_.back()) {
            if (true_ie == i) {
                return e;
            }
            true_ie ++;
        }
    }
    assert(!"never come here");
    return nullptr;
}

bool LoopNode::IsConsistOf(const SgNode* n) const {
    LOG(INFO)("---- %s = %s", n->Id().ToString().c_str(), backward_.front()->InNode()->Id().ToString().c_str());
    return backward_.front()->InNode() == n;
}

void LoopNode::ReduceSub() {
    assert(forward_.size() == 0);
    for (auto e : backward_) {
        static_cast<PathEdge*>(e)->Reduce("Contain");
    }
}


SgNode::ID LoopNode::CreateID() const {
    SgNode::ID aid = backward_.back()->OutNode()->Id();
    SgNode::ID bid = backward_.front()->InNode()->Id();
    assert(aid == bid);

    return SgNode::ID(ID_NODE_TYPE_LOOP, aid);
}

}