#include "cross.hpp"

namespace fsa {


bool CrossSimplifier::PreCondition() {
    return graph_.GetAsmData().GetInconsistentOverlaps() != nullptr;
}

void CrossSimplifier::Running() {
    auto cands = graph_.CollectNodes([this](SgNode* n) {
        bool r = IsCross(n);
        Debug("Test cands: %s(IsCross=%d)\n", ToString(n).c_str(), r);
        return r;
    });

    std::vector<CrossNode*> crosses(cands.size(), nullptr);

    LOG(INFO)("Detecting cross structures: %zd %zd", cands.size(), crosses.size());

    MultiThreadMap(1, cands, crosses, [&](SgNode* n) {
        return DetectCross(n);
    });

    LOG(INFO)("Found cross structures %zd", std::count_if(crosses.begin(), crosses.end(), [](const CrossNode* n) {
        return n != nullptr;
    }));

    // 有重复
    std::sort(crosses.begin(), crosses.end(), [](const CrossNode* a, const CrossNode* b) {
        if (a != nullptr && b != nullptr) {
            return *a < *b;
        } else {
            return (size_t)a < (size_t)b;
        }
    });

    std::unordered_set<SgEdge*> removed;
    for (auto n : crosses) {
        if (n == nullptr) continue; 
        bool sol = TrySoluteCrossNode(n, removed);
        if (!sol && graph_.QueryNode(n->Id()) == nullptr && n->IsRaw()) {
            Debug("Insert cross node: %s %d\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(), n);
  
            auto rn = n->Reverse(graph_);
            graph_.InsertCrossNode(n);
            graph_.InsertCrossNode(rn);

        } else {
            if (!n->IsRaw()) {
                Debug("CrossNode IsRaw 0\n");
            }
            delete n;
        }
    }

    
	LOG(INFO)("resolved cross: %zd", removed.size());
	for (auto r : removed) {
		Debug("remove: %s\n", ToString(r).c_str());
		auto pe = static_cast<PathEdge*>(const_cast<SgEdge*>(r));
		if (!pe->IsReduced()) {
			pe->Reduce("resolved_crossnode", true);
			graph_.ReverseEdge(pe)->Reduce("resolved_crossnode", true);
		}
	}

}


bool CrossSimplifier::IsLinkedReversedNode(const SgNode* front, const SgNode* back, size_t max_depth) {

    // collect reversed nodes
    auto curr = graph_.SgGraph::ReverseNode(front);
    std::unordered_set<const SgNode*> rnodes { curr };
    for (size_t i = 0; i < max_depth; ++i) {
        for (size_t i = 0; i < curr->OutDegree(); ++i) {
            rnodes.insert(curr->OutNode(i));
        }
        
        if (curr->OutDegree() == 1) {
            curr = curr->OutNode(0);
        } else {
            break;
        }
    }
    
    // check 
    curr = back;
    if (rnodes.find(curr) != rnodes.end()) {
        return true;
    }
    for (size_t i = 0; i < max_depth; ++i) {
        for (size_t i = 0; i < curr->OutDegree(); ++i) {
            if (rnodes.find(curr->OutNode(i)) != rnodes.end()) {
                return true;
            }
        }
        if (curr->OutDegree() == 1) {
            curr = curr->OutNode(0);
        } else {
            break;
        }
    }

    return false;
}


bool CrossSimplifier::TrySoluteCrossNode(CrossNode* node, std::unordered_set<SgEdge*>& removed) {
    Debug("Solve cross: %s %zd %zd\n", ToString(node).c_str(), node->OriginInDegree(), node->OriginOutDegree());
    if (node->OriginInDegree() != 2 || node->OriginOutDegree() != 2) {
        return false;
    }
    Debug("Solve cross: %s\n", ToString(node).c_str());
    auto in0 = node->OriginInEdge(0)->OutNode();
    auto in1 = node->OriginInEdge(1)->OutNode();
    auto out0 = node->OriginOutEdge(0)->InNode();
    auto out1 = node->OriginOutEdge(1)->InNode();

    auto lnk00 = IsLinkedReversedNode(in0, out0, 5);
    auto lnk01 = IsLinkedReversedNode(in0, out1, 5);
    auto lnk10 = IsLinkedReversedNode(in1, out0, 5);
    auto lnk11 = IsLinkedReversedNode(in1, out1, 5);

    Debug("Solve cross link: %d %d %d %d\n", lnk00, lnk01, lnk10, lnk11);
    if ((lnk00 || lnk11) && !lnk01 && !lnk10) {
        // 01 and 10 is right
        Debug("Solute CrossNode(01): (%s %s) -> (%s %s)\n", ToString(in0).c_str(), ToString(in1).c_str(), ToString(out0).c_str(), ToString(out1).c_str());
        auto path = node->GetInternalPath(0, 0);
        removed.insert(path.begin(), path.end());
        path = node->GetInternalPath(1, 1);
        removed.insert(path.begin(), path.end());
        return true;
    } else if ((lnk01 || lnk10) && !lnk00 && !lnk11) {
        // 00 and 11 is right
        Debug("Solute CrossNode(00): (%s %s) -> (%s %s)\n", ToString(in0).c_str(), ToString(in1).c_str(), ToString(out0).c_str(), ToString(out1).c_str());
        auto path = node->GetInternalPath(1, 0);
        removed.insert(path.begin(), path.end());
        path = node->GetInternalPath(0, 1);
        removed.insert(path.begin(), path.end());
        return true;

    } else {
        return false;
    }
}

bool CrossSimplifier::IsInconsistent(const SgNode* n0, const SgNode* n1) const {
    if (n0->IsType("path") && n1->IsType("path")) {

        const PathNode* pn0 = (const PathNode*)n0;
        const PathNode* pn1 = (const PathNode*)n1;
        if (pn0->GetBaseNode() != nullptr && pn1->GetBaseNode() != nullptr) {
            auto inconsist = graph_.GetAsmData().GetInconsistentOverlaps();

            int id0 = pn0->GetBaseNode()->ReadId();
            int id1 = pn1->GetBaseNode()->ReadId();

            return inconsist->Contain(id0, id1);

        }
    }

    return false;
}

CrossNode* CrossSimplifier::DetectCross(SgNode* n) const {
    assert(IsCross(n));

    auto in_p0 = GetOutLinearPath(n->OutEdge(0));
    Debug("in_p0 %s(%s), %s(%s)\n", 
        ToString(in_p0.front()).c_str(), in_p0.front()->Id().ToString().c_str(), 
        ToString(in_p0.back()).c_str(), in_p0.back()->Id().ToString().c_str());
    auto in_p1 = GetOutLinearPath(n->OutEdge(1));
    Debug("in_p1 %s(%s), %s(%s)\n", 
        ToString(in_p1.front()).c_str(), in_p1.front()->Id().ToString().c_str(), 
        ToString(in_p1.back()).c_str(), in_p1.back()->Id().ToString().c_str());
    std::vector<SgNode*> bs = { in_p0.back()->OutNode(), in_p1.back()->OutNode() };
    assert(bs[0]->InDegree() == 2 && bs[1]->InDegree() == 2);
    
    auto out_p0 = GetInLinearPath(bs[0]->InEdge(0));
    Debug("out_p0 %s(%s), %s(%s)\n", 
        ToString(out_p0.front()).c_str(), out_p0.front()->Id().ToString().c_str(), 
        ToString(out_p0.back()).c_str(), out_p0.back()->Id().ToString().c_str());
    auto out_p1 = GetInLinearPath(bs[0]->InEdge(1));
    Debug("out_p1 %s(%s), %s(%s)\n", 
        ToString(out_p1.front()).c_str(), out_p1.front()->Id().ToString().c_str(), 
        ToString(out_p1.back()).c_str(), out_p1.back()->Id().ToString().c_str());
    std::vector<SgNode*> as = { out_p0.front()->InNode(), out_p1.front()->InNode() };
    assert(as[0]->OutDegree() == 2 && as[1]->OutDegree() == 2);
    Debug("Nodes: %s(%s),%s(%s) -> %s(%s),%s(%s)\n", 
        ToString(as[0]).c_str(), as[0]->Id().ToString().c_str(),
        ToString(as[1]).c_str(), as[1]->Id().ToString().c_str(),
        ToString(bs[0]).c_str(), bs[0]->Id().ToString().c_str(),
        ToString(bs[1]).c_str(), bs[1]->Id().ToString().c_str());
    if (IsInconsistent(as[0], as[1]) || IsInconsistent(bs[0], bs[1])) {
        return new CrossNode(as, bs);
    }

    return nullptr;
}


} // namespace fsa