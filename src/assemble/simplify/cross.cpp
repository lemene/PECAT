#include "cross.hpp"

namespace fsa {


bool CrossSimplifier::PreCondition() {
    return graph_.GetAsmData().GetInconsistentOverlaps() != nullptr;
}

void CrossSimplifier::Running() {
    auto cands = graph_.CollectNodes([this](SgNode* n) {
        return IsCross(n);
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

    for (auto n : crosses) {
        if (n == nullptr) continue; 
        if (graph_.QueryNode(n->Id()) == nullptr && n->IsRaw()) {
            printf("Insert cross node: %s %d\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(), n);
  
            auto rn = n->Reverse(graph_);
            graph_.InsertCrossNode(n);
            graph_.InsertCrossNode(rn);

        } else {
            if (!n->IsRaw()) {
                printf("CrossNode IsRaw 0\n");
            }
            delete n;
        }
    }

}

bool CrossSimplifier::IsCross(const SgNode* in_0) const {
    if (in_0->OutDegree() == 2) {
        const SgNode* out_0 = in_0->OutNode(0);
        const SgNode* out_1 = in_0->OutNode(1);
        if (out_0 != out_1 && out_0->InDegree() == 2 && out_1->InDegree() == 2) {
            const SgNode* in_1 = nullptr;
            for (size_t i = 0; i < out_0->InDegree(); ++i) {
                if (out_0->InNode(i) != in_0) {
                    in_1 = out_0->InNode(i);
                    break;
                }
            }
            if (in_1 != nullptr && in_1->OutDegree() == 2) {
                if ((in_1->OutNode(0) == out_0 &&  in_1->OutNode(1) == out_1) ||
                    (in_1->OutNode(0) == out_1 &&  in_1->OutNode(1) == out_0)) {

                    return TestLength(in_0, in_1) && TestExtends(in_0, in_1, out_0, out_1);
                }
                
            }
        }
    }
    return false;
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

    std::vector<SgNode*> bs = { n->OutNode(0), n->OutNode(1) };
    assert(bs[0]->InDegree() == 2 && bs[1]->InDegree() == 2);
    std::vector<SgNode*> as = { bs[0]->InNode(0), bs[0]->InNode(1) };
    assert(as[0]->OutDegree() == 2 && as[1]->OutDegree() == 2);

    if (IsInconsistent(as[0], as[1]) && IsInconsistent(bs[0], bs[1])) {
        return new CrossNode(as, bs);
    }

    return nullptr;
}

bool CrossSimplifier::TestLength(const SgNode* n0, const SgNode* n1) const {
    for (size_t i = 0; i < n0->OutDegree(); ++i) {
        auto e = n0->OutEdge(i);
        if (e->NodeSize() > max_node_size) {
            return false;
        }
    }

    for (size_t i = 0; i < n1->OutDegree(); ++i) {
        auto e = n1->OutEdge(i);
        if (e->NodeSize() > max_node_size) {
            return false;
        }
    }

    return true;
}

bool CrossSimplifier::TestExtends(const SgNode* in0, const SgNode* in1, const SgNode* out0, const SgNode* out1) const {
    if (in0->InDegree() != 1 || in1->InDegree() != 1 || out0->OutDegree() != 1 || out1->OutDegree() != 1) return false;

    const int NODESIZE = 30;
    if (in0->InEdge(0)->NodeSize() <= NODESIZE && (in0->InNode(0) == out0 || in0->InNode(0) == out1)) return false;
    if (in1->InEdge(0)->NodeSize() <= NODESIZE && (in1->InNode(0) == out0 || in1->InNode(0) == out1)) return false;
    
    if (out0->OutEdge(0)->NodeSize() <= NODESIZE && (out0->OutNode(0) == in0 || out0->OutNode(0) == in1)) return false;
    if (out1->OutEdge(0)->NodeSize() <= NODESIZE && (out1->OutNode(0) == in0 || out1->OutNode(0) == in1)) return false;

    auto test_in_nodesize = [](const SgNode* n) {
        assert(n->OutDegree() == 2 && n->InDegree() == 1);
        size_t nodesize = std::max<size_t>(n->OutEdge(0)->NodeSize(), n->OutEdge(1)->NodeSize());
        return n->InEdge(0)->NodeSize() >= nodesize;
    };

    auto test_out_nodesize = [](const SgNode* n) {
        assert(n->InDegree() == 2 && n->OutDegree() == 1);
        size_t nodesize = std::max<size_t>(n->InEdge(0)->NodeSize(), n->InEdge(1)->NodeSize());
        return n->OutEdge(0)->NodeSize() >= nodesize;
    };

    if (!test_in_nodesize(in0)) return false;
    if (!test_in_nodesize(in1)) return false;

    if (!test_out_nodesize(out0)) return false;
    if (!test_out_nodesize(out1)) return false;

    return true;

}

} // namespace fsa