#include "spur.hpp"

namespace fsa {


bool SpurSimplifier::ParseParameters(const std::vector<std::string> &params) {
    assert(params[0] == "spur");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "n") {
            min_nodesize = std::stoi(it[1]);
        } else if (it[0] == "l") {
            min_length = std::stoi(it[1]);
        } else {
            return false;
        }
    }
    return true;
}

std::string SpurSimplifier::GetParameters() const {
    std::ostringstream oss;
    oss << "n=" << min_nodesize << ":l=" << min_length;
    return oss.str();
}

//
//         ->
//        /
// ------>------------->
// 
void SpurSimplifier::Running() {
    std::unordered_set<BaseEdge*> removed;
    std::mutex mutex;
    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() > 1;
    });


    std::atomic<size_t> index { 0 };

    auto combine_func = [&removed, &mutex](std::unordered_set<BaseEdge*> &rs) {
        std::lock_guard<std::mutex> lock(mutex);

        removed.insert(rs.begin(), rs.end());
    };

    auto work_func = [&](size_t tid) {
        std::unordered_set<BaseEdge*> rs;

        for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
            auto n = nodes[i];
            assert(n->OutDegree() > 1);

            size_t del = 0;           
            for (auto e : n->GetOutEdges()) {
                assert(!e->IsReduce());

                if (IsTrivialEdge(e, min_nodesize, min_length)) {
                    rs.insert(e);
                    Debug("spur: %s", ToString(e).c_str());
                    for (auto n = e->OutNode(); n->OutDegree() == 1; n = n->OutNode<BaseNode>(0)) {
                        assert(n->OutDegree() == 1 && n->InDegree() == 1);
                        rs.insert(n->OutEdge<BaseEdge>(0));
                        Debug("-> %s", ToString(n->OutEdge<BaseEdge>(0)).c_str());
                    }
                    Debug("\n");
                    del++;
                }
            
                if (del + 1 == n->OutDegree()) break;
            }
        }
        combine_func(rs);
    };
    

    MultiThreadRun((size_t)graph_.Options().thread_size, work_func);

    for (auto e : removed) {
        if (!e->IsReduce()) {
            e->Reduce(BaseEdge::RT_SPUR);
            graph_.ReverseEdge(e)->Reduce(BaseEdge::RT_SPUR);
        }
    }

}

bool SpurSimplifier::IsTrivialEdge(const SgEdge* e, int min_nodesize, int min_length) {
    int nodesize = 1;
    int length = e->Length();

    auto n = e->OutNode();
    while (n->OutDegree() == 1 && n->InDegree() == 1 && (nodesize < min_nodesize || length < min_length)) {
        nodesize += 1;
        length += n->OutEdge(0)->Length();
        n = n->OutNode(0);
    }
    return n->OutDegree() == 0 && n->InDegree() == 1 && (nodesize < min_nodesize || length < min_length);

}
} // namespace fsa