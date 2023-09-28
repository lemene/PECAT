#include "best.hpp"

namespace fsa {


bool BestOverlapsSimplifier::ParseParameters(const std::vector<std::string> &params) {
    assert(params[0] == "best");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "cmp") {
            auto ps = SplitStringByChar(it[1], ',');
            if (ps.size() >= 1 && ps[0].size() > 0) {
                th_count = std::stoi(ps[0]);
            }
            if (ps.size() >= 2 && ps[1].size() > 0) {
                th_rate[0] = std::stod(ps[1]);
            }
            if (ps.size() >= 3 && ps[2].size() > 0) {
                th_rate[1] = std::stod(ps[2]);
            }
        } else {
            if (!ParseParameter(it)) {
                return false;
            }
        }
    }
    return true;
}


std::string BestOverlapsSimplifier::GetParameters() const {
    std::ostringstream oss;
    oss << "cmp=" << th_count << "," << th_rate[0] << "," << th_rate[1];
    return oss.str();
}


auto BestOverlapsSimplifier::FindBestOutEdge(BaseNode* n) -> std::vector<BestItem> {
    std::vector<BestItem> bests;
    for (size_t i = 0; i < n->OutDegree(); ++i) {
        auto e = (BaseEdge*)n->OutEdge(i);
        bool done = false;
        auto s = GetEdgeScore(e);
        for (size_t i = 0; i < bests.size(); ++i) {

            if (!IsOutEdgeInconsistent(e, bests[i].e, 3, bad_ols_)) {
                if (EdgeScoreGreater(s, bests[i].s)) {
                    bests[i] = {e, s};
                }         
                done = true;
            }        
        }
        if (!done) {
            bests.push_back({e, s});
        }
    }

    // remove duplications
    std::vector<BestItem> bests_nodup;
    std::unordered_set<BaseEdge*> added;
    for (auto &b : bests) {
        if (added.find(b.e) == added.end()) {
            added.insert(b.e);
            bests_nodup.push_back(b);
        }
    }
    Debug("begin sort %zd\n", bests.size());
    for (auto b : bests) {
        Debug("sorting: e=%s s=(%d %d %d %d)\n", ToString(b.e).c_str(), b.s[0], b.s[1], b.s[2], b.s[3]);
    }
    std::sort(bests_nodup.begin(), bests_nodup.end(), [this](const BestItem &a, const BestItem &b) {
        return EdgeScoreGreater(a.s, b.s);
    });
    //assert(bests.size() == bests_nodup.size());
    Debug("end sort %zd\n", bests_nodup.size());
    for (auto b : bests_nodup) {
        Debug("sorted: e=%s s=(%d %d %d %d)\n", ToString(b.e).c_str(), b.s[0], b.s[1], b.s[2], b.s[3]);
    }


    return bests_nodup;
}


void BestOverlapsSimplifier::FindCandidateBestEdges() {
    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() > 0;
    });

    cands_.resize(nodes.size());

    std::atomic<size_t> index { 0 };

    auto work_func = [&](size_t tid) {
        for (size_t i = index.fetch_add(1); i < nodes.size(); i = index.fetch_add(1)) {
            auto n = nodes[i];
            assert(n->OutDegree() > 0);
            cands_[i][0] = FindBestOutEdge(n);;
        }
    };

    MultiThreadRun((size_t)graph_.Options().thread_size, work_func);
}


void BestOverlapsSimplifier::ComfirmCandidateBestEdges() {

    std::unordered_map<BaseEdge*, int> counts;
    std::unordered_set<BaseEdge*> abandoned;
    std::unordered_set<BaseNode*> ok_out_nodes;
    
    for (auto &i : cands_) {
        auto &out_bests = i[0];

        if (out_bests.size() > 0) {

            for (size_t i = 0; i < out_bests.size(); ++i) {
                auto e = out_bests[i].e;
                auto s = out_bests[i].s;

                Debug("cand best type: %s -> %s, %d %d %d %d, (%zd/%zd)\n", graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
                    graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str(), s[0], s[1], s[2], s[3], i, out_bests.size());
 
                if (!EdgeScoreSignificantlyGreater(out_bests[0].s, out_bests[i].s)) {
                    Debug("cand best type ok: %s -> %s, %d %d %d %d, | %zd %d\n", graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
                        graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str(), s[0], s[1], s[2], s[3], out_bests[i].e, counts[out_bests[i].e]);
                    counts[out_bests[i].e] += 1;
                    counts[graph_.ReverseEdge(out_bests[i].e)] += 1;
                    Debug("cand best type ok: %s -> %s, %d %d %d %d, | %zd %d\n", graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
                        graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str(), s[0], s[1], s[2], s[3], out_bests[i].e, counts[out_bests[i].e]);
                } else {
                    Debug("cand discard: %s -> %s\n", graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
                        graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str());
                    auto ve = graph_.ReverseEdge(e);
                    Debug("cand discard: %s -> %s\n", graph_.GetAsmData().QueryNameById(ve->InNode()->ReadId()).c_str(),
                        graph_.GetAsmData().QueryNameById(ve->OutNode()->ReadId()).c_str());
                    abandoned.insert(e);
                    abandoned.insert(graph_.ReverseEdge(e));
                    e->subject_ = false;
                    graph_.ReverseEdge(e)->subject_ = false;
                    //break;
                }
            }
        }

        if (EdgeScoreSignificantlyGreater(out_bests[0].s, EdgeScore({0, 0, 0}))) {
            ok_out_nodes.insert(graph_.ReverseNode(out_bests[0].e->InNode()));
        }
    }

    std::unordered_set<SgNode*> s_out_nodes;
    for (auto &i : counts) {
        Debug("cand best count: %s -> %s, %d\n", graph_.GetAsmData().QueryNameById(i.first->InNode()->ReadId()).c_str(),
            graph_.GetAsmData().QueryNameById(i.first->OutNode()->ReadId()).c_str(), i.second);
        if (i.second >= 2) {
            best_edges.insert(i.first);
            s_out_nodes.insert(i.first->OutNode());
            s_out_nodes.insert(graph_.ReverseNode(i.first->InNode()));
        }
    }

    for (auto &i : cands_) {
        auto &out_bests = i[0];

        if (out_bests.size() > 0) {
            bool has = false;
            for (size_t i = 0; i < out_bests.size(); ++i) {

                if (best_edges.find(out_bests[i].e) != best_edges.end()) {
                    has = true;
                    break;
                }
            }

            if (!has) {
                bool added = false;
                for (size_t i = 0; i < out_bests.size(); ++i) {
                    auto e = out_bests[i].e;
                    if (s_out_nodes.find(e->OutNode()) == s_out_nodes.end()) {
                        Debug("cand add: %s -> %s\n", graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
                            graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str());
                        best_edges.insert(e);
                        
                        s_out_nodes.insert(e->OutNode());
                        s_out_nodes.insert(graph_.ReverseNode(e->InNode()));
                        added = true;
                        break;
                    }
                }

                if (!added) {
                    for (size_t i = 0; i < out_bests.size(); ++i) {
                        auto e = out_bests[i].e;
                        if (true || ok_out_nodes.find(e->OutNode()) == ok_out_nodes.end()) {
                            Debug("cand add: %s -> %s\n", graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
                                graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str());
                            best_edges.insert(e);
                            break;
                        }
                    }
                }
            }
        }
    }
}

void BestOverlapsSimplifier::Running() {
    FindCandidateBestEdges();
    ComfirmCandidateBestEdges();
    LOG(INFO)("best_edges size: %zd", best_edges.size());
    graph_.ReduceOtherEdges(best_edges, BaseEdge::RT_NO_BEST);
}


} // namespace fsa  