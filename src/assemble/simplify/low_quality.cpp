#include "low_quality.hpp"

namespace fsa {


bool LowQualitySimplifier::ParseParameters(const std::vector<std::string> &params) {
    //assert(params[0] == "quality");

    CreateDebugFile();
    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "s") {
        } else {
            return false;
        }
    }
    return true;
}

void LowQualitySimplifier::Running() {
    assert(PreCondition());

    auto threshold = ComputeThreshold();
    RemoveLowQuality(threshold);
    ReactiveEdges(threshold);

}

double LowQualitySimplifier::ComputeThreshold() {


    nodes_ = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() >= 1;
    });

    quals_.resize(nodes_.size());
    for (size_t i = 0; i < nodes_.size(); ++i) {
        assert(nodes_[i]->OutDegree() >= 1);
        quals_[i].assign(nodes_[i]->OutDegree(), 0.0);
    }

    std::atomic<size_t> index { 0 };
    MultiThreadRun((size_t)graph_.Options().thread_size, [&index, this](size_t id) {
        for (size_t i = index.fetch_add(1); i < quals_.size(); i = index.fetch_add(1)) {
            auto n = nodes_[i];
            auto& qs = quals_[i];
            for (size_t i = 0; i < qs.size(); i++) {
                auto e = n->GetOutEdge(i);
                qs[i] = graph_.GetOverlapQuality(*e->ol_);
                Debug("low-quality: %s %s %f\n", graph_.GetAsmData().QueryNameById(e->ol_->a_.id).c_str(), 
                    graph_.GetAsmData().QueryNameById(e->ol_->b_.id).c_str(), qs[i]);
            }
        }
    });

    std::vector<double> qs;
    for (const auto& c : quals_) {
        std::vector<double> qv = c;
        std::sort(qv.begin(), qv.end(), std::greater<double>());
        for (size_t i=0; i < std::min<size_t>(4, qv.size()); ++i) {
            qs.push_back(qv[i]);
        }
    }

    auto m = ComputeMedianAbsoluteDeviation(qs);
    auto auto_threshold = m[0] - 6*1.4826*m[1];
    auto threshold = graph_.Options().min_identity < 0 ? auto_threshold : graph_.Options().min_identity;
    LOG(INFO)("Low-quality threshold %.02f%(%.02f%)", threshold*100, auto_threshold*100);
    return threshold;
}

void LowQualitySimplifier::RemoveLowQuality(double threshold) {

    std::unordered_set<BaseEdge*> removed;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        auto n = nodes_[i];
        const auto& qs = quals_[i];
        for (size_t i = 0; i < n->OutDegree(); ++i) {
            if (qs[i] < threshold) {
                removed.insert(n->GetOutEdge(i));
            }
        }
    }

    LOG(INFO)("Remove low-quality edges: %zd", removed.size()*2);
    graph_.ReduceEdges(removed, BaseEdge::RT_LOW_QUALITY);

    RepairRmoved(threshold, removed);
}


void LowQualitySimplifier::RepairRmoved(double threshold, std::unordered_set<BaseEdge*> &removed) {
    auto is_transitive = [this](BaseEdge *e0, BaseEdge *e1) {
        return graph_.GetEdge(e0->OutNode()->Id(), e1->OutNode()->Id()) != nullptr ||
               graph_.GetEdge(e1->OutNode()->Id(), e0->OutNode()->Id()) != nullptr ;

    };
    std::unordered_set<BaseEdge*> recovery;
    for (auto e0 : removed) {
        auto n = e0->InNode();
        if (n->OutDegree() == 0) {
            std::vector<BaseEdge*> works;
            std::vector<BaseEdge*> www;
            for (auto e1 : n->GetReducedOutEdge()) {
                if (e1->reduce_type_ == BaseEdge::RT_TRANSITIVE) {
                    if (is_transitive(e0, e1)) {
                        works.push_back(e1);
                    }
                }
                www.push_back(e1);
            }

            std::sort(works.begin(), works.end(), [](const BaseEdge* a, const BaseEdge* b) {
                return a->ol_->AlignedLength() > b->ol_->AlignedLength();
            });

            std::vector<BaseEdge*> rcv(e0->InNode()->GetOutEdges());
            for (size_t i = 0; i < works.size(); ++i) {
                auto qual = graph_.GetOverlapQuality(*works[i]->ol_);
                Debug("low-quality repair: %s %s %f\n", graph_.GetAsmData().QueryNameById(works[i]->ol_->a_.id).c_str(), 
                    graph_.GetAsmData().QueryNameById(works[i]->ol_->b_.id).c_str(), qual);
                if (qual >= threshold) {
                    bool done = false;
                    auto a1 = works[i];
                    for (auto a0 : rcv) {
                        if (is_transitive(a0, a1)) {
                            done = true;
                            break;
                        }
                    }
                    if (!done) rcv.push_back(works[i]);
                }

            }

            for (size_t i = e0->InNode()->OutDegree(); i < rcv.size(); ++i)  {
                recovery.insert(rcv[i]);

            }
        }
    }
    LOG(INFO)("Recovery size0: %zd", recovery.size()*2);
    graph_.ReactiveEdges(recovery);
}

void LowQualitySimplifier::ReactiveEdges(double threshold) {

    std::unordered_set<BaseEdge*> recovery;
    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        return n->OutDegree() == 0;
    });

    for (auto n : nodes) {
        if (n->OutDegree() == 0) {
            std::vector<BaseEdge*> works;
            for (auto e : n->GetReducedOutEdge()) {
                if (e->reduce_type_ == BaseEdge::RT_TRANSITIVE) {
                    works.push_back(e);
                }
            }

            std::sort(works.begin(), works.end(), [](const BaseEdge* a, const BaseEdge* b) {
                return a->ol_->AlignedLength() > b->ol_->AlignedLength();
            });

            for (size_t i = 0; i < works.size(); ++i) {
                auto qual = graph_.GetOverlapQuality(*works[i]->ol_);
                Debug("low-quality reactive: %s %s %f\n", graph_.GetAsmData().QueryNameById(works[i]->ol_->a_.id).c_str(), 
                    graph_.GetAsmData().QueryNameById(works[i]->ol_->b_.id).c_str(), qual);
                if (qual >= threshold) {
                    recovery.insert(works[i]);
                    break;
                }

            }
        }
    }

    LOG(INFO)("Recovery size1: %zd", recovery.size()*2);
    graph_.ReactiveEdges(recovery);
}

} // namespace fsa