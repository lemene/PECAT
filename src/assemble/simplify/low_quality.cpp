#include "low_quality.hpp"

namespace fsa {


bool QualitySimplifier::ParseParameters(const std::vector<std::string> &params) {
    //assert(params[0] == "quality");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "s") {
        } else {
            return false;
        }
    }
    return true;
}

void QualitySimplifier::Running() {
    assert(PreCondition());

    auto threshold = ComputeThreshold();
    RemoveLowQuality(threshold);
    ReactiveEdges(threshold);
    
    //ReactiveContainedEdges(threshold);
}

double QualitySimplifier::ComputeThreshold() {


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

    auto m0 = ComputeMeanAbsoluteDeviation(qs);
    auto auto_threshold0 = m0[0] -  6*1.253*m0[1];
    auto m1 = ComputeMedianAbsoluteDeviation(qs);
    auto auto_threshold1 = m1[0] - 6*1.4826*m1[1];
    auto auto_threshold = (auto_threshold0 + auto_threshold1*2) / 3;
    graph_.desired_edge_quality = auto_threshold;
    auto threshold = graph_.Options().min_identity < 0 ? auto_threshold : graph_.Options().min_identity;
    LOG(INFO)("Low-quality threshold %.02f%(%.02f%, %.02f%)", threshold*100, auto_threshold1*100, auto_threshold0*100);
    graph_.actual_min_identity = threshold;
    return threshold;
}

void QualitySimplifier::RemoveLowQuality(double threshold) {

    std::unordered_set<BaseEdge*> removed;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        auto n = nodes_[i];
        const auto& qs = quals_[i];
        for (size_t i = 0; i < n->OutDegree(); ++i) {
            if (qs[i] < threshold) {
                auto e = n->GetOutEdge(i);
                bool has_ok_dup = false;
                auto nol = ReplaceHighQualityOverlap(e->ol_, threshold);
                if (nol.first != nullptr) {
                    graph_.GetAsmData().ReplaceOverlapInGroup(nol.first, e->ol_);
                    e->Replace(nol.first);
                    graph_.ReverseEdge(e)->Replace(nol.first);
                    has_ok_dup = true;

                }
                if (!has_ok_dup) removed.insert(n->GetOutEdge(i));
            }
        }
    }

    LOG(INFO)("Remove low-quality edges: %zd", removed.size()*2);
    graph_.ReduceEdges(removed, BaseEdge::RT_LOW_QUALITY);

    RepairRmoved(threshold, removed);
}





void QualitySimplifier::RepairRmoved(double threshold, std::unordered_set<BaseEdge*> &removed) {
    auto is_transitive = [this](BaseEdge *e0, BaseEdge *e1) {
        return graph_.GetEdge(e0->OutNode()->Id(), e1->OutNode()->Id()) != nullptr ||
               graph_.GetEdge(e1->OutNode()->Id(), e0->OutNode()->Id()) != nullptr ;

    };
    std::unordered_set<BaseEdge*> recovery;
    for (auto e0 : removed) {
        auto n = e0->InNode();
        Debug("repair cand: %s\n", ToString(e0).c_str());
        if (n->OutDegree() == 0) {
            std::vector<BaseEdge*> works;
            for (auto e1 : n->GetReducedOutEdge()) {
                Debug("test: %s\n", ToString(e1).c_str());
                if (e1->reduce_type_ == BaseEdge::RT_TRANSITIVE) {
                    {//if (is_transitive(e0, e1)) {
                        works.push_back(e1);
                        Debug("test add: %s\n", ToString(e1).c_str());
                    }
                }
            }

            std::sort(works.begin(), works.end(), [](const BaseEdge* a, const BaseEdge* b) {
                return a->ol_->AlignedLength() > b->ol_->AlignedLength();
            });

            assert(e0->InNode()->GetOutEdges().size() == 0);
            std::vector<BaseEdge*> rcv(e0->InNode()->GetOutEdges()); 
            for (size_t i = 0; i < works.size(); ++i) {
                auto qual = graph_.GetOverlapQuality(*works[i]->ol_);
                Debug("low-quality repair: %s %s %f\n", graph_.GetAsmData().QueryNameById(works[i]->ol_->a_.id).c_str(), 
                    graph_.GetAsmData().QueryNameById(works[i]->ol_->b_.id).c_str(), qual);
                if (qual >= threshold) {
                    auto done = std::find_if(rcv.begin(), rcv.end(), [&works, i, is_transitive](BaseEdge* a0) { 
                        return is_transitive(a0, works[i]);
                    });
                    if (done == rcv.end()) {
                        rcv.push_back(works[i]);
                    }
                } else {
                    auto nol = ReplaceHighQualityOverlap(works[i]->ol_, threshold);
                    if (nol.first != nullptr) {
                        graph_.GetAsmData().ReplaceOverlapInGroup(nol.first, works[i]->ol_);
                        works[i]->Replace(nol.first);
                        graph_.ReverseEdge(works[i])->Replace(nol.first);
                    
                        auto done = std::find_if(rcv.begin(), rcv.end(), [&works, i, is_transitive](BaseEdge* a0) { 
                            return is_transitive(a0, works[i]);
                        });
                        if (done == rcv.end()) {
                            rcv.push_back(works[i]);
                        }
                    }
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

void QualitySimplifier::ReactiveEdges(double threshold) {

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


void QualitySimplifier::ReactiveContainedEdges(double threshold) {

    std::unordered_set<BaseEdge*> recovery;
    auto nodes = graph_.CollectNodes([](BaseNode* n) {
        if (n->OutDegree() == 0 && n->InDegree()) {
            const int N = 10;
            SgNode* curr = n->InNode(0);
            for (size_t i = 0; i < N; ++i) {
                if (curr->InDegree() != 1 || curr->OutDegree() != 1) {
                    return false;
                }
                curr = curr->InNode(0);
            }
            return true;
        }
        return false;
    });


    auto &asmdata = graph_.GetAsmData();
    auto get_extend = [&asmdata, this, threshold](Seq::EndId start, const std::unordered_set<Seq::Id> &added) {

        Seq::Id sid = Seq::EndIdToId(start);
        int send = Seq::End(start);

        auto ols= asmdata.GetExtendOverlaps(sid, send);
    
        bool linked = false;
        std::pair<Seq::EndId, const Overlap*> result = { Seq::NID, nullptr};
        for (auto ol : ols) {
            auto& alt = ol->GetOtherRead(sid);
            if (added.find(alt.id) != added.end()) continue;

            auto rend = send == 0 ? (ol->SameDirect() ? 0 : 1) : (ol->SameDirect() ? 1 : 0);
            auto it = asmdata.read_infos_.find(alt.id);
            auto n = graph_.GetNode(Seq::IdToEndId(alt.id, rend));

            
            //LOG(INFO)("check: %s : %lld", graph_.GetAsmData().GetStringPool().QueryStringById(alt.id).c_str(), n);

            if (!linked || (linked && n != nullptr && result.second->AlignedLength() < ol->AlignedLength())) {
                if ((n == nullptr && it != asmdata.read_infos_.end() && it->second.filtered.type == RdReason::RS_CONTAINED) ||
                    (n != nullptr && n->OutDegree() > 0 && (n->InDegree() == 0 || (n->InDegree()==1 && n->InNode(0)->OutDegree() > 1)))) { 
                    auto qual = graph_.GetOverlapQuality(*ol);
                    
                    if (qual < threshold) {
                        auto nol = ReplaceHighQualityOverlap(ol, threshold);
                        if (nol.first != nullptr) {
                            qual = nol.second;
                            ol = nol.first;
                        }
                    }

                    if (qual >= threshold) {

                        if ((linked && n != nullptr && result.second->AlignedLength() < ol->AlignedLength()) ||
                            (!linked && (n != nullptr || result.second == nullptr  || result.second->AlignedLength() < ol->AlignedLength()))) {

                                //LOG(INFO)("   check1: %f",  qual);
                                result.first = Seq::IdToEndId(alt.id, rend);
                                result.second = ol;
                                if (n != nullptr) {
                                    linked = true;
                                }
                            }
                    }
                }
            }
        }

        return result;
    };

    std::vector<std::pair<SgEdgeID, std::vector<const Overlap*>>> paths;

    LOG(INFO)("Contained: %zd", nodes.size());
    for (auto n : nodes) {
        assert(n->OutDegree() == 0);
        Debug("reactive contained cand: %s\n", ToString(n).c_str());
            
        std::unordered_set<Seq::Id> added;

        std::vector<std::pair<Seq::EndId, const Overlap*>> path;
        Seq::EndId curr = n->Id().Value(0);
        for (size_t i = 0; i < 15; ++i) {
            added.insert(Seq::EndIdToId(curr));
            auto ext = get_extend(curr, added);

            if (ext.first == Seq::NID) {
                // find no path
                break;
            } else if (graph_.GetNode(ext.first) != nullptr) {
                // find a path
                paths.push_back(std::pair<SgEdgeID, std::vector<const Overlap*>>());
                paths.back().first = SgEdgeID(ID_EDGE_TYPE_MULT, n->Id(), SgNodeID(ext.first));
                for (auto p : path) {
                    paths.back().second.push_back(p.second);
                }
                paths.back().second.push_back(ext.second);
                
                Debug("reactive contained end:%s %zd\n", n->Id().ToString(graph_.GetAsmData().GetStringPool()).c_str(), paths.back().second.size());
                break;

            } else {
                // continue
                path.push_back(ext);
                curr = ext.first;
                Debug("reactive contained ols:%d, %s\n", path.size(), graph_.GetAsmData().GetStringPool().QueryStringById(Seq::EndIdToId(ext.first)).c_str());
                
            }
        }
        
    }

    LOG(INFO)("recover contained: %zd", paths.size());

    std::unordered_set<const Overlap*> repair;
    std::unordered_set<SgEdgeID, SgEdgeID::Hash> done;
    for (const auto& p : paths) {
        Debug("add path: %s\n", p.first.ToString(graph_.GetAsmData().GetStringPool()).c_str());
        if (done.find(p.first) == done.end()) {
            done.insert(p.first);
            done.insert(SgEdgeID::Reverse(p.first));
            
            Debug("add path1: %s\n", p.first.ToString(graph_.GetAsmData().GetStringPool()).c_str());
            repair.insert(p.second.begin(), p.second.end());
        }
    }
    
    LOG(INFO)("recover contained2: %zd", repair.size()*2);
    Repair(repair);
}
    

void QualitySimplifier::Repair(const std::unordered_set<const Overlap*> &ols) {

    auto get_edge = [this](const Overlap* ol) {
        if (ol->SameDirect()) {
            Seq::EndId fB = Seq::IdToEndId(ol->a_.id, 0);
            Seq::EndId gB = Seq::IdToEndId(ol->b_.id, 0);
            return graph_.GetEdge(fB, gB);
        } else {
            Seq::EndId fB = Seq::IdToEndId(ol->a_.id, 0);
	        Seq::EndId gE = Seq::IdToEndId(ol->b_.id, 1);
            return graph_.GetEdge(fB, gE);
        }
    };
    
    for (auto ol : ols) {

        auto e = get_edge(ol);
        if (e != nullptr) {
            graph_.ReactiveEdges(std::vector<BaseEdge*>({e}));

        } else {
            graph_.AddOverlap(ol);
        }

    }
}


} // namespace fsa
