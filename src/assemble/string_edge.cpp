#include "string_edge.hpp"

#include "string_graph.hpp"

namespace fsa {

template<typename T>
static void PathToString(std::ostream &oss, const T& path, const StringPool &sp) {
    oss << path.front()->InNode()->Id().ToString(sp);
    for (auto e : path) {
        oss << "->" << e->OutNode()->Id().ToString(sp);
    }
}

std::string SgEdgeID::ToString(const StringPool &sp) const {
    std::ostringstream oss;
    oss << type_;
    for (auto v : values_) {
        oss << "~" << v.ToString(sp);
    }
    return oss.str();
}

std::string SgEdgeID::ToString() const {
    std::ostringstream oss;
    oss << type_;
    for (auto v : values_) {
        oss << "~" << v.ToString();
    }
    return oss.str();
}

Seq::Tile BaseEdge::GetTile() const {
    
    int start = 0;
    int end = 0;
    if (read_ == ol_->a_.id) {
        if (ol_->a_.start == 0) {
            start = ol_->a_.end;
            end = ol_->a_.len;
        } else {
            start = 0;
            end = ol_->a_.start;
        }
    } else {
        if (ol_->b_.start == 0) {
            start = ol_->b_.end;
            end = ol_->b_.len;
        } else {
            start = 0;
            end = ol_->b_.start;
        }

    }
    return {OutNode()->ReadId(), OutNode()->Id().End() == 1 ? 0 : 1, start, end};
}


void BaseEdge::Reactivate() {
    
    OutNode()->ReactivateInEdge(this);
    InNode()->ReactivateOutEdge(this);
    reduce_type_ = RT_ACTIVE;
}

double BaseEdge::Identity() const {
    assert(ol_ != nullptr);
    return ol_->identity_;
}


size_t BaseEdge::Length() const {
    auto tile = GetTile();
    return (tile.end - tile.start);
}

bool BaseEdge::TestOutExtend(int minlen, int minnode, bool both) {
    int len = Length();
    int node = 1;
    const BaseNode* curr = OutNode();

    while ((len < minlen || node < minnode) && curr->InDegree() == 1 && curr->OutDegree() == 1) {
        len += curr->GetOutEdge(0)->Length();
        curr = curr->GetOutEdge(0)->OutNode();
        node += 1; 
    }
    if (both) return  len >= minlen && node >= minnode;
    else      return  len >= minlen || node >= minnode;

}

bool BaseEdge::TestInExtend(int minlen, int minnode, bool both) {
    int len = Length();
    int node = 1;
    const BaseNode* curr = InNode();

    while ((len < minlen || node < minnode) && curr->InDegree() == 1 && curr->OutDegree() == 1) {
        len += curr->GetInEdge(0)->Length();
        curr = curr->GetInEdge(0)->InNode();
        node += 1; 
    }

    if (both) return  len >= minlen && node >= minnode;
    else      return  len >= minlen || node >= minnode;
}


std::string SimplePathEdge::ToString(const StringPool &sp) const {
    assert(path_.size() > 0);
    std::ostringstream oss;

    PathToString(oss, path_, sp);

    return oss.str();
}

std::string BubbleEdge::ToString(const StringPool &sp) const {
    assert(simple_paths_.size() > 0);

    std::ostringstream oss;
    for (auto& p : simple_paths_) {
        if (oss.str().empty()) {
            oss << p->Id().ToString(sp);
        } else {
            oss << "|" << p->Id().ToString(sp);
        }
    }
    return oss.str();
}


void BubbleEdge::IdentifySimplePaths(StringGraph& string_graph) {
    if (string_edges_.size() == 0) {
        assert(simple_paths_.size() > 0);

        // LOG(INFO)("MaxFlow %s, %s", in_node_->Id().ToString(string_graph.GetReadStore().GetStringPool()).c_str(), 
        //     out_node_->Id().ToString(string_graph.GetReadStore().GetStringPool()).c_str());
        std::unordered_set<const SgEdge*> doable(simple_paths_.begin(), simple_paths_.end());
        // LOG(INFO)("MaxFlow doable size, %zd", doable.size());
        // for (auto e : doable) {
        //     LOG(INFO)("MaxFlow doable %zd", e);
        //     LOG(INFO)("MaxFlow doable %s", e->Id().ToString(string_graph.GetReadStore().GetStringPool()).c_str());
        // }
        
        std::vector<std::vector<SgEdge*>> subpaths = SgGraph::MaximumFlow(in_node_, out_node_, doable);
        // LOG(INFO)("MaxFlow %zd, %zd", subpaths.size(), doable.size());
        for (auto & subpath : subpaths) {
            std::list<BaseEdge*> s;
            for (auto & e : subpath) {
                assert(e->IsType("path"));
                auto pe = static_cast<PathEdge*>(e);
                pe->IdentifySimplePaths(string_graph);
                auto sp0 = pe->GetSimplePath(0);
                s.insert(s.end(), sp0.begin(), sp0.end());
            }
            string_edges_.push_back(s);
        }

        if (string_edges_.size() > 0) {
            std::unordered_map<BaseEdge*, std::array<int,3>> scores;
            for (size_t i = 0; i < string_edges_.size(); ++i) {
                auto &sg = string_graph;
                auto &p = string_edges_[i];
            
                std::array<int,3> score = {0, 0, 0};
                auto rvr = sg.GetAsmData().GetReadVariants();
                if ( rvr != nullptr) {
                    auto s = rvr->Test(*p.front()->ol_);
                    score[0] += -s[1];
                    score[1] += s[0];

                    s = rvr->Test(*p.back()->ol_);
                    score[0] += -s[1];
                    score[1] += s[0];
                }
                score[2] += p.front()->Score();
                score[2] += p.back()->Score();
                scores[p.front()] = score;
            }
            

            std::sort(string_edges_.begin(), string_edges_.end(), [&scores](const std::list<BaseEdge*> &a, const std::list<BaseEdge*> &b) {
                return CompareArray<3,int>(scores[a.front()], scores[b.front()]) > 0;
            });
        }
        if (string_edges_.size() == 0) {
            LOG(INFO)("ASSERT ERROR: %s", Id().ToString().c_str());
        }
        assert(string_edges_.size() > 0);

        length_ = 0;
        for (auto& p : string_edges_) {
            length_ = std::max(length_, std::accumulate(p.begin(), p.end(), 0, [](int a, BaseEdge* b) { return a + b->Length(); }));
        }
    }
}



bool ValidateDipolid(std::list<BaseEdge*> &path0, std::list<BaseEdge*> &path1, PhaseInfoFile &ignored, const ReadStore &rs) {
    int count = 0;
    for (auto p0 : path0) {
        auto name0 = rs.QueryNameById(p0->OutRead()->id);
        for (auto p1 : path1) {
            auto name1 = rs.QueryNameById(p1->OutRead()->id);
            if (ignored.Contain(name0, name1)) {
                count++;
                break;
            }
        }
    }
    // return count*1.0 / std::min(path0.size(), path1.size()) > 0.3;
    return count*1.0 / path1.size() > 0.3;
}

bool BubbleEdge::Validate(PhaseInfoFile *ignored, const ReadStore &rs, StringGraph *sg, int max_length) {
    assert(sg != nullptr );

    IdentifySimplePaths(*sg);

    if (string_edges_.size() >= 2) {
        if (ignored != nullptr) {
            for (size_t i=0; i<string_edges_.size(); ++i) {
                for (size_t j=i+1; j<string_edges_.size(); ++j) {
                    validated_ = ValidateDipolid(string_edges_[i], string_edges_[j], *ignored, rs);
                    if (validated_) {
                        break;
                    }
                }
            }
        }
        if (!validated_) {
            for (const auto &path : string_edges_) {
                int length = std::accumulate(path.begin(), path.end(), 0, [](int a, BaseEdge* e) {
                    return a + e->Length();
                });
                if (length > max_length) {
                    return false;
                }
            }
            return true;
        }
        return true;

    } else {
        return true;
    }

}

std::string SemiBubbleEdge::ToString(const StringPool &sp) const {
    assert(paths_.size() > 0);

    std::ostringstream oss;

    for (const auto &path : paths_) {
        PathToString(oss, path, sp);
        oss << "|";
    }
    return oss.str();
}

void SemiBubbleEdge::ReduceStringEdge() {
    
    for (const auto &path : paths_) {
        for (const auto &p : path) {
            p->ReduceStringEdge();
        }
    }
}

bool SemiBubbleEdge::Contain(const BaseEdge *e) const {
    for (const auto &path : paths_) {
        for (const auto &p : path) {
            if (p->Contain(e)) return true;
        }
    }

    return false;
}

std::unordered_set<Seq::Id> SemiBubbleEdge::GetReads() {
    std::unordered_set<Seq::Id> reads;
    for (const auto &path : paths_) {
        for (const auto &p : path) {
            auto r = p->GetReads();
            reads.insert(r.begin(), r.end());
        }
    }
    return reads;
}

SemiBubbleEdge* SemiBubbleEdge::Reverse(PathGraph &graph) {
    std::vector<std::vector<PathEdge*>> vpaths(paths_.size());

    for (size_t i=0; i<paths_.size(); ++i) {
        for (auto &p : paths_[i]) {
            vpaths[i].push_back(graph.ReverseEdge(p));
        }
        std::reverse(vpaths[i].begin(), vpaths[i].end());
    }

    return new SemiBubbleEdge(vpaths);

}


void SemiBubbleEdge::IdentifySimplePaths(StringGraph& sg) {
    if (string_edges_.size() == 0) {
        std::list<std::list<BaseEdge*>> others;
        for (const auto &path : paths_) {
            std::list<BaseEdge*> edges; 
            for (const auto &p : path) {
                p->IdentifySimplePaths(sg);
                assert(p->SimplePathSize() >= 1);

                auto& es = p->GetSimplePath(0);
                edges.insert(edges.end(), es.begin(), es.end());
                for (size_t i=1; i < p->SimplePathSize(); ++i) {
                    others.push_back(p->GetSimplePath(i));
                }
            }
            string_edges_.push_back(std::move(edges));
        }

        assert(string_edges_.size() >= 2);
    }
}

bool SemiBubbleEdge::Validate(PhaseInfoFile *ignored, const ReadStore &rs, StringGraph *sg, int max_length) {
    assert(sg != nullptr );

    IdentifySimplePaths(*sg);
        
    auto path_len = [](const std::list<BaseEdge*>& path) {
        size_t len = 0;
        for (auto p : path) {
            len += p->Length();
        }
        return len;
    };

    if (string_edges_.size() >= 2 && ignored != nullptr) {
        bool r = true;
        auto len0 = path_len(string_edges_[0]);
        auto b = Branchs();
        bool check_len = ((b&1) != 0 && string_edges_[0].back()->OutNode()->OutDegree() > 0 ) || 
                         ((b&2) != 0 && string_edges_[0].front()->InNode()->InDegree() > 0);
        for (size_t j=1; j<string_edges_.size(); ++j) {
            if ((check_len && path_len(string_edges_[j]) >= len0 * 1.2) || !ValidateDipolid(string_edges_[0], string_edges_[j], *ignored, rs)) {
                r = false;
                break;
            }
        }
        return r;

    } else {
        return false;
    }

    

}

std::string LoopEdge::ToString(const StringPool &sp) const {
    std::ostringstream oss;
    
    PathToString(oss, forward_, sp);
    oss << "|";
    PathToString(oss, backward_, sp);
    oss << "|";

    return oss.str();

}


void LoopEdge::ReduceStringEdge() {
    
    for (const auto &p : forward_) {
        static_cast<PathEdge*>(p)->ReduceStringEdge();
    }
    for (const auto &p : backward_) {
        static_cast<PathEdge*>(p)->ReduceStringEdge();
    }
}

bool LoopEdge::Contain(const BaseEdge *e) const {
    for (const auto &p : forward_) {
        if (static_cast<PathEdge*>(p)->Contain(e)) return true;
    }
    for (const auto &p : backward_) {
        if (static_cast<PathEdge*>(p)->Contain(e)) return true;
    }

    return false;
}

std::unordered_set<Seq::Id> LoopEdge::GetReads() {
    std::unordered_set<Seq::Id> reads;

    for (const auto &p : forward_) {
        auto r = static_cast<PathEdge*>(p)->GetReads();
        reads.insert(r.begin(), r.end());
    }
    for (const auto &p : backward_) {
        auto r = static_cast<PathEdge*>(p)->GetReads();
        reads.insert(r.begin(), r.end());
    }
    return reads;
}

void LoopEdge::IdentifySimplePaths(StringGraph& string_graph) {
    if (string_edges_.size() == 0) {
        string_edges_.push_back(std::list<BaseEdge*> ());
        std::list<BaseEdge*> &edges = string_edges_.back(); 

        std::vector<std::vector<SgEdge*>*> paths { &forward_, &backward_, &forward_} ;
        for (auto path : paths) {
            for (const auto &p : *path) {
                static_cast<PathEdge*>(p)->IdentifySimplePaths(string_graph);
                assert(static_cast<PathEdge*>(p)->SimplePathSize() >= 1);
                auto& es = static_cast<PathEdge*>(p)->GetSimplePath(0);
                edges.insert(edges.end(), es.begin(), es.end());
            }

        }
    }
}

LoopEdge* LoopEdge::Reverse(PathGraph &graph) {
    
    std::vector<SgEdge*> rforward;
    std::vector<SgEdge*> rbackward;

    
    for (auto &p : forward_) {
        rforward.push_back(graph.ReverseEdge(static_cast<PathEdge*>(p)));
        std::reverse(rforward.begin(), rforward.end());

    }

    for (auto &p : backward_) {
        rbackward.push_back(graph.ReverseEdge(static_cast<PathEdge*>(p)));
        std::reverse(rbackward.begin(), rbackward.end());
    }
    
    return new LoopEdge(graph.ReverseNode(OutNode()), graph.ReverseNode(InNode()), rforward, rbackward);

}

void LoopEdge::AddToGraph(PathGraph &path) {

}

void AAAEdge::IdentifySimplePaths(StringGraph& sg) {
    if (string_edges_.size() == 0) {
        origin->IdentifySimplePaths(sg);
        assert(origin->SimplePathSize() >= 1);
        string_edges_.assign(1, std::list<BaseEdge*>());

        if (in_node_->IsType("loop") && in_node_->InDegree() == 0) {
            auto n = static_cast<LoopNode*>(in_node_);
            auto back = n->GetBackward();
            for (auto e : back) {
                if (e->IsType("path")) {
                    static_cast<PathEdge*>(e)->IdentifySimplePaths(sg);
                    auto p = static_cast<PathEdge*>(e)->GetSimplePath(0);
                    string_edges_[0].insert(string_edges_[0].end(), p.begin(), p.end());
                }
            }
        }
        auto p = origin->GetSimplePath(0);
        string_edges_[0].insert(string_edges_[0].end(), p.begin(), p.end());

        if (out_node_->IsType("loop") ) {
            auto n = static_cast<LoopNode*>(out_node_);
            auto back = n->GetBackward();
            for (auto e : back) {
                if (e->IsType("path")) {
                    static_cast<PathEdge*>(e)->IdentifySimplePaths(sg);
                    auto p = static_cast<PathEdge*>(e)->GetSimplePath(0);
                    string_edges_[0].insert(string_edges_[0].end(), p.begin(), p.end());
                }
            }
        }

        for (size_t i = 1; i < origin->SimplePathSize(); ++i) {
            string_edges_.push_back(origin->GetSimplePath(i));
        }
    }
}

}