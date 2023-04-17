#include "simplifier.hpp"


namespace fsa {

EdgeScore GetEdgeScore(const BaseEdge* e, ReadVariants* rvs) {
    std::array<int, 3> score = { (int)e->Score(), 0, 0 };     //  length, XIANGTONG BUTONG SNP
    if (rvs != nullptr) {
        auto mm = rvs->Test(*e->ol_);
        score[1] = mm[0]; 
        score[2] = mm[1];
    } 
    return score;
};


bool EdgeScoreSignificantlyGreater(const EdgeScore &a, const EdgeScore &b, int count, int rate) {

    // if (((a[1] - a[2]) - (b[1] - b[2]))  / 2>= threshold_count) {
    //     double ra = (a[1] + a[2]) > 0 ? (a[1] - a[2])*1.0 / (a[1] + a[2]) : 0.0;
    //     double rb = (b[1] + b[2]) > 0 ? (b[1] - b[2])*1.0 / (b[1] + b[2]) : 0.0;
    //     if ((ra - rb)/2 >= threshold_rate) {
    //         return true;
    //     }
    // }
    return (a[1] - a[2]) - (b[1] - b[2]) >= std::max<int>((a[1] + a[2]) * 0.5 + (b[1]+b[2]) * 0.5, 6);
    if (((a[1] - a[2]) - (b[1] - b[2]))  / 2>= count) {
        double ra = (a[1] + a[2]) > 0 ? (a[1] - a[2])*1.0 / (a[1] + a[2]) : 0.0;
        double rb = (b[1] + b[2]) > 0 ? (b[1] - b[2])*1.0 / (b[1] + b[2]) : 0.0;
        if ((ra - rb)/2 >= rate) {
            return true;
        }
    }

    return false;
}

Simplifier::~Simplifier() {
}


bool Simplifier::ParseParameter(const std::vector<std::string> &param) {
    assert(param.size() >= 1);
    return false;
}


void Simplifier::Debug(const char* const format, ...) const {
    if (DUMPER.IsWorking()) {
        va_list arglist;
        va_start(arglist, format);
        DUMPER[name_].Write(format, arglist);
        va_end(arglist);
    }
}

void Simplifier::DebugEdges(const std::string& msg, const std::unordered_set<BaseEdge*>& edges) {
    for (auto e : edges) {
        Debug("%s: %s\n", msg.c_str(), ToString(e).c_str());
    }
}

std::list<SgEdge*> Simplifier::GetInLinearPath(SgEdge* e) {
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


std::list<const SgEdge*> Simplifier::GetInLinearPath(const SgEdge* e) {
    std::list<const SgEdge*> edges {e};
    std::unordered_set<const SgEdge*> done {e};

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


std::list<SgEdge*> Simplifier::GetOutLinearPath(SgEdge* e) {
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


std::list<const SgEdge*> Simplifier::GetOutLinearPath(const SgEdge* e) {
    std::list<const SgEdge*> edges {e};
    std::unordered_set<const SgEdge*> done {e};

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


bool Simplifier::IsCross(SgNode* in0) const {
    Debug("in0: %s, %zd\n", ToString(in0).c_str(), in0->OutDegree());
   if (in0->OutDegree() == 2) {
        std::list<SgEdge*> in0_p0 = GetOutLinearPath(in0->OutEdge(0));
        std::list<SgEdge*> in0_p1 = GetOutLinearPath(in0->OutEdge(1));
        SgNode* out0 = in0_p0.back()->OutNode();
        SgNode* out1 = in0_p1.back()->OutNode();
        Debug("IsCross: OutNode = %s:%zd, %s:%zd\n", ToString(out0).c_str(),in0_p0.size(), ToString(out1).c_str(),in0_p1.size());
        if (out0 != out1 && out0->InDegree() == 2 && out1->InDegree() == 2) {
            SgNode* in1 = nullptr;
            for (size_t i = 0; i < out0->InDegree(); ++i) {
                std::list<SgEdge*> ip0 = GetInLinearPath(out0->InEdge(i));
                Debug("find in1: %s, %s \n", ToString(ip0.front()->InNode()).c_str(), ToString(out0->InEdge(i)).c_str());
                if (ip0.front()->InNode() != in0) {
                    in1 = ip0.front()->InNode();
                    break;
                }
            }
            if (in1 != nullptr && in1->OutDegree() == 2) {
                Debug("IsCross: InNode = %s, %s\n", ToString(in0).c_str(), ToString(in1).c_str());
                std::list<SgEdge*> in1_p0 = GetOutLinearPath(in1->OutEdge(0));
                std::list<SgEdge*> in1_p1 = GetOutLinearPath(in1->OutEdge(1));
                if ((in1_p0.back()->OutNode() == out0 &&  in1_p1.back()->OutNode() == out1) ||
                    (in1_p0.back()->OutNode() == out1 &&  in1_p1.back()->OutNode() == out0)) {
                    
                    bool tl_00 = TestLength(in0_p0);
                    bool tl_01 = TestLength(in0_p1);
                    bool tl_10 = TestLength(in1_p0);
                    bool tl_11 = TestLength(in1_p1);
                    bool te = TestExtends(in0, in1, out0, out1);
                    Debug("TestLengthExtend: TestLength(%d, %d, %d, %d), TestExtend(%d)\n", 
                        tl_00, tl_01, tl_10, tl_11, te);
                    return tl_00 && tl_01 && tl_10 && tl_11 && te;
                }
                
            }
        }
    }
    return false;
}



int Simplifier::PathCoreLength(const std::list<const SgEdge*> &path) const {
    
    int length = std::accumulate(path.begin(), path.end(), 0, [](int a, const SgEdge* b) {
        return a + (int)b->Length();
    });
    assert(path.back()->OutNode()->IsType(PathNode::TypeName()));
    Seq::Id id = path.back()->OutNode<PathNode>()->GetBaseNode()->ReadId();
    length -= ori_graph_.GetAsmData().GetReadStore().GetSeqLength(id);
    return length;
}

bool Simplifier::TestLength(const std::list<SgEdge*> &path) const {

    size_t max_length_ = 1000000;
    size_t max_node_size_ = 10;
    size_t node_size = SgGraph::PathNodeSize(path);;
    int length = 0;
    for (auto &p : path) {
        length += p->Length();
    }
    assert(path.back()->OutNode()->IsType(PathNode::TypeName()));
    Seq::Id id = path.back()->OutNode<PathNode>()->GetBaseNode()->ReadId();
    length -= ori_graph_.GetAsmData().GetReadStore().GetSeqLength(id);

    Debug("TestLength: Length(%d <= %zd) && NodeSize(%zd <= %zd)\n", length, max_length_, node_size, max_node_size_);
    return length <= (int)max_length_;
}

bool Simplifier::TestExtends(SgNode* in0, SgNode* in1, SgNode* out0, SgNode* out1) const {
    if (in0->InDegree() != 1 || in1->InDegree() != 1 || out0->OutDegree() != 1 || out1->OutDegree() != 1) return false;

    const int NODESIZE = 5;

    auto path_size = [](const std::list<const SgEdge*>& path) {
        size_t ps = 0;
        for (auto &p : path) {
            ps += p->NodeSize();
        }
        return ps;
    };

    auto in0_path = GetInLinearPath((const SgEdge*)in0->InEdge(0));
    auto in1_path = GetInLinearPath((const SgEdge*)in1->InEdge(0));
    auto out0_path = GetOutLinearPath((const SgEdge*)out0->OutEdge(0));
    auto out1_path = GetOutLinearPath((const SgEdge*)out1->OutEdge(0));
    
    Debug("nodesize1: %zd, %zd  %zd %zd\n", in0_path.size(), in1_path.size(), out0_path.size(), out1_path.size());
    Debug("nodesize2: %zd, %zd  %zd %zd\n", path_size(in0_path), path_size(in1_path), path_size(out0_path), path_size(out1_path));
    if (path_size(in0_path) <= NODESIZE && (in0_path.front()->InNode() == out0 || in0_path.front()->InNode() == out1)) return false;
    if (path_size(in1_path) <= NODESIZE && (in1_path.front()->InNode() == out0 || in1_path.front()->InNode() == out1)) return false;
    
    if (path_size(out0_path) <= NODESIZE && (out0_path.back()->OutNode() == in0 || out0_path.back()->OutNode() == in1)) return false;
    if (path_size(out1_path) <= NODESIZE && (out1_path.back()->OutNode() == in0 || out1_path.back()->OutNode() == in1)) return false;

    auto test_in_nodesize = [this, path_size](const SgNode* n) {
        assert(n->OutDegree() == 2 && n->InDegree() == 1);
        auto out_path0 = GetOutLinearPath(n->OutEdge(0));
        auto out_path1 = GetOutLinearPath(n->OutEdge(1));
        auto in_path = GetInLinearPath(n->InEdge(0));
        Debug("inedge lenght: %d > max(%d,%d)\n",  PathCoreLength(in_path), PathCoreLength(out_path0), PathCoreLength(out_path1));
        return PathCoreLength(in_path) > std::max(PathCoreLength(out_path0), PathCoreLength(out_path1)) &&
               path_size(in_path) > std::max(path_size(out_path0), path_size(out_path1));
        
    };

    auto test_out_nodesize = [this, path_size](const SgNode* n) {
        assert(n->InDegree() == 2 && n->OutDegree() == 1);
        auto in_path0 = GetInLinearPath(n->InEdge(0));
        auto in_path1 = GetInLinearPath(n->InEdge(1));
        auto out_path = GetOutLinearPath(n->OutEdge(0));
        Debug("outedge lenght: %d > max(%d,%d)\n",  PathCoreLength(out_path), PathCoreLength(in_path0), PathCoreLength(in_path1));
        return PathCoreLength(out_path) > std::max(PathCoreLength(in_path0), PathCoreLength(in_path1)) &&
               path_size(out_path) > std::max(path_size(in_path0), path_size(in_path1));
        
    };

    //if (!test_in_nodesize(in0)) return false;
    //if (!test_in_nodesize(in1)) return false;

    //if (!test_out_nodesize(out0)) return false;
    //if (!test_out_nodesize(out1)) return false;

    auto t_in_0 = test_in_nodesize(in0);
    auto t_in_1 = test_in_nodesize(in1);
    auto t_out_0 = test_out_nodesize(out0);
    auto t_out_1 = test_out_nodesize(out1);
    

    return t_in_0 && t_in_1 && t_out_0 && t_out_1;

}


bool Simplifier::TestCrossExtends(const std::list<const SgEdge*> in0_2_out0, const std::list<const SgEdge*> in0_2_out1, 
                        const std::list<const SgEdge*> in1_2_out0, const std::list<const SgEdge*> in1_2_out1) {

    const SgNode* in0 = in0_2_out0.front()->InNode();
    const SgNode* in1 = in1_2_out0.front()->InNode();
    const SgNode* out0 = in0_2_out0.back()->OutNode();
    const SgNode* out1 = in0_2_out1.back()->OutNode();
    // 它以及是一个cross
    assert(in0->InDegree() == 1 && in1->InDegree() == 1 && out0->OutDegree() == 1 && out1->OutDegree() == 1);
    assert( in0 == in0_2_out1.front()->InNode() &&  in1 == in1_2_out1.front()->InNode() &&
           out0 == in1_2_out0.back()->OutNode() && out1 == in1_2_out1.back()->OutNode());


    auto other_2_in0 = GetInLinearPath(in0->InEdge(0));
    auto other_2_in1 = GetInLinearPath(in1->InEdge(0));
    auto out0_2_other = GetInLinearPath(out0->OutEdge(0));
    auto out1_2_other = GetInLinearPath(out1->OutEdge(0));
    Debug("Other2Cross: NodeSize(->in0=%zd|%zd, ->in1=%zd|%zd, ->out0=%zd|%zd, ->out1=%zd|%zd\n", 
        other_2_in0.size(), SgGraph::PathNodeSize(other_2_in0),
        other_2_in1.size(), SgGraph::PathNodeSize(other_2_in1),
        out0_2_other.size(), SgGraph::PathNodeSize(out0_2_other),
        out1_2_other.size(), SgGraph::PathNodeSize(out1_2_other));
    //Debug("nodesize1: %zd, %zd  %zd %zd\n", in0_path.size(), in1_path.size(), out0_path.size(), out1_path.size());
    //Debug("nodesize2: %zd, %zd  %zd %zd\n", path_size(in0_path), path_size(in1_path), path_size(out0_path), path_size(out1_path));
   
    // 检查是不是有循环
    const int NODESIZE = 5;
    if ((SgGraph::PathNodeSize(other_2_in0) <= NODESIZE && (other_2_in0.front()->InNode() == out0 || other_2_in0.front()->InNode() == out1)) ||
        (SgGraph::PathNodeSize(other_2_in1) <= NODESIZE && (other_2_in1.front()->InNode() == out0 || other_2_in1.front()->InNode() == out1)) ||
        (SgGraph::PathNodeSize(out0_2_other) <= NODESIZE && (out0_2_other.back()->OutNode() == in0 || out0_2_other.back()->OutNode() == in1)) ||
        (SgGraph::PathNodeSize(out1_2_other) <= NODESIZE && (out1_2_other.back()->OutNode() == in0 || out1_2_other.back()->OutNode() == in1))) {
        Debug("Detect circle error\n");
        return false;
    }

    int len_in0_2_out0 = PathCoreLength(in0_2_out0);
    int len_in0_2_out1 = PathCoreLength(in0_2_out1);
    int len_in1_2_out0 = PathCoreLength(in1_2_out0);
    int len_in1_2_out1 = PathCoreLength(in1_2_out1);

    int len_other_2_in0 = PathCoreLength(other_2_in0);
    int len_other_2_in1 = PathCoreLength(other_2_in1);
    int len_out0_2_other = PathCoreLength(out0_2_other);
    int len_out1_2_other = PathCoreLength(out1_2_other);
    
    Debug("TestCrossExtends(length): ->in0=%d, ->in1=%d, out0->=%d, out1->=%d, in0->out0=%d, in0->out1=%d, in1->out0=%d, in1->out1=%d\n", 
        len_other_2_in0, len_other_2_in1, len_out0_2_other, len_out1_2_other,
        len_in0_2_out0, len_in0_2_out1, len_in1_2_out0, len_in1_2_out1);

    if (len_other_2_in0 < len_in0_2_out0 || len_other_2_in0 < len_in0_2_out1 || 
        len_other_2_in1 < len_in1_2_out0 || len_other_2_in1 < len_in1_2_out1 ||
        len_out0_2_other < len_in0_2_out0 || len_out0_2_other < len_in1_2_out0 ||
        len_out1_2_other < len_in0_2_out1 || len_out1_2_other < len_in1_2_out1) {
        return false;
    }

    return true;
}


std::pair<const Overlap*, double> Simplifier::ReplaceHighQualityOverlap(const Overlap* ol, double threshold) {

    std::pair<const Overlap*, double> nol = {nullptr, 0.0};
    auto dup0 = ori_graph_.GetAsmData().dup_groups_.find(ol->a_.id);
    if (dup0 != ori_graph_.GetAsmData().dup_groups_.end()) {
        auto dup00 = dup0->second.find(ol->b_.id);
        if (dup00 != dup0->second.end()) {
            
            Debug("dup: %s %s, %zd\n", QueryStringById(ol->a_.id).c_str(), 
                QueryStringById(ol->b_.id).c_str(), dup00->second.size());
            auto loc = ol->Location(0);
            assert(loc == Overlap::Loc::Right || loc == Overlap::Loc::Left);
            for (auto o : dup00->second) {
                Debug("dup loc: %s %s %d, %d\n", QueryStringById(o->a_.id).c_str(), 
                    QueryStringById(o->b_.id).c_str(), loc, o->Location(0));
                if (o->Location(0) == loc) {
                    auto qual = static_cast<StringGraph&>(ori_graph_).GetOverlapQuality(*o);
                    Debug("dup qual: %s (%d, %d) - %s (%d, %d) %f > %f\n", QueryStringById(o->a_.id).c_str(), o->a_.start, o->a_.end,
                        QueryStringById(o->b_.id).c_str(), o->b_.start, o->b_.end, qual, threshold);
                    if (qual >= threshold) {
                        nol.first = o;
                        nol.second = qual;
                        break;
                    }
                }
            }
        }
    }
    return nol;
}

/** 如果两个节点的有多条长度（子节点数）<=3的的路径，则只保留第一个（排序根据节点名称） */
void DuplicateSimplifier::Running() {

    std::unordered_map<BaseEdge::ID, std::vector<PathEdge*>, BaseEdge::ID::Hash> dup_edges;
    std::vector<SgEdge*> cands = graph_.CollectEdges([&dup_edges](SgEdge* e) {
        if (e->IsType("simple")) {
            SimplePathEdge *se = static_cast<SimplePathEdge*>(e);
            
            if (!se->IsReduced() && se->path_.size() < 3) {
                auto id = BaseEdge::ID(0, e->InNode()->Id(), e->OutNode()->Id());
                auto iter = dup_edges.find(id);
                if (iter != dup_edges.end()) {
                    iter->second.push_back(se);
                }
                else {
                    dup_edges[id] = std::vector<PathEdge*>{ se };
                }
                
            }
        }
        return false;
    });

    std::unordered_set<PathEdge*> done;     // 判断对偶路径是否已经处理完成
    for (auto &i : dup_edges) {
        std::vector<PathEdge*> &dups = i.second;
        if (dups.size() > 1 && done.find(dups[0]) == done.end()) {
            done.insert(dups[0]);
            done.insert(graph_.ReverseEdge(dups[0]));
            for (size_t i = 1; i < dups.size(); ++i) {
                auto e = dups[i];
                e->Reduce("simple_dup", true);

                auto re = graph_.ReverseEdge(dups[i]);
                re->Reduce("simple_dup", true);
                
                done.insert(re);
                done.insert(e);
            }
        }
    }
}



bool IsOutEdgeInconsistent(const BaseEdge* e0, const BaseEdge* e1, size_t N, const PhaseInfoFile* pif) {
    
    auto get_nodes = [N](const BaseEdge* e) {
        std::vector<const BaseNode*> nodes;
        const BaseEdge* curr = e;
        for (size_t i = 0; i < N && curr != nullptr; ++i) {
            nodes.push_back(curr->OutNode());
            curr = nodes.back()->OutDegree() > 0 ? (const BaseEdge*)nodes.back()->OutEdge(0) : nullptr;
        }
        return nodes;
    };

    if (pif != nullptr) {
        std::vector<const BaseNode*> nodes0 = get_nodes(e0);
        std::vector<const BaseNode*> nodes1 = get_nodes(e1);

        for (auto n0 : nodes0) {
            for (auto n1 : nodes1) {
                if (pif->Contain(n0->ReadId(), n1->ReadId())) {
                    return true;
                }
            }
        }
    }

    return false;
}


bool IsInEdgeInconsistent(const BaseEdge* e0, const BaseEdge* e1, size_t N, const PhaseInfoFile* pif) {
    
    auto get_nodes = [N](const BaseEdge* e) {
        std::vector<const BaseNode*> nodes;
        const BaseEdge* curr = e;
        for (size_t i = 0; i < N && curr != nullptr; ++i) {
            nodes.push_back(curr->InNode());
            curr = nodes.back()->InDegree() > 0 ? (const BaseEdge*)nodes.back()->InEdge(0) : nullptr;
        }
        return nodes;
    };

    if (pif != nullptr) {
        std::vector<const BaseNode*> nodes0 = get_nodes(e0);
        std::vector<const BaseNode*> nodes1 = get_nodes(e1);

        for (auto n0 : nodes0) {
            for (auto n1 : nodes1) {
                if (pif->Contain(n0->ReadId(), n1->ReadId())) {
                    return true;
                }
            }
        }
    }

    return false;
}

} // namespace fsa