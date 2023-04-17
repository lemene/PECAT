#include "phase.hpp"

#include "phase_path.hpp"

namespace fsa {


bool PhaseCrossSimplifier::ParseParameters(const std::vector<std::string> &params) {
    assert(params[0] == "phase");

    for (size_t i = 1; i < params.size(); ++i) {
        auto it = SplitStringByChar(params[i], '=');
        if (it[0] == "sc") {
            opts_.min_snp_count = std::stoi(it[1]);
        } else if (it[0] == "sr") {
            opts_.min_snp_rate = std::stof(it[1]);
        } else if (it[0] == "rs") {
            opts_.min_read_support = std::stoi(it[1]);
        } else if (it[0] == "vc") {
            opts_.min_valid_count = std::stoi(it[1]);
        } else if (it[0] == "vr") {
            opts_.min_valid_rate = std::stof(it[1]);

        } else {
            if (!ParseParameter(it)) {
                return false;
            }
        }
    }
    return true;
}

std::string PhaseCrossSimplifier::GetParameters() const {
    std::ostringstream oss;
    oss << "sc=" << opts_.min_snp_count << ":"
        << "sr=" << opts_.min_snp_rate << ":" 
        << "rs=" << opts_.min_read_support << ":"
        << "vc=" << opts_.min_valid_count << ":"
        << "vr=" << opts_.min_valid_rate;
    return oss.str();
}

std::vector<std::vector<BaseNode*>> PhaseCrossSimplifier::CollectCross() {
    std::vector<std::vector<BaseNode*>> cands;
    std::unordered_set<BaseNode::ID,BaseNode::ID::Hash> done;

    auto check_extend = [this](BaseNode* n, int minlen, int minnode) {
        if (n->OutDegree() != 2) return false;
        if (!n->OutEdge<BaseEdge>(0)->TestOutExtend(minlen, minnode)) return false;
        if (!n->OutEdge<BaseEdge>(1)->TestOutExtend(minlen, minnode)) return false;
        return true;
    };
    auto nodes = graph_.CollectNodes([this](BaseNode* n) {
        return n->InDegree() == 2;
    });


    for (auto n : nodes) {
        if (done.find(n->Id()) != done.end()) continue;

        std::vector<BaseNode*> cand = {n};

        BaseNode* curr = n;
        int length = 0;
        while (curr->OutDegree() == 1 && curr->GetOutEdge(0)->OutNode()->InDegree() == 1) {
            auto e = curr->OutEdge<BaseEdge>(0);
            cand.push_back(e->OutNode());
            length += (e->Length() + graph_.ReverseEdge(e)->Length()) / 2;
            curr = e->OutNode();
        }


        if (curr->OutDegree() == 2) {
            auto rn = graph_.ReverseNode(n);
            if (check_extend(curr, length*2, cand.size()*2) && check_extend(rn, length*2, cand.size()*2) &&
                IsOutEdgeInconsistent(curr->OutEdge<BaseEdge>(0), curr->OutEdge<BaseEdge>(1), 3, graph_.GetAsmData().GetInconsistentOverlaps()) &&
                IsOutEdgeInconsistent(rn->OutEdge<BaseEdge>(0), rn->OutEdge<BaseEdge>(1), 3, graph_.GetAsmData().GetInconsistentOverlaps()))

            cands.push_back(cand);
            done.insert(cand.front()->Id());
            done.insert(BaseNode::ReverseId(cand.front()->Id()));
            done.insert(cand.back()->Id());
            done.insert(BaseNode::ReverseId(cand.back()->Id()));
        }
    }

    LOG(INFO)("Find PhaseCross %zd", cands.size());
    return cands;
}

void PhaseCrossSimplifier::Running() {
    assert(rvs_ != nullptr);

    for (size_t _ = 0; _ < 3; _++) {

    std::vector<std::vector<BaseNode*>> cands = CollectCross();

    for (auto &cand : cands) {
        CrossPhaser phaser(*this, cand);
        if (phaser.Phase()) {
            ReplaceCross(cand, phaser.paths);
        }
    }
    } // for (size_t _ = 0; _ < 3; _++) {
}


void PhaseCrossSimplifier::ReplaceCross(std::vector<BaseNode*>& cand, std::vector<PhasePath>& paths) {
    for (auto c : cand) {
        for (auto e : c->GetInEdges()) {
            if (!e->IsReduce()) {

                e->Reduce(BaseEdge::RT_PHASED);
                graph_.ReverseEdge(e)->Reduce(BaseEdge::RT_PHASED);
            }
        }
        for (auto e : c->GetOutEdges()) {
            if (!e->IsReduce()) {
                e->Reduce(BaseEdge::RT_PHASED);
                graph_.ReverseEdge(e)->Reduce(BaseEdge::RT_PHASED);
            }
        }
    }

    //
    for (auto &path : paths) {
        Debug("phasepath path:\n");
        auto & pathnode = path.tips;
        for (size_t i = 1; i < pathnode.size(); ++i) {
            auto eid = BaseEdge::CreateID(pathnode[i-1], pathnode[i]);
            auto eiter = graph_.edges_.find(eid);
            if (eiter != graph_.edges_.end()) {
                if (eiter->second->IsReduce()) {
                    eiter->second->Reactivate();
                    graph_.ReverseEdge(eiter->second)->Reactivate();
                }
            } else {
                auto o = graph_.GetAsmData().QueryOverlap(Seq::EndIdToId(pathnode[i-1]), Seq::EndIdToId(pathnode[i]));
                assert (o != nullptr);
                graph_.AddOverlap(o);
            }
        }

    }
}


PhasePath::PhasePath(AsmDataset &ad, Seq::EndId start, ReadVariants *r, Seq::Id alt, const PhaseCrossSimplifier::Options& _opts)
 : asmdata(ad), rvs(r), opts(_opts) {

    tips.push_back(start);
    reads.insert(Seq::EndIdToId(start));

    Seq::Id rid = Seq::EndIdToId(tips.back());
    int rend = Seq::End(tips.back());
    std::unordered_set<const Overlap*> ols = asmdata.GetBackOverlaps(rid, rend);

    AddVariants(rvs->GetVariants(rid));
    for (auto ol : ols) {
        auto qread = ol->GetOtherRead(rid);
        auto altol = asmdata.QueryOverlap(qread.id, alt);
        Debug("add ol %s (%d)\n", asmdata.QueryNameById(qread.id).c_str(), altol==nullptr);
        if (altol == nullptr)  {
            reads.insert(qread.id);
            AddVariants(rvs->GetVariants(qread.id));
        } 
    }

    //
    auto inconsist = asmdata.GetInconsistentOverlaps();
    assert(inconsist != nullptr);
    auto alt_inconsist = inconsist->Get(alt);
    for (auto r : alt_inconsist) {
        Debug("add alt %s\n", asmdata.QueryNameById(r).c_str());
        AddVariants(rvs->GetVariants(r));
    }
    PrintVariants();
}

void PhasePath::Extend(const std::vector<BaseNode*> cand, const std::vector<BaseNode*>& ends) {
    for (size_t i = 0; i < cand.size()*6; ++i) {

        std::unordered_set<const Overlap*> extend;
        const Overlap* best = FindBestExtendOverlap(tips.back(), extend);
        PrintVariants();
        if (Check(best, extend)) {
            Seq::Id rid = Seq::EndIdToId(tips.back());
            for (auto ol : extend) {
                auto qread = ol->GetOtherRead(rid);
                reads.insert(qread.id);
                AddVariants(rvs->GetVariants(qread.id));
            }

            MoveToNextTip(tips.back(), best, extend);
            if (ReachEnds(ends)) {
                Debug("ReachEnds\n");
                break;
            }
        } else {
            Debug("no check\n");
            break;
        }
    }
}

void PhasePath::ExtendWithOtherPath(const std::vector<BaseNode*> cand, const std::vector<BaseNode*>& ends, PhasePath& other) {
    tips.erase(tips.begin()+1, tips.end());
    reads.clear();
    visited.clear();
    dst.clear();

    std::unordered_set<Seq::Id> exclude;
    for (auto t : other.tips) {
        exclude.insert(Seq::EndIdToId(t));
    }

    for (size_t i = 0; i < cand.size()*6; ++i) {
        Seq::Id rid = Seq::EndIdToId(tips.back());
        int rend = Seq::End(tips.back());

        std::unordered_set<const Overlap*> extols;

        std::unordered_set<const Overlap*> ols = asmdata.GetExtendOverlaps(rid, rend);

        const Overlap* best = nullptr;
        int best_rr = 0;
        int best_extlen = 0;
        for (const auto &o : ols) {
            auto qread = o->GetOtherRead(rid);
            if (exclude.find(qread.id) != exclude.end()) continue;

            auto rr = rvs->Test(other.reads, qread.id);

            if (best == nullptr) {
                best = o;
                best_rr = rr[0] - rr[1];
                best_extlen = o->Extension(rid, rend);
            } else {
                if (rr[0] - rr[1] < best_rr || (rr[0] - rr[1] == best_rr && o->Extension(rid, rend) > best_extlen)) {
                    best_rr = rr[0] - rr[1];
                    best = o;
                    best_extlen = o->Extension(rid, rend);
                }
            }
            extols.insert(o);

        }

        if (best != nullptr) {
            MoveToNextTip(tips.back(), best, extols);
            if (ReachEnds(ends)) {
                break;
            }
        } else {
            break;
        }
    }

}


const Overlap* PhasePath::FindBestExtendOverlap(Seq::EndId start, std::unordered_set<const Overlap*> &extols) {
    Seq::Id rid = Seq::EndIdToId(start);
    int rend = Seq::End(start);
    std::unordered_set<const Overlap*> ols = asmdata.GetExtendOverlaps(rid, rend);

    const Overlap* best = nullptr;
    int best_rr = 0;
    int best_extend = 0;
    auto comfirm = GetComfirmVariants();
    for (const auto &o : ols) {
        auto rr = TestVariants(comfirm, o->GetOtherRead(rid).id);
        Debug("check %d %d, %s %s\n",rr[0], rr[1], asmdata.QueryNameById(o->a_.id).c_str(), asmdata.QueryNameById(o->b_.id).c_str());

        if (rr[0] >= std::max<int>(opts.min_snp_count, std::ceil(opts.min_snp_rate *(rr[0]+rr[1])))) {
            if (best == nullptr) {
                best = o;
                best_rr = rr[0] - rr[1];
                best_extend = o->Extension(rid, rend);
            } else {
                if (rr[0] - rr[1] > best_rr || (rr[0] - rr[1] == best_rr && o->Extension(rid, rend) > best_extend)) {
                    best_rr = rr[0] - rr[1];
                    best = o;
                    best_extend = o->Extension(rid, rend);
                }
            }
            extols.insert(o);
        }
    }
    if (best != nullptr)    Debug("best %zd %s %s\n", extols.size(), asmdata.QueryNameById(best->a_.id).c_str(), asmdata.QueryNameById(best->b_.id).c_str());
    else Debug("best 0\n");
    return best;
}

bool PhasePath::Check(const Overlap* best,  std::unordered_set<const Overlap*> &extols) {
    return best != nullptr && (int)extols.size() >= opts.min_read_support;
}

void PhasePath::MoveToNextTip(Seq::EndId tip, const Overlap* best,  std::unordered_set<const Overlap*> &extols) {
    Seq::Id rid = Seq::EndIdToId(tip);
    int rend = Seq::End(tip);

    for (auto o : extols) {
        auto qread = o->GetOtherRead(rid);
        int end = 0;
        if (rend == 0) {
            end = o->SameDirect() ? 0 : 1;
        } else {
            end = o->SameDirect() ? 1 : 0;
        }
        visited.insert(Seq::IdToEndId(qread.id, end));
        reads.insert(qread.id);
    }

    auto qread = best->GetOtherRead(rid);
    rid = qread.id;
    if (rend == 0) {
        rend = best->SameDirect() ? 0 : 1;
    } else {
        rend = best->SameDirect() ? 1 : 0;
    }
    tips.push_back(Seq::IdToEndId(rid, rend));

}

bool PhasePath::ReachEnds(const std::vector<BaseNode*> &ends) {
    for (auto end : ends) {
        if (visited.find(end->Id().MainNode()) != visited.end()) {
            dst.push_back(end->Id().MainNode());
        }
    }

    if (dst.size() == 1) {
        tips.back() = dst[0];
    }

    return dst.size() > 0;
}

bool PhasePath::IsIndependentPath(const std::vector<PhasePath> &path) {
    assert(path.size() == 2);
    return path[0].dst.size() == 1 && path[1].dst.size() == 1 && path[0].dst[0] != path[1].dst[0];
}


void PhasePath::AddVariants(const std::vector<ReadVariants::Variants>* vars) {
    if (vars == nullptr) return;

    for (auto &v : *vars) {
        auto &ctg = variants[v.contig];
        for (auto &i : v.vars) {
            auto &s = ctg[i.first];
            bool done = false;
            for (auto &is : s) {
                if (is[0] == i.second) {
                    is[1] += 1;
                    done = true;
                    break;
                }
            }
            if (!done) {
                s.push_back({i.second, 1});
            }
        }
    }
}

void PhasePath::PrintVariants() const {
    for (auto &ctg : variants) {

        std::vector<int> keys(ctg.second.size());
        std::transform(ctg.second.begin(), ctg.second.end(), keys.begin(), [](const std::pair<int, std::vector<std::array<int, 2>>> &pair) { return pair.first;});

        std::sort(keys.begin(), keys.end());
        for (auto k : keys) {
            auto &v = ctg.second.find(k)->second;
            Debug("%d", k);
            for (const auto &iv : v) {
                Debug("-(%d,%d) ", iv[0], iv[1]);
            }
            Debug(" ");

        }
        Debug("\n");
    }
}

std::array<int,2> PhasePath::TestVariants(const std::unordered_map<int, std::unordered_map<int, int>>& avars, const std::vector<ReadVariants::Variants>* bvars) const {

    std::array<int,2> result  {0, 0};
    if  (bvars != nullptr) {
        for (auto & ibv : *bvars) {
            std::array<int,2> rs = {0, 0};
            auto actg = avars.find(ibv.contig);
            if (actg != avars.end()) {
                for (auto ib : ibv.vars) {
                    auto ia = actg->second.find(ib.first);
                    if (ia != actg->second.end()) {
                        if (ia->second == ib.second) {
                            rs[0] ++;
                        } else {
                            rs[1] ++;
                        }

                    }
                }
            }
            if (result[0] + result[1] < rs[0] + rs[1]) {
                result = rs;
            }
        }
    }
    return result;
}

std::array<int,2> PhasePath::TestVariants(const std::unordered_map<int, std::unordered_map<int, int>>& avar, int name) const {
    return TestVariants(avar, rvs->GetVariants(name));

}

std::unordered_map<int, std::unordered_map<int, int>> PhasePath::GetComfirmVariants() const {
    auto important = [this](std::vector<std::array<int, 2>> d) {


        std::sort(d.begin(), d.end(), [](const std::array<int,2> &a, const std::array<int,2> &b) {
            return (a[0] >= 0 && b[0] >= 0) ? a[1] > b[1] : a[0] > b[0];
        });

        int sum = std::accumulate(d.begin(), d.end(), 0, [](int a, const std::array<int,2>& b) {
            return a + (b[0] >= 0 ? b[1] : 0);        // remove seq error (-1)
        });

        assert(d.size() > 0);
        if (d[0][1] >= sum*opts.min_valid_rate && d[0][1] >= opts.min_valid_count) { // TODO PARAM NAME
            return d[0][0];
        } else {
            return -1;
        }

    };

    std::unordered_map<int, std::unordered_map<int, int>> vars;

    for (auto &i : variants) {
        auto &cv = vars[i.first];

        for (auto &j : i.second) {
            auto b = important(j.second);
            Debug("comfirm: %d %d\n", j.first, b);
            if (b >= 0) {
                cv[j.first] = b;
            }
        }

    }

    return vars;
}

void PhasePath::Debug(const char* const format, ...) const {
    if (DUMPER.IsWorking()) {
        va_list arglist;
        va_start(arglist, format);
        DUMPER["phase"].Write(format, arglist);
        va_end(arglist);
    }
}

CrossPhaser::CrossPhaser(PhaseCrossSimplifier& owner, const std::vector<BaseNode*> cand)
 : owner_(owner), cand_(cand) {
}

bool CrossPhaser::Phase() {

    auto& cand = cand_;

    if ((int)cand.size() > owner_.max_cand_size) return false;
    Debug("cand(%zd): %s->%s\n", cand.size(), owner_.graph_.GetAsmData().QueryNameById(cand.front()->ReadId()).c_str(),
        owner_.graph_.GetAsmData().QueryNameById(cand.back()->ReadId()).c_str());

    for (auto e : cand.back()->GetOutEdges()) {
        ends.push_back(e->OutNode());
    }

    Debug("start dophase\n");
    for (auto e : cand.front()->GetInEdges()) {
        Debug(" -- edge: %s %s\n", owner_.graph_.GetAsmData().QueryNameById(e->InNode()->ReadId()).c_str(),
            owner_.graph_.GetAsmData().QueryNameById(e->OutNode()->ReadId()).c_str());

        auto alt = std::find_if(cand.front()->GetInEdges().begin(), cand.front()->GetInEdges().end(), [e](BaseEdge* a) {
            return a != e;
        });
        if (alt == cand.front()->GetInEdges().end()) return false;
        assert(alt != cand.front()->GetInEdges().end());

        paths.push_back(PhasePath(owner_.graph_.GetAsmData(), e->InNode()->Id().MainNode(), owner_.rvs_, (*alt)->InNode()->ReadId(), owner_.opts_));
        PhasePath &path = paths.back();

        path.Extend(cand, ends);
    }

    Debug("end find %zd %zd\n", paths[0].dst.size(), paths[1].dst.size());
    if (paths.size() != 2) {
        LOG(WARNING)("phasepath != 2 indegree = %zd", cand.front()->InDegree());
        return false;
    }
    assert(paths.size() == 2);
    if (paths[0].dst.size() == 1 && paths[1].dst.size() != 1) {
        paths[1].ExtendWithOtherPath(cand, ends, paths[0]);
        Debug("phasepath amb 1\n");
    } else if (paths[1].dst.size() == 1 && paths[0].dst.size() != 1) {
        paths[0].ExtendWithOtherPath(cand, ends, paths[1]);
        Debug("phasepath amb 0\n");
    }

    return PhasePath::IsIndependentPath(paths);
}

void CrossPhaser::Debug(const char* const format, ...) const {
    if (DUMPER.IsWorking()) {
        va_list arglist;
        va_start(arglist, format);
        DUMPER["phase"].Write(format, arglist);
        va_end(arglist);
    }
}


} // namespace fsa