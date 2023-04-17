#include "contig_generator.hpp"

#include <vector>

#include "../utility.hpp"
#include "../phase/hic_read_infos.hpp"
namespace fsa {

void ContigGenerator::IdentifyContigs() {

    const auto& paths = path_graph_.GetPaths();

    for (size_t i = 0; i < paths.size(); i += 2) {
        contigs.push_back(Contig(this, i/2, paths[i], string_graph_));
        contigs.back().PhasedReads(dataset_.GetInconsistentOverlaps());
    }
    std::sort(contigs.begin(), contigs.end(), [](const Contig& a, const Contig &b) { return a.Weight() > b.Weight(); });
    LOG(INFO)("Contig size: %zd", contigs.size());
}

void ContigGenerator::FindInconsistentOverlap(const Contig &a, const Contig &b) {

}

bool ContigGenerator::IsTrivialContig(const Contig& ctg) {

    size_t threshold = opts_.max_trivial_length;
    auto &cluster = path_graph_.clusters_[path_graph_.edge2cluster_[ctg.path_->front()]];
    DUMPER["asm"]("Contig(%d,%d): Cluster(%zd): size %zd,%zd\n", ctg.id_, ctg.Weight(), path_graph_.edge2cluster_[ctg.path_->front()], cluster.Size(), cluster.Length());
    if (cluster.Length() < threshold) {
        DUMPER["asm"]("Cluster SetTrivial0 %zd", ctg.id_);
        return true;
    } else {
        if (ctg.Length() < threshold) {
            if (cluster.nontrivial.find(ctg.path_->front()) == cluster.nontrivial.end()) {
                DUMPER["asm"]("Cluster SetTrivial1 %zd", ctg.id_);
                return true;
            }
        }
    }
    return false;
}

void ContigGenerator::GroupContigs() {
    // contigs are sorted by length
    for (size_t i = 0; i < contigs.size(); ++i) {

        if (IsTrivialContig(contigs[i])) {
            contigs[i].SetTrivial();
        } else {
            bool assigned = false;
            for (auto &g : groups_) {
                if (!opts_.phased.empty()) {
                    for (auto c : g.ctgs) {
                        assert(!c->IsTrivial());
                        if (contigs[i].IsDiploid(*c, opts_)) {
                            g.ctgs.push_back(&contigs[i]);
                            assigned = true;
                            break;
                        }
                    }
                }
                if (assigned) break;
            }
            if (!assigned) {
                groups_.push_back(Group());
                groups_.back().ctgs.push_back(&contigs[i]);
            }

        }

        if (!opts_.phased.empty()) {
            int pri = -1;
            double score = 0;
            for (size_t j = i; j > 0; --j) {
                if (contigs[j-1].IsTrivial()) continue;
                if (contigs[i].IsDiploid(contigs[j-1], opts_)) {
                    auto s = contigs[i].DiploidScore(contigs[j-1]);
                    if (s > score) {
                        pri = j-1;
                        score = s;
                    }
                }
            }
                    
            if (pri >= 0) {
                contigs[i].SetPrimary(contigs[pri]);
                contigs[pri].AddAlternate(contigs[i]);
            }
            
        }
    }
    
    DUMPER["asm"]("ctggroup size: %zd\n", groups_.size());
    for (size_t i = 0; i < groups_.size(); ++i) {
        auto &g = groups_[i];
        for (size_t j = 0; j < g.ctgs.size(); ++j) {
            DUMPER["asm"]("ctggroup: g(%zd): id=%d, bubbles=%zd\n", i, g.ctgs[j]->id_, g.ctgs[j]->acontigs.size());
        }
    }
}

void ContigGenerator::Save() {
    IdentifyContigs();
    GroupContigs();


    if (opts_.contig_format.find("prialt") != std::string::npos || opts_.contig_format.empty()) {
        SaveContigs(contigs, "prialt");
    }

    if (opts_.contig_format.find("dual") != std::string::npos) {
        PhaseBubbles(contigs);
        SaveContigs(contigs, "dual");
    }

}

void ContigGenerator::PhaseBubbles(std::vector<Contig> &contigs) {
    
    if (!opts_.hic_info.empty()) {
        HicReadInfos hic_infos(dataset_.string_pool_);
        hic_infos.Load(opts_.hic_info);
        if (dataset_.GetReadVariants() != nullptr) {
            for (auto &c : contigs) {
                c.PhaseBubbles(hic_infos, *dataset_.GetReadVariants());
            }
        }
    } else {
        if (dataset_.GetReadVariants() != nullptr) {
            for (auto &c : contigs) {
                c.PhaseBubbles();
            }
        }
    }
}

void ContigGenerator::SaveContigs(std::vector<Contig>& contigs, const std::string& format) {

    std::mutex mutex;
    std::string seqs0_name = format == "prialt" ? "primary.fasta" : "haplotype_1.fasta";
    std::string tiles0_name = format == "prialt" ? "primary_tiles" : "haplotype_1_tiles";
    std::string seqs1_name = format == "prialt" ? "alternate.fasta" : "haplotype_2.fasta";
    std::string tiles1_name = format == "prialt" ? "alternate_tiles" : "haplotype_2_tiles";
    GzFileWriter fcontig_seqs(OutputPath(seqs0_name));
    GzFileWriter fcontig_tiles(OutputPath(tiles0_name));
    GzFileWriter fcontig_seqs1(OutputPath(seqs1_name));
    GzFileWriter fcontig_tiles1(OutputPath(tiles1_name));
    GzFileWriter fbubble_seqs(OutputPath("rest.fasta"));
    GzFileWriter fbubble_tiles(OutputPath("rest_tiles"));
    

    auto combine_func = [&](std::ostringstream& oss_pri, std::ostringstream& oss_pri_tiles, 
                            std::ostringstream& oss_alt, std::ostringstream& oss_alt_tiles, 
                            std::ostringstream& oss_oth, std::ostringstream& oss_oth_tiles) {
        std::lock_guard<std::mutex> lock(mutex);
        fcontig_seqs << oss_pri.str();          oss_pri.str("");
        fcontig_tiles << oss_pri_tiles.str();   oss_pri_tiles.str("");
        fcontig_seqs1 << oss_alt.str();         oss_alt.str("");
        fcontig_tiles1 << oss_alt_tiles.str();  oss_alt_tiles.str("");
        fbubble_seqs << oss_oth.str();          oss_oth.str("");
        fbubble_tiles << oss_oth_tiles.str();   oss_oth_tiles.str("");
    };

    std::atomic<size_t> index { 0 };
    auto work_func = [&contigs, &index, this, combine_func, &format](size_t id) {
        std::ostringstream oss_pri;
        std::ostringstream oss_pri_tiles;
        std::ostringstream oss_alt;
        std::ostringstream oss_alt_tiles;
        std::ostringstream oss_oth;
        std::ostringstream oss_oth_tiles;

        for (size_t curr = index.fetch_add(1); curr < contigs.size(); curr = index.fetch_add(1)) {
            if (contigs[curr].IsTrivial())  {
                if(this->ConstructContig(contigs[curr].pcontig).size() < (size_t)opts_.max_trivial_length) {
                    contigs[curr].SavePrimary(oss_oth, oss_oth_tiles);
                    contigs[curr].SaveBubbles(oss_oth, oss_oth_tiles, oss_oth, oss_oth_tiles, format);
                    continue;
                }
            }
            if (contigs[curr].IsPrimary()) {
                contigs[curr].SavePrimary(oss_pri, oss_pri_tiles);
                contigs[curr].SaveBubbles(oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles, format);
                
                // if (format == "dual" && contigs[curr].acontigs.size() == 0 && contigs[curr].alts_.size() == 0) {
                //     contigs[curr].SavePrimary(oss_alt, oss_alt_tiles);
                //     contigs[curr].SaveBubbles(oss_pri, oss_pri_tiles, oss_oth, oss_oth_tiles, format);
                // }
            } else {
                // if (contigs[curr].GetPrimary()->IsPrimary()) {
                //     contigs[curr].SavePrimary(oss_alt, oss_alt_tiles);
                //     contigs[curr].SaveBubbles(oss_pri, oss_pri_tiles, oss_oth, oss_oth_tiles, format);
                // } else {
                //     if ( // !contigs[curr].GetPrimary()->GetPrimary()->IsPrimary() ||
                //         contigs[curr].IsDiploid(*contigs[curr].GetPrimary()->GetPrimary(), opts_) ||
                //         contigs[curr].IsOverlapped(*contigs[curr].GetPrimary()->GetPrimary(), opts_)) {

                //         contigs[curr].SavePrimary(oss_oth, oss_oth_tiles);
                //         contigs[curr].SaveBubbles(oss_oth, oss_oth_tiles, oss_oth, oss_oth_tiles, format);
                //     } else {
                //         contigs[curr].SavePrimary(oss_pri, oss_pri_tiles);
                //         contigs[curr].SaveBubbles(oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles, format);
                //     }
                // }
                size_t count = 1;
                auto pp = contigs[curr].GetPrimary();
                while (!pp->IsPrimary()) {
                    count ++;
                    pp = pp->GetPrimary();
                }

                if (count % 2) {
                    contigs[curr].SavePrimary(oss_alt, oss_alt_tiles);
                    contigs[curr].SaveBubbles(oss_pri, oss_pri_tiles, oss_oth, oss_oth_tiles, format);
                } else {
                    contigs[curr].SavePrimary(oss_pri, oss_pri_tiles);
                    contigs[curr].SaveBubbles(oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles, format);
                }

            }
        }
        combine_func(oss_pri, oss_pri_tiles, oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles);
    };

    MultiThreadRun((size_t)opts_.thread_size, work_func);
    
}



std::string ContigGenerator::ConstructContigS(const std::list<BaseEdge*> &contig) {
    std::string seq;
    seq.reserve(1.1*std::accumulate(contig.begin(), contig.end(), 0, [](int a, const BaseEdge*b) {
        return a + b->Length();
    }));

    for (auto e : contig) {
        if (e == contig.front()) {
            if (e->InNode()->OutDegree() == 1) {
                int read = e->InNode()->ReadId();
                std::string readseq = *dataset_.GetReadStore().GetSeq(read).ToString();
                seq += e->InNode()->Id().End() == 1 ? readseq : Seq::ReverseComplement(readseq);
                seq += EdgeToSeq(e);
            } else {
                
                int read = e->OutNode()->ReadId();
                std::string readseq = *dataset_.GetReadStore().GetSeq(read).ToString();
                seq += e->OutNode()->Id().End() == 1 ? readseq : Seq::ReverseComplement(readseq);
            }
        } else if (e == contig.back()) {
            if (e->OutNode()->InDegree() == 1) {
                seq += EdgeToSeq(e);
            }
        } else {
            seq += EdgeToSeq(e);
        } 
    }

    return seq;
}

std::string ContigGenerator::ConstructContigStraight(const std::list<BaseEdge*> &contig) {
    std::string seq;
    
    auto first = contig.front()->InNode();
    int read = contig.front()->InNode()->ReadId();
    std::string readseq = *dataset_.GetReadStore().GetSeq(read).ToString();
    seq += first->Id().End() == 1 ? readseq : Seq::ReverseComplement(readseq);
    
    for (auto e : contig) {
        seq += EdgeToSeq(e);
    }
    return seq;
}


std::string ContigGenerator::ConstructContig(const std::list<BaseEdge*> &contig) {
    std::string seq;

    // if (contig.size() == 1) return seq;
    auto first = contig.front()->InNode();
    // first->OutDegree() == 0: never happen
    // first->OutDegree() == 1: read of first node should be add to this contig
    // first->OutDegree() >  1: read of first node are not determined how to add
    if (first->OutDegree() == 1) {  
        // first->InDegree() == 0: Add the whole read
        // first->InDegree() == 1: Add the whole read. It is a circle, so it is add when dealing last
        // first->InDegree() >  1: Add the shortest in_edge.
        if (first->InDegree() == 0) {
            int read = first->ReadId();
            std::string readseq = *dataset_.GetReadStore().GetSeq(read).ToString();
            seq += first->Id().End() == 1 ? readseq : Seq::ReverseComplement(readseq);

        } else if (first->InDegree() > 1) {
            auto pid = PathNode::CreateId(first->Id());
            auto pnode = path_graph_.QueryNode(pid);
            if (pnode != nullptr && pnode->InDegree() == 0) {
                int read = first->ReadId();
                std::string readseq = *dataset_.GetReadStore().GetSeq(read).ToString();
                seq += first->Id().End() == 1 ? readseq : Seq::ReverseComplement(readseq);
                DUMPER["asm"]("ConstructContig:PathNode:InDegree=0: %s\n", pid.ToString(dataset_.GetStringPool()).c_str());
            }
            // auto m = first->InEdge(0);
            // for (size_t i = 1; i < first->InDegree(); ++i) {
            //     if (first->InEdge(i)->Length() < m->Length()) {
            //         m = first->InEdge(i);
            //     }
            // }

            // seq += EdgeToSeq(static_cast<BaseEdge*>(m));
        }
    }

    for (auto e : contig) {
        if ( e != contig.back()) {
            seq += EdgeToSeq(e);
        }
    }

    auto last = contig.back()->OutNode();
    // last->InDegree() == 0: never happen
    // last->InDegree() == 1: Add the contig.back()
    // last->InDegree() >  1: contig.back() is repeat area, and should be converted to independent ctg. 
    
    auto m = last->InEdge(0);
    for (size_t i = 1; i < last->InDegree(); ++i) {
        if (last->InEdge(i)->Length() < m->Length()) {
            m = last->InEdge(i);
        }
    }

    if (m == contig.back()) {
        // last->OutDegree() == 0: It is a end point
        // last->OutDegree() == 1: It is a circle, 
        // last->OutDegree() >  1: It has some branches
        seq += EdgeToSeq(contig.back());
    }

    return seq;
}


std::vector<Seq::Tile> ContigGenerator::EdgesToTiles(const std::vector<BaseEdge*> &contig) {
    std::vector<Seq::Tile> tile;

    auto first = contig.front()->InNode();
    // first->OutDegree() == 0: never happen
    // first->OutDegree() == 1: read of first node should be add to this contig
    // first->OutDegree() >  1: read of first node are not determined how to add
    if (first->OutDegree() == 1) {  
        // first->InDegree() == 0: Add the whole read
        // first->InDegree() == 1: Add the whole read. It is a circle, so it is add when dealing last
        // first->InDegree() >  1: Add the shortest in_edge.
        if (first->InDegree() == 0) {
            int read = contig.front()->InNode()->ReadId();

            tile.push_back({read, first->Id().End() == 1 ? 0 : 1, 0, contig.front()->InRead()->len});

        } else if (first->InDegree() > 1) {
            auto m = std::min_element(first->GetInEdges().begin(), first->GetInEdges().end(), [](const BaseEdge* a, const BaseEdge*b) {
                return a->Length() < b->Length();
            });
            tile.push_back((*m)->GetTile());
        }
    }

    for (auto e : contig) {
        if ( e != contig.back()) {
            tile.push_back(e->GetTile());
        }
    }

    auto last = contig.back()->OutNode();
    // last->InDegree() == 0: never happen
    // last->InDegree() == 1: Add the contig.back()
    // last->InDegree() >  1: contig.back() is repeat area, and should be converted to independent ctg. 
    if (last->InDegree() == 1) {
        // last->OutDegree() == 0: It is a end point
        // last->OutDegree() == 1: It is a circle, 
        // last->OutDegree() >  1: It has some branches
        tile.push_back(contig.back()->GetTile());
    }

    return tile;
}

std::string ContigGenerator::EdgeToSeq(const BaseEdge *e) {

    int read = e->OutNode()->ReadId();
    std::string readseq = *dataset_.GetReadStore().GetSeq(read).ToString();

    Seq::Tile tile = e->GetTile();

    std::string a0, a1;
    if (tile.strand == 0) {
        a0 = readseq.substr(tile.start, tile.end-tile.start);
    } else {
        a0 = Seq::ReverseComplement(readseq.substr(tile.start, tile.end-tile.start));
    }

    return a0;
}

ContigGenerator::Contig::Contig(ContigGenerator* owner, size_t id, const std::list<PathEdge*> &path, StringGraph &sg) 
 : id_(id), path_(&path), owner_(owner) {
    
    auto path_len = [](const std::list<BaseEdge*>& path) {
        size_t len = 0;
        for (auto p : path) {
            len += p->Length();
        }
        return len;
    };

    for (auto p : path) {
        p->IdentifySimplePaths(sg); // TODO It may not be necessary to call
        assert(p->SimplePathSize() >= 1);

        std::vector<std::list<BaseEdge*>> ctgs;
        if (p->SimplePathSize() >= 1) {
            ctgs.push_back(p->GetSimplePath(0));

            if (p->IsType("semi")) {
                auto semi = static_cast<SemiBubbleEdge*>(p);
                std::vector<std::array<bool,2>> pairs;
                
                DUMPER["asm"]("ctgpathsiz: %zd, %zd\n", p->SimplePathSize(), semi->paths_.size());
                if (p->SimplePathSize() != semi->paths_.size()) {
                    for (auto &ip : semi->paths_) {
                        DUMPER["asm"]("ctgpsize ---");
                        for (auto &ipp : ip) {
                            DUMPER["asm"]("ctgpsize: %s\n", ipp->Id().ToString(owner->dataset_.GetStringPool()).c_str());
                        }
                    }
                    DUMPER["asm"]("ctgpathsiz: %zd, %zd\n", p->SimplePathSize(), semi->paths_.size());
                }
                assert(p->SimplePathSize() == semi->paths_.size());
                for (size_t i = 1; i < p->SimplePathSize(); ++i) {
                    auto &c = p->GetSimplePath(i);

                    int loc = -1;
                    if (semi->paths_[i].front()->InNode() == semi->paths_[0].front()->InNode()) {
                        loc = 0;
                    } else if (semi->paths_[i].back()->OutNode() == semi->paths_[0].back()->OutNode()) {
                        loc = 1;
                    }
                    assert(loc == 0 || loc == 1);

                    bool added = false;
                    for (size_t j = 0; j < pairs.size(); ++j) {
                        if (!pairs[j][loc]) {
                            if (loc == 0) {
                                ctgs[j+1].insert(ctgs[j+1].begin(), c.begin(), c.end());
                            } else {
                                assert(loc == 1);
                                ctgs[j+1].insert(ctgs[j+1].end(), c.begin(), c.end());
                            } 
                            pairs[j][loc] = true;
                            added = true;
                        }
                    }
                    if (!added) {
                        std::array<bool, 2> pair = { 0, 0};
                        pair[loc] = true;
                        pairs.push_back(pair);
                        ctgs.push_back(c);
                    }
                }
            } else {
                for (size_t i = 1; i < p->SimplePathSize(); ++i) {
                    ctgs.push_back(p->GetSimplePath(i));
                }
            }

        }
        if (ctgs.size() > 1) {
            if ( p->IsType("bubble") || 
                (p == path.front() && p->IsType("semi") && static_cast<SemiBubbleEdge*>(p)->Branchs() == 2) ||
                (p == path.back()  && p->IsType("semi") && static_cast<SemiBubbleEdge*>(p)->Branchs() == 1)) {

                DUMPER["asm"]("ctgswap -- : %zd, %s\n", id_, p->Type().c_str());
                auto mx = std::max_element(ctgs.begin(), ctgs.end(), [path_len](const std::list<BaseEdge*>& a, const std::list<BaseEdge*> &b) {
                    return path_len(a) < path_len(b);
                });
                if (ctgs.begin() != mx) {
                    std::swap(ctgs[0], *mx);
                }
            }
        }

        if (pcontig.size() > 0 && ctgs[0].size() > 0) {
            if (pcontig.back()->OutNode() != ctgs[0].front()->InNode()) {
                auto s = sg.GetEdge(pcontig.back()->OutNode()->Id(), ctgs[0].front()->InNode()->Id());
                //assert(s != nullptr);
                if (s != nullptr) pcontig.push_back(s);
            }
        }
        pcontig.insert(pcontig.end(), ctgs[0].begin(), ctgs[0].end()); 
        segs_.push_back(ctgs);       
        if (ctgs.size() > 1) {
            acontigs.push_back(std::make_pair(p, std::move(ctgs)));
        }
    }
    weight_ = path_len(pcontig);
    phasing_.assign(acontigs.size(), 0);
}

// c is primary
bool ContigGenerator::Contig::IsDiploid(const Contig& c, const AsmOptions& opts) const {
    size_t count = 0;
    for (auto r : reads) {
        if (c.vreads.find(r) != c.vreads.end()) {
            count ++;
        } else {
        }
    }
    DUMPER["asm"]("alt: %d -> %d: %zd/%zd\n", id_, c.id_, count, reads.size());
    return (count*1.0 / reads.size() >= opts.diploid_rate && (int)count >= opts.diploid_count) ||  count*1.0 / reads.size() >= 1.0;
}

double ContigGenerator::Contig::DiploidScore(const Contig& c) const {
    size_t count = 0;
    for (auto r : reads) {
        if (c.vreads.find(r) != c.vreads.end()) {
            count ++;
        } else {
        }
    }
    return count*1.0 / reads.size() ;
}



// c is primary
bool ContigGenerator::Contig::IsOverlapped(const Contig& c, const AsmOptions& opts) const {
    size_t count = 0;
    for (auto r : vreads) {
        if (c.vreads.find(r) != c.vreads.end()) {
            count ++;
        } 
    }
    
    DUMPER["asm"]("overlaped: %zd -> %zd: %zd/%zd\n", id_, c.id_, count, vreads.size());
    return count > 0 && count*1.0 / vreads.size() >= 0.5; // TODO Parameterization 

}

void ContigGenerator::Contig::PhasedReads(PhaseInfoFile *phased) {
    if (phased != nullptr) {
        for (auto p : pcontig) {
            // if (p == pcontig.front()) {
            //     if (p->InNode()->InDegree() == 0) {
            //         reads.insert(p->InRead()->id);
            //         auto rs = phased->Get(p->OutRead()->id);
            //         vreads.insert(rs.begin(), rs.end());
            //     }
            // }
            reads.insert(p->OutRead()->id);
            auto rs = phased->Get(p->OutRead()->id);
            vreads.insert(rs.begin(), rs.end());
        }

        for (auto & lsa: acontigs) {
            for (auto &as : lsa.second) {   
                for (auto a : as) {
                    //reads.insert(a->OutRead()->id);
                    auto rs = phased->Get(a->OutRead()->id);
                    vreads.insert(rs.begin(), rs.end());
                }
            }
        }
    }
}

void ContigGenerator::Contig::SavePrimary(std::ostream& fctg, std::ostream& ftile, uint8_t hap) {
    auto pcontig = ConstructPrimaryPath(hap);
    if (pcontig.size() == 0) return;

    std::vector<std::list<BaseEdge*>> pctgs;
    for (auto p : pcontig) {
        if (pctgs.size() == 0 || pctgs.back().back()->OutNode() != p->InNode()) {
            pctgs.push_back({p});
        } else {
            pctgs.back().push_back(p);
        }
    }

    for (size_t i = 0; i < pctgs.size(); i++) {
        auto& p = pctgs[i];
        std::vector<std::string> seqs = {owner_->ConstructContig(p)};
        assert(seqs.size() >= 1);
        std::ostringstream name;
        name << MainName() << (hap==0 ? "" : "_alt");
        if (i != 0) name << "_" << i;
        if (seqs[0].length() >= (size_t)owner_->opts_.min_contig_length) {
            fctg << ">" << name.str() << " " << Description(seqs[0].size()) << "\n" << seqs[0] << "\n";

            
            const StringGraph& sg = owner_->string_graph_;
            for (auto e : p) {
                ftile << name.str() << " edge=" << sg.IdToString(e->InNode()->Id()) << 
                        "~" << sg.IdToString(e->OutNode()->Id()) << "\n";
            }
        }

    }
}



void ContigGenerator::Contig::SaveBubbles(std::ostream &fctg,  std::ostream &ftile) {
    int ibubble = 1;

    for (const auto &bubble : acontigs) {
        const auto &paths = bubble.second;
        assert(paths.size() > 0);
        SaveBubbles(fctg, ftile, ibubble, paths);

        ibubble++;
    }
}


void ContigGenerator::Contig::SaveBubbles(std::ostream &fctg,  std::ostream &ftile, int id, const std::vector<std::list<BaseEdge*>>& paths) {
    assert(paths.size() > 0);

    // for (size_t i = 1; i < paths.size(); ++i) {
    //     std::string seq = ass.ConstructContigS(paths[i]);
    //     fctg << ">" << SubName(id, i) << " length=" << seq.size() << "\n" 
    //          << seq << "\n";

    //     for (auto p : paths[i]) {
    //         ftile << SubName(id, i) <<
    //                 " edge=" <<  ass.string_graph_.IdToString(p->InNode()->Id()) << 
    //                 "~" << ass.string_graph_.IdToString(p->OutNode()->Id()) << "\n";

    //     }
    // }
    assert(paths.size() >= 2);
    const auto &path = phasing_[id-1] == 0 ? paths[1] : paths[0];
    std::string seq = owner_->ConstructContigS(path);
    fctg << ">" << SubName(id, 1) << " length=" << seq.size() << "\n" 
            << seq << "\n";

    for (auto p : path) {
        ftile << SubName(id, 1) <<
                " edge=" <<  owner_->string_graph_.IdToString(p->InNode()->Id()) << 
                "~" << owner_->string_graph_.IdToString(p->OutNode()->Id()) << "\n";

    }
}

void ContigGenerator::Contig::SaveBubbles(std::ostream& fctg0, std::ostream& ftile0, std::ostream& fctg1, std::ostream& ftile1, const std::string& format) {
    if (format == "dual") {
        std::array<size_t, 2> count = {0, 0};
        for (const auto &bubble : acontigs) {
            const auto &paths = bubble.second;
            assert(paths.size() > 0);
            
            bool is_covered = IsCovered(paths);
            if (is_covered) {
                count[1] += paths.size();
            }
            count[0] += paths.size();

        }
        DUMPER["asm"]("ctg(%zd): overlap: %zd %zd", id_, count[0], count[1]);
        if (count[1] < 0.5 * count[0]) {
            SavePrimary(fctg0, ftile0, 1);
        }
    } else {
        int ibubble = 1;

        for (const auto &bubble : acontigs) {
            DUMPER["asm"]("CONTIG BUBBLE %d %d", id_, ibubble);
            const auto &paths = bubble.second;
            assert(paths.size() > 0);
            
            bool is_covered = IsCovered(paths);
            std::ostream &fctg =  !is_covered ? fctg0  : fctg1;
            std::ostream &ftile = !is_covered ? ftile0 : ftile1;

            SaveBubbles(fctg, ftile, ibubble, paths);

            ibubble++;
        }
    }
}

void ContigGenerator::Contig::PhaseBubbles(const HicReadInfos &infos, const ReadVariants& rvs) {
    if (acontigs.size() < 2) return;

    DUMPER["asm"]("phasing: ctg=%zd\n", id_);

    std::vector<size_t> scores(acontigs.size()*2*acontigs.size()*2, 0);

    std::unordered_map<SnpAllele, size_t, SnpAllele::Hash> alleles;
    std::unordered_set<SnpAllele, SnpAllele::Hash> dups;


    size_t index = 0;
    for (auto& actg : acontigs) {
        assert(actg.second.size() >= 2);

        std::array<std::unordered_map<SnpSite, std::array<uint8_t, 4>, SnpSite::Hash>, 2> ctg_vars;
        for (size_t i = 0; i < 2; ++i) {
            auto &path = actg.second[i];
            for (auto &p : path) {
                if (&p != &path.back()) {
                    const auto r = p->OutRead();
                    auto vars = rvs.GetVariants(r->id);
                    if (vars != nullptr) {
                        for (auto &vg : *vars) {
                            for (auto v : vg.vars) {                
                                if (v.second >= 0 && v.second < 4) {
                                    assert(vg.contig >= 0);
                                    SnpSite site = {(uint32_t)vg.contig, (uint32_t)v.first};
                                    ctg_vars[i][site][v.second]++;
                                }
                            }
                        }
                    }
                }
            }
        }

        // 保留相同的
        std::unordered_set<SnpSite, SnpSite::Hash> comm;
        for (auto &v : ctg_vars[0]) {
            if (ctg_vars[1].find(v.first) != ctg_vars[1].end()) {
                comm.insert(v.first);
            }
        }

        for (auto &c : comm) {
            std::vector<SnpAllele> pair;
            for (size_t i = 0; i < 2; ++i) {
                auto& v = ctg_vars[i][c];
                auto sum = std::accumulate(v.begin(), v.end(), 0);
                DUMPER["asm"]("ctg(%zd,%zd): %zd - %zd - %zd,%zd,%zd,%zd\n", index, i, c.ctg, c.offset, v[0], v[1], v[2], v[3]);
                for (size_t ib = 0; ib < v.size(); ++ib) {
                    if (v[ib] > sum*0.66) {
                        SnpAllele al = {c, (uint8_t)ib};
                        pair.push_back(al);
                    }
                }
            }
            if (pair.size() == 2 && pair[0] != pair[1]) {
                if (alleles.find(pair[0]) == alleles.end() && alleles.find(pair[1]) == alleles.end()) {
                    alleles[pair[0]] = index + 0;
                    alleles[pair[1]] = index + 1;
                } 
                // TODO DUP
            }
        }
        
        index += 2;
    }
    for (const auto& d : dups) {
        alleles.erase(d);
    }
    DUMPER["asm"]("alleles size = %zd\n", alleles.size());
    // for (auto &a: alleles) {
    //     printf("al: (%u,%u,%u)->%zd\n", a.first.site.ctg, a.first.site.offset, a.first.base, a.second);
    // }
    // fflush(stdout);
    DUMPER["asm"]("actg.size = %zd, index=%zd\n", acontigs.size(), index);
    for (auto& info : infos.GetInfos()) {
        std::array<std::vector<size_t>, 2> links; 
        for (size_t i = 0; i < 2; ++i) {
            for (auto& h : info.second.hic[i]) {
                for (auto al : h.alleles) {
                    auto iter = alleles.find(al);
                    if (iter != alleles.end()) {
                        links[i].push_back(iter->second);
                    }
                }
            }
        }

        for (auto s : links[0]) {
            for (auto t : links[1]) {
                //if (s != t) {
                    assert(s*acontigs.size()*2 + t < scores.size());
                    scores[s*acontigs.size()*2 + t] ++;
                    scores[t*acontigs.size()*2 + s] ++;
                //}
            }
        }
    }

    DUMPER["asm"]("print scores\n");
    for (size_t i = 0; i < acontigs.size()*2; ++i) {
        for (size_t j = 0; j < acontigs.size()*2; ++j) {
            DUMPER["asm"]("%zd, ", scores[i*acontigs.size()*2+j]);
        }
        DUMPER["asm"]("\n");
    }

    for (size_t _ = 0; _ < 1000; ++_) {
        DUMPER["asm"]("phasing(%zd), size=%zd\n", _, phasing_.size());
        std::vector<double> link_scores(phasing_.size(), 0.0);
        for (size_t  i = 0; i < link_scores.size(); ++i) {
            std::array<size_t, 2> count = {0, 0};
            for (size_t j = 0; j < link_scores.size(); ++j) {
                if (i == j) continue;
                if (phasing_[i] == phasing_[j]) {
                    count[0] += scores[i*2*acontigs.size()*2 + j*2];
                    count[0] += scores[(i*2+1)*acontigs.size()*2 + j*2+1];

                    count[1] += scores[i*2*acontigs.size()*2 + j*2+1];
                    count[1] += scores[(i*2+1)*acontigs.size()*2 + j*2];
                } else {
                    count[0] += scores[i*2*acontigs.size()*2 + j*2+1];
                    count[0] += scores[(i*2+1)*acontigs.size()*2 + j*2];
                    
                    count[1] += scores[i*2*acontigs.size()*2 + j*2];
                    count[1] += scores[(i*2+1)*acontigs.size()*2 + j*2+1];
                }
                DUMPER["asm"]("phasing(%zd), count index(%zd, %zd, %zd, %zd)\n", _, 
                    i*2*acontigs.size()*2 + j*2, (i*2+1)*acontigs.size()*2 + j*2+1, i*2*acontigs.size()*2 + j*2+1, (i*2+1)*acontigs.size()*2 + j*2);
                DUMPER["asm"]("phasing(%zd), count sub(%zd, %zd) = (%zd, %zd)\n", _, i, j, count[0], count[1]);
            }
            
            DUMPER["asm"]("phasing cout =(%zd, %zd)\n", count[0], count[1]);
            auto s = count[0] + count[1];
            link_scores[i] = s > 0 ? (count[0]*1.0 - count[1]) / s : 0.0;
            DUMPER["asm"]("phasing link_scores(%zd) = %f\n", i, link_scores[i]);
        }
        auto mm = std::min_element(link_scores.begin(), link_scores.end());
        DUMPER["asm"]("phasing mm = %zd\n", (size_t)(mm - link_scores.begin()));
        if (*mm < 0) {
            phasing_[mm - link_scores.begin()] = (phasing_[mm - link_scores.begin()] + 1) % 2;
        } else {
            break;
        }
        for (size_t i = 0; i < phasing_.size(); ++i) {
            DUMPER["asm"]("phasing result[%zd] = %d\n", i, phasing_[i]);
        }
    }

}

void ContigGenerator::Contig::PhaseBubbles() {

    size_t iphasing = 0;
    for (auto curr = acontigs.begin(); curr != acontigs.end(); ++curr, ++iphasing) {
        auto next = curr;
        next++;
        if (next == acontigs.end()) break;

        auto in_edge_0 = curr->second[0].back();
        auto in_edge_1 = curr->second[1].back();

        auto out_edge_0 = next->second[0].front();
        auto out_edge_1 = next->second[1].front();

        DUMPER["asm"]("start_phase_cross %s %s %s %s\n", 
            in_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            in_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            out_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            out_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str());

        auto r = PhaseCross(in_edge_0, in_edge_1, out_edge_0, out_edge_1);
        if (r < 0) {       
            DUMPER["asm"]("phase_cross %s %s %s %s\n", 
                in_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
                out_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
                in_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
                out_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str());
            phasing_[iphasing+1] = (phasing_[iphasing]+1) % 2;
        } else if (r > 0) {
            DUMPER["asm"]("phase_cross %s %s %s %s\n", 
                in_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
                out_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
                in_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
                out_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str());
            phasing_[iphasing+1] = phasing_[iphasing];
        } else {
            // pass
        }

    }

}

void ContigGenerator::Contig::PhaseHeadTail() {
    size_t first_bubble = segs_.size();
    size_t last_bubble = segs_.size();
    for (size_t i = 0; i < segs_.size(); ++i) {
        if (segs_.size() > 1) {
            if (first_bubble == segs_.size()) {
                first_bubble = i;
            }
            last_bubble = i;
        }
    }

    if (first_bubble > 0 && first_bubble < segs_.size())  {

        auto in_edge = segs_[first_bubble-1][0].back();

        auto out_edge_0 = segs_[first_bubble][0].front();
        auto out_edge_1 = segs_[first_bubble][1].front();

        DUMPER["asm"]("start_phase_head %s %s %s\n", 
            in_edge->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            out_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            out_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str());

//        PhaseHead(in_edge, out_edge_0, out_edge_1);
    }

    if (last_bubble + 1 < segs_.size()) {
        
        auto in_edge_0 = segs_[last_bubble][0].back();
        auto in_edge_1 = segs_[last_bubble][1].back();

        auto out_edge = segs_[last_bubble+1][0].front();

        DUMPER["asm"]("start_phase_tail %s %s %s\n", 
            in_edge_0->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            in_edge_1->Id().ToString(owner_->dataset_.GetStringPool()).c_str(),
            out_edge->Id().ToString(owner_->dataset_.GetStringPool()).c_str());


    }

}

std::unordered_map<int, std::unordered_map<int, int>> ComfirmVariants(const ContigGenerator::Contig::Variants &vars, double rate, size_t count);
std::unordered_map<int, std::unordered_set<int>> PairComfirmPosition(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm0, const std::unordered_map<int, std::unordered_map<int, int>> &comfirm1);
std::unordered_map<int, std::unordered_set<int>> UnionComfirmPosition(const std::unordered_map<int, std::unordered_set<int>> &pos0, const std::unordered_map<int, std::unordered_set<int>>& pos1);
std::array<size_t,2> Similarity(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm0, const std::unordered_map<int, std::unordered_map<int, int>> &comfirm1, const std::unordered_map<int, std::unordered_set<int>> &position) ;
std::array<size_t,2> Similarity(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm0, const std::unordered_map<int, std::unordered_map<int, int>> &comfirm1) ;

void PrintComfirm(const std::unordered_map<int, std::unordered_set<int>>& comfirm, const std::string& msg);
std::unordered_map<int, std::unordered_set<int>> IntersectComfirmPosition(const std::unordered_map<int, std::unordered_set<int>> &pos0, const std::unordered_map<int, std::unordered_set<int>>& pos1) ;
int ContigGenerator::Contig::PhaseCross(BaseEdge* in_edge_0, BaseEdge *in_edge_1, BaseEdge *out_edge_0, BaseEdge *out_edge_1) {
    auto in_node_0 = in_edge_0->OutNode() == in_edge_1->OutNode() ? in_edge_0->InNode() : in_edge_0->OutNode();
    auto in_node_1 = in_edge_0->OutNode() == in_edge_1->OutNode() ? in_edge_1->InNode() : in_edge_1->OutNode();

    auto out_node_0 = out_edge_0->InNode() == out_edge_1->InNode() ? out_edge_0->OutNode() : out_edge_0->InNode();
    auto out_node_1 = out_edge_0->InNode() == out_edge_1->InNode() ? out_edge_1->OutNode() : out_edge_1->InNode();

    auto in_vars = ExtendOutSnps(in_node_0, in_node_1, out_node_0, out_node_1);
    auto out_vars = ExtendInSnps(in_node_0, in_node_1, out_node_0, out_node_1);

    auto in_comfirm0 = ComfirmVariants(in_vars[0], 0.66, 1);
    auto in_comfirm1 = ComfirmVariants(in_vars[1], 0.66, 1);
    auto in_comfirm_position = PairComfirmPosition(in_comfirm0, in_comfirm1);
    PrintComfirm(in_comfirm_position, "in_pair");

    auto out_comfirm0 = ComfirmVariants(out_vars[0], 0.66, 1);
    auto out_comfirm1 = ComfirmVariants(out_vars[1], 0.66, 1);
    auto out_comfirm_position = PairComfirmPosition(out_comfirm0, out_comfirm1);
    PrintComfirm(out_comfirm_position, "out_pair");


    auto comfirm_position1 = UnionComfirmPosition(in_comfirm_position, out_comfirm_position);
    auto comfirm_position = IntersectComfirmPosition(in_comfirm_position, out_comfirm_position);
    PrintComfirm(comfirm_position1, "union");
    PrintComfirm(comfirm_position, "intersect");

    
    // test similariy
    auto r00 = Similarity(in_comfirm0, out_comfirm0, comfirm_position);
    auto r01 = Similarity(in_comfirm0, out_comfirm1, comfirm_position);
    auto r10 = Similarity(in_comfirm1, out_comfirm0, comfirm_position);
    auto r11 = Similarity(in_comfirm1, out_comfirm1, comfirm_position);
                
    DUMPER["asm"]("phase_cross_similary_intersect r00=(%zd,%zd), r01(%zd %zd), r10=(%zd,%zd), r11=(%zd,%zd)\n", 
        r00[0], r00[1], r01[0], r01[1], r10[0], r10[1], r11[0], r11[1]);

    if (r00[1] + r11[1] < r01[1] + r10[1]) {
        return 1;
    } else if (r00[1] + r11[1] > r01[1] + r10[1]) {
        return -1;
    } else {
        auto r00 = Similarity(in_comfirm0, out_comfirm0, comfirm_position1);
        auto r01 = Similarity(in_comfirm0, out_comfirm1, comfirm_position1);
        auto r10 = Similarity(in_comfirm1, out_comfirm0, comfirm_position1);
        auto r11 = Similarity(in_comfirm1, out_comfirm1, comfirm_position1);
        DUMPER["asm"]("phase_cross_similary_union r00=(%zd,%zd), r01(%zd %zd), r10=(%zd,%zd), r11=(%zd,%zd)\n", 
            r00[0], r00[1], r01[0], r01[1], r10[0], r10[1], r11[0], r11[1]);
            
        if (r00[1] + r11[1] < r01[1] + r10[1]) {
            return 1;
        } else if (r00[1] + r11[1] > r01[1] + r10[1]) {
            return -1;
        } else {
            auto r00 = Similarity(in_comfirm0, out_comfirm0);
            auto r01 = Similarity(in_comfirm0, out_comfirm1);
            auto r10 = Similarity(in_comfirm1, out_comfirm0);
            auto r11 = Similarity(in_comfirm1, out_comfirm1);
            DUMPER["asm"]("phase_cross_similary_all r00=(%zd,%zd), r01(%zd %zd), r10=(%zd,%zd), r11=(%zd,%zd)\n", 
                r00[0], r00[1], r01[0], r01[1], r10[0], r10[1], r11[0], r11[1]);
                
            if (r00[1] + r11[1] < r01[1] + r10[1]) {
                return 1;
            } else if (r00[1] + r11[1] > r01[1] + r10[1]) {
                return -1;
            } else {
                return 0;
            }
        }
    }
    
}
/*
int ContigGenerator::Contig::PhaseHead(BaseEdge* in_edge, BaseEdge *out_edge_0, BaseEdge *out_edge_1) {
    auto in_node_0 = in_edge_0->InNode();

    auto out_node_0 = out_edge_0->InNode() == out_edge_1->InNode() ? out_edge_0->OutNode() : out_edge_0->InNode();
    auto out_node_1 = out_edge_0->InNode() == out_edge_1->InNode() ? out_edge_1->OutNode() : out_edge_1->InNode();

    auto in_vars = ExtendOutSnps(in_node, nullptr, out_node_0, out_node_1);
    LOG(INFO)("ExtendOutSnps");
    auto out_vars = ExtendInSnps(in_node, nullptr, out_node_0, out_node_1);
    LOG(INFO)("ExtendInSnps");

    auto in_comfirm0 = ComfirmVariants(in_vars[0], 0.66, 2);
    auto in_comfirm1 = ComfirmVariants(in_vars[1], 0.66, 2);
    auto in_comfirm_position = PairComfirmPosition(in_comfirm0, in_comfirm1);
    PrintComfirm(in_comfirm_position, "in_pair");

    auto out_comfirm0 = ComfirmVariants(out_vars[0], 0.66, 2);
    auto out_comfirm1 = ComfirmVariants(out_vars[1], 0.66, 2);
    auto out_comfirm_position = PairComfirmPosition(out_comfirm0, out_comfirm1);
    PrintComfirm(out_comfirm_position, "out_pair");


    auto comfirm_position1 = UnionComfirmPosition(in_comfirm_position, out_comfirm_position);
    auto comfirm_position = IntersectComfirmPosition(in_comfirm_position, out_comfirm_position);
    PrintComfirm(comfirm_position1, "union");
    PrintComfirm(comfirm_position, "intersect");

    
    // test similariy
    auto r00 = Similarity(in_comfirm0, out_comfirm0, comfirm_position);
    auto r01 = Similarity(in_comfirm0, out_comfirm1, comfirm_position);
    auto r10 = Similarity(in_comfirm1, out_comfirm0, comfirm_position);
    auto r11 = Similarity(in_comfirm1, out_comfirm1, comfirm_position);
                
    DUMPER["asm"]("phase_cross_similary_new r00=(%zd,%zd), r01(%zd %zd), r10=(%zd,%zd), r11=(%zd,%zd)\n", 
        r00[0], r00[1], r01[0], r01[1], r10[0], r10[1], r11[0], r11[1]);

    if (r00[1] + r11[1] < r01[1] + r10[1]) {
        return 1;
    } else if (r00[1] + r11[1] > r01[1] + r10[1]) {
        return -1;
    } else {
        return 0;
    }
    
}
*/

// std::array<size_t,2> ContigGenerator::Contig::Similarity(const Variants &vars0, const Variants &vars1) {
   
//     auto comfirm0 = ComfirmVariants(vars0, 0.66, 2);
//     auto comfirm1 = ComfirmVariants(vars1, 0.66, 2);

//     std::array<size_t, 2> count = {0, 0};
//     for (const auto& i0 : comfirm0) {
//         auto i1 = comfirm1.find(i0.first);
//         if (i1 != comfirm1.end()) {
//             for (const auto& ii0 : i0.second) {
//                 auto ii1 = i1->second.find(ii0.first);
//                 if (ii1 != i1->second.end()) {
//                     if (ii0.second == ii1->second) {
//                         count[0]++;
//                     } else {
//                         count[1]++;
//                     }
//                 }
//             }
//         }
//     }
//     return count;
// }

std::array<size_t,2> Similarity(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm0, const std::unordered_map<int, std::unordered_map<int, int>> &comfirm1, const std::unordered_map<int, std::unordered_set<int>> &position) {
    
    std::array<size_t, 2> count = {0, 0};
    for (auto &p : position) {
        auto i0 = comfirm0.find(p.first);
        auto i1 = comfirm1.find(p.first);
        if (i0 != comfirm0.end() && i1 != comfirm1.end()) {
            for (auto & pp : p.second) {
                auto ii0 = i0->second.find(pp);
                auto ii1 = i1->second.find(pp);

                if (ii0 != i0->second.end() && ii1 != i1->second.end()) {
                    if (ii0->second == ii1->second) {
                        count[0] ++;
                    } else {
                        count[1] ++;
                    }
                }

            }
        }
    }
    return count;
}

std::array<size_t,2> Similarity(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm0, const std::unordered_map<int, std::unordered_map<int, int>> &comfirm1) {
    
    std::array<size_t, 2> count = {0, 0};
    for (auto &i0 : comfirm0) {
        auto i1 = comfirm1.find(i0.first);
        if (i1 != comfirm1.end()) {
            for (auto & ii0 : i0.second) {
                auto ii1 = i1->second.find(ii0.first);

                if (ii1 != i1->second.end()) {
                    if (ii0.second == ii1->second) {
                        count[0] ++;
                    } else {
                        count[1] ++;
                    }
                }

            }
        }
    }
    return count;
}

void AddVariants(ContigGenerator::Contig::Variants &allvars, const std::vector<ReadVariants::Variants>* vars, bool add_ctg=false) {
    if (vars == nullptr) return;

    for (auto &v : *vars) {
        if (!add_ctg && allvars.find(v.contig) == allvars.end()) continue;

        auto &ctg = allvars[v.contig];
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


void PrintVariants(const std::unordered_map<int, std::unordered_map<int, std::vector<std::array<int, 2>>>>& variants) {
    for (auto &ctg : variants) {

        DUMPER["asm"]("Variant(%d):", ctg.first);
        std::vector<int> keys(ctg.second.size());
        std::transform(ctg.second.begin(), ctg.second.end(), keys.begin(), [](const std::pair<int, std::vector<std::array<int, 2>>> &pair) { return pair.first;});

        std::sort(keys.begin(), keys.end());
        for (auto k : keys) {
            auto &v = ctg.second.find(k)->second;
            DUMPER["asm"]("%d", k);
            for (const auto &iv : v) {
                DUMPER["asm"]("-(%d,%d) ", iv[0], iv[1]);
            }
            DUMPER["asm"](" ");

        }
        DUMPER["asm"]("\n");
    }
}


void PrintComfirm(const std::unordered_map<int, std::unordered_set<int>>& comfirm, const std::string& msg) {

    for (auto &ctg : comfirm) {
        
        DUMPER["asm"]("comfirm(%s, %d):", msg.c_str(), ctg.first);
        std::vector<int> keys(ctg.second.begin(), ctg.second.end());

        std::sort(keys.begin(), keys.end());
        for (auto k : keys) {
            DUMPER["asm"]("%d, ", k);
        }
        DUMPER["asm"]("\n");
    }
}

std::unordered_map<int, std::unordered_map<int, int>> ComfirmVariants(const ContigGenerator::Contig::Variants &vars, double rate, size_t count) {
    auto important = [rate, count](std::vector<std::array<int, 2>> d) {

        std::sort(d.begin(), d.end(), [](const std::array<int,2> &a, const std::array<int,2> &b) {
            return (a[0] >= 0 && b[0] >= 0) ? a[1] > b[1] : a[0] > b[0];
        });

        int sum = std::accumulate(d.begin(), d.end(), 0, [](int a, const std::array<int,2>& b) {
            return a + (b[0] >= 0 ? b[1] : 0);        // remove seq error (-1)
        });

        assert(d.size() > 0);
        if (d[0][1] >= sum*rate && d[0][1] >= (int)count) { // TODO PARAM NAME
            return d[0][0];
        } else {
            return -1;
        }

    };

    auto get_comfirm = [important](const ContigGenerator::Contig::Variants &vars) {
        std::unordered_map<int, std::unordered_map<int, int>> vs;    
        for (auto &i : vars) {
            auto &cv = vs[i.first];

            for (auto &j : i.second) {
                auto b = important(j.second);
                if (b >= 0) {
                    cv[j.first] = b;
                }
            }

        }

        return vs;
    };

    return get_comfirm(vars);
}

std::unordered_map<int, std::unordered_set<int>> PairComfirmPosition(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm0, const std::unordered_map<int, std::unordered_map<int, int>> &comfirm1) {
    std::unordered_map<int, std::unordered_set<int>> positions;

    for (const auto &i0 : comfirm0) {
        auto i1 = comfirm1.find(i0.first);
        if (i1 != comfirm1.end()) {
            for (const auto& ii0 : i0.second) {
                auto ii1 = i1->second.find(ii0.first);
                if (ii1 != i1->second.end() && ii1->second != ii0.second) {
                    positions[i0.first].insert(ii0.first);
                }
            }
        }
    }
    return positions;
}

std::unordered_map<int, std::unordered_set<int>> UnionComfirmPosition(const std::unordered_map<int, std::unordered_set<int>> &pos0, const std::unordered_map<int, std::unordered_set<int>>& pos1) {
    std::unordered_map<int, std::unordered_set<int>> pos(pos0);

    for (const auto &i : pos1) {
        pos[i.first].insert(i.second.begin(), i.second.end());
    }
    return pos;
}

std::unordered_map<int, std::unordered_set<int>> IntersectComfirmPosition(const std::unordered_map<int, std::unordered_set<int>> &pos0, const std::unordered_map<int, std::unordered_set<int>>& pos1) {
    std::unordered_map<int, std::unordered_set<int>> pos;

    for (const auto &i0 : pos0) {
        auto i1 = pos1.find(i0.first);
        if (i1 != pos1.end()) {
            for (const auto &ii0 : i0.second) {
                auto ii1 = i1->second.find(ii0);
                if (ii1 != i1->second.end()) {
                    pos[i0.first].insert(ii0);
                }
            }
        } 
    }
    return pos;
}


std::array<size_t,2> CompareVariants(const std::unordered_map<int, std::unordered_map<int, int>> &comfirm, const std::vector<ReadVariants::Variants>* vars) {
    std::array<size_t, 2> count = {0, 0};
    if (vars == nullptr) return count;

    for (auto &v : *vars) {
        auto it = comfirm.find(v.contig);
        if (it != comfirm.end()) {

            for (auto &i : v.vars) {
                auto it1 = it->second.find(i.first);
                if (it1 != it->second.end()) {
                    if (it1->second == i.second) {
                        count[0] ++;
                    } else {
                        count[1] --;
                    }
                }
            }
        }
    }
    
    return count;
}

auto ContigGenerator::Contig::ExtendOutSnps(BaseNode* in_node_0, BaseNode* in_node_1, BaseNode* out_node_0, BaseNode* out_node_1) -> std::array<Variants,2> {

    auto & asmdata = owner_->dataset_;

    std::array<Variants,2> variants;
    auto reads_0 = CollectOutSnps(variants[0], in_node_0, in_node_1);
    auto reads_1 = CollectOutSnps(variants[1], in_node_1, in_node_0);

    auto comfirm_0 = ComfirmVariants(variants[0], 0.66, 1);
    auto comfirm_1 = ComfirmVariants(variants[1], 0.66, 1);


    std::unordered_set<Seq::Id> done;
    done.insert(reads_0.begin(), reads_0.end());
    done.insert(reads_1.begin(), reads_1.end());

    Seq::Id in_id_0 = Seq::EndIdToId(in_node_0->Id().MainNode());
    auto overlap_reads_0 = asmdata.GetOverlapReads(in_id_0);

    Seq::Id in_id_1 = Seq::EndIdToId(in_node_1->Id().MainNode());
    auto overlap_reads_1 = asmdata.GetOverlapReads(in_id_1);

    auto rvs = owner_->dataset_.GetReadVariants();
    assert(rvs != nullptr);


    
    for (auto r : overlap_reads_0) {
        if (done.find(r) != done.end()) continue;
        auto vs = rvs->GetVariants(r);
        if (vs != nullptr) {

            auto diff0 = CompareVariants(comfirm_0, vs);
            auto diff1 = CompareVariants(comfirm_1, vs);

            if (diff0[1] < diff1[1] && diff0[0] > diff0[1]) {
                DUMPER["asm"]("phase_cross add out ol %s to 0 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[0], vs);
            } else if (diff1[1] < diff0[1] && diff1[0] > diff1[1]) {
                DUMPER["asm"]("phase_cross add out ol %s to 1 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[1], vs);
            }
        }

        done.insert(r);
    }
    
    for (auto r : overlap_reads_1) {
        if (done.find(r) != done.end()) continue;
        auto vs = rvs->GetVariants(r);
        if (vs != nullptr) {

            auto diff0 = CompareVariants(comfirm_0, vs);
            auto diff1 = CompareVariants(comfirm_1, vs);

            if (diff0[1] < diff1[1] && diff0[0] > diff0[1]) {
                DUMPER["asm"]("phase_cross add out ol %s to 0 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[0], vs);
            } else if (diff1[1] < diff0[1] && diff1[0] > diff1[1]) {
                DUMPER["asm"]("phase_cross add out ol %s to 1 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[1], vs);
            }
        }

        done.insert(r);
    }
    
    PrintVariants(variants[0]);
    PrintVariants(variants[1]);
    return variants;
}


std::unordered_set<Seq::Id> ContigGenerator::Contig::CollectOutSnps(Variants& vars, BaseNode* node, BaseNode* altnode) {
    
    auto & asmdata = owner_->dataset_;


    std::unordered_set<Seq::Id> reads;

    Seq::Id id = Seq::EndIdToId(node->Id().MainNode());
    int end = Seq::End(node->Id().MainNode());
    std::unordered_set<const Overlap*> ols = asmdata.GetBackOverlaps(id, end);
    for (auto ol : ols) {
        auto qread = ol->GetOtherRead(id);
        auto altol = asmdata.QueryOverlap(qread.id, Seq::EndIdToId(altnode->Id().MainNode()));
        DUMPER["asm"]("phase_cross add out ol %s (%d)\n", asmdata.QueryNameById(qread.id).c_str(), altol==nullptr);
        if (altol == nullptr)  {
            reads.insert(qread.id);
        } 
    }

    auto pif = asmdata.GetInconsistentOverlaps();
    assert(pif != nullptr);

    auto p_node = node->InDegree() == 1 ? node->InNode(0) : nullptr;
    Seq::Id p_id = p_node == nullptr ? -1 : Seq::EndIdToId(p_node->Id().MainNode());
    
    auto p_altnode = altnode->InDegree() == 1 ? altnode->InNode(0) : nullptr;

    auto alt_inconsist = pif->Get(Seq::EndIdToId(altnode->Id().MainNode()));
    for (auto r : alt_inconsist) {
        if (asmdata.QueryOverlap(id, r) != nullptr || (p_id != -1 && asmdata.QueryOverlap(p_id, r) != nullptr)) {
            DUMPER["asm"]("phase_cross add out0 alt %s\n", asmdata.QueryNameById(r).c_str());
            reads.insert(r);
        }
    }

    if (p_altnode != nullptr) {

        auto p_alt_inconsist = pif->Get(Seq::EndIdToId(p_altnode->Id().MainNode()));
        for (auto r : p_alt_inconsist) {
            if (asmdata.QueryOverlap(id, r) != nullptr || (p_id != -1 && asmdata.QueryOverlap(p_id, r) != nullptr)) {
                DUMPER["asm"]("phase_cross add out1 alt %s\n", asmdata.QueryNameById(r).c_str());
                reads.insert(r);
            }
        }
    }

    auto rvs = owner_->dataset_.GetReadVariants();
    AddVariants(vars, rvs->GetVariants(id), true);
    for (auto r : reads) {
        AddVariants(vars, rvs->GetVariants(r));
    }
    PrintVariants(vars);
    return reads;
}

auto ContigGenerator::Contig::ExtendInSnps(BaseNode* in_node_0, BaseNode* in_node_1, BaseNode* out_node_0, BaseNode* out_node_1) -> std::array<Variants,2> {

    auto & asmdata = owner_->dataset_;

    std::array<Variants,2> variants;
    auto reads_0 = CollectInSnps(variants[0], out_node_0, out_node_1);
    auto reads_1 = CollectInSnps(variants[1], out_node_1, out_node_0);

    auto comfirm_0 = ComfirmVariants(variants[0], 0.66, 1);
    auto comfirm_1 = ComfirmVariants(variants[1], 0.66, 1);


    std::unordered_set<Seq::Id> done;
    done.insert(reads_0.begin(), reads_0.end());
    done.insert(reads_1.begin(), reads_1.end());

    Seq::Id out_id_0 = Seq::EndIdToId(out_node_0->Id().MainNode());
    auto overlap_reads_0 = asmdata.GetOverlapReads(out_id_0);

    Seq::Id out_id_1 = Seq::EndIdToId(out_node_1->Id().MainNode());
    auto overlap_reads_1 = asmdata.GetOverlapReads(out_id_1);

    auto rvs = owner_->dataset_.GetReadVariants();
    assert(rvs != nullptr);


    for (auto r : overlap_reads_0) {
        if (done.find(r) != done.end()) continue;
        auto vs = rvs->GetVariants(r);
        if (vs != nullptr) {

            auto diff0 = CompareVariants(comfirm_0, vs);
            auto diff1 = CompareVariants(comfirm_1, vs);

            if (diff0[1] < diff1[1] && diff0[0] > diff0[1]) {
                DUMPER["asm"]("phase_cross add in ol %s to 0 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[0], vs);
            } else if (diff1[1] < diff0[1] && diff1[0] > diff1[1]) {
                DUMPER["asm"]("phase_cross add in ol %s to 1 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[1], vs);
            }
        }

        done.insert(r);
    }
    
    for (auto r : overlap_reads_1) {
        if (done.find(r) != done.end()) continue;
        auto vs = rvs->GetVariants(r);
        if (vs != nullptr) {

            auto diff0 = CompareVariants(comfirm_0, vs);
            auto diff1 = CompareVariants(comfirm_1, vs);

            if (diff0[1] < diff1[1] && diff0[0] > diff0[1]) {
                DUMPER["asm"]("phase_cross add in ol %s to 0 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[0], vs);
            } else if (diff1[1] < diff0[1] && diff1[0] > diff1[1]) {
                DUMPER["asm"]("phase_cross add in ol %s to 1 (%zd, %zd, %zd, %zd)\n", asmdata.QueryNameById(r).c_str(), diff0[0], diff0[1], diff1[0], diff1[1]);
                AddVariants(variants[1], vs);
            }
        }

        done.insert(r);
    }
    
    PrintVariants(variants[0]);
    PrintVariants(variants[1]);
    return variants;
}


// void ContigGenerator::Contig::CollectInSnps(Variants& vars, BaseNode* node, BaseNode* altnode) {
    
//     auto & asmdata = owner_->dataset_;
//     auto rvs = owner_->dataset_.GetReadVariants();

//     Seq::Id rid = Seq::EndIdToId(node->Id().MainNode());
//     int rend = Seq::End(node->Id().MainNode());
//     std::unordered_set<const Overlap*> ols = asmdata.GetExtendOverlaps(rid, rend);

//     AddVariants(vars, rvs->GetVariants(rid));
//     for (auto ol : ols) {
//         auto qread = ol->GetOtherRead(rid);
//         auto altol = asmdata.QueryOverlap(qread.id, Seq::EndIdToId(altnode->Id().MainNode()));
//         DUMPER["asm"]("phase_cross add in ol %s (%d)\n", asmdata.QueryNameById(qread.id).c_str(), altol==nullptr);
//         if (altol == nullptr)  {
//             reads.insert(qread.id);
//             AddVariants(vars, rvs->GetVariants(qread.id));
//         } 
//     }

//     //
//     auto pif = asmdata.GetInconsistentOverlaps();
//     assert(pif != nullptr);
//     auto inconsist = pif->Get(Seq::EndIdToId(node->Id().MainNode()));
//     auto alt_inconsist = pif->Get(Seq::EndIdToId(altnode->Id().MainNode()));
//     for (auto r : alt_inconsist) {
//         if (inconsist.find(r) == inconsist.end()) {
//             DUMPER["asm"]("phase_cross add in alt %s\n", asmdata.QueryNameById(r).c_str());
//             AddVariants(vars, rvs->GetVariants(r));
//         }
//     }
//     PrintVariants(vars);
// }

std::unordered_set<Seq::Id> ContigGenerator::Contig::CollectInSnps(Variants& vars, BaseNode* node, BaseNode* altnode) {
    
    auto & asmdata = owner_->dataset_;

    std::unordered_set<Seq::Id> reads;

    Seq::Id id = Seq::EndIdToId(node->Id().MainNode());
    int end = Seq::End(node->Id().MainNode());

    std::unordered_set<const Overlap*> ols = asmdata.GetExtendOverlaps(id, end);
    for (auto ol : ols) {
        auto qread = ol->GetOtherRead(id);
        auto altol = asmdata.QueryOverlap(qread.id, Seq::EndIdToId(altnode->Id().MainNode()));
        DUMPER["asm"]("phase_cross add in ol %s (%d)\n", asmdata.QueryNameById(qread.id).c_str(), altol==nullptr);
        if (altol == nullptr)  {
            reads.insert(qread.id);
        } 
    }

    auto pif = asmdata.GetInconsistentOverlaps();
    assert(pif != nullptr);

    auto p_node = node->OutDegree() == 1 ? node->OutNode(0) : nullptr;
    Seq::Id p_id = p_node == nullptr ? -1 : Seq::EndIdToId(p_node->Id().MainNode());

    auto alt_inconsist = pif->Get(Seq::EndIdToId(altnode->Id().MainNode()));
    for (auto r : alt_inconsist) {
        if (asmdata.QueryOverlap(id, r) != nullptr || (p_id != -1 && asmdata.QueryOverlap(p_id, r) != nullptr)) {
            DUMPER["asm"]("phase_cross add in0 alt %s\n", asmdata.QueryNameById(r).c_str());
            reads.insert(r);
        }
    }

    auto p_altnode = altnode->OutDegree() == 1 ? altnode->OutNode(0) : nullptr;
    if (p_altnode != nullptr) {

        auto p_alt_inconsist = pif->Get(Seq::EndIdToId(p_altnode->Id().MainNode()));
        for (auto r : p_alt_inconsist) {
            if (asmdata.QueryOverlap(id, r) != nullptr || (p_id != -1 && asmdata.QueryOverlap(p_id, r) != nullptr)) {
                DUMPER["asm"]("phase_cross add in1 alt %s\n", asmdata.QueryNameById(r).c_str());
                reads.insert(r);
            }
        }
    }

    auto rvs = owner_->dataset_.GetReadVariants();
    AddVariants(vars, rvs->GetVariants(id), true);
    for (auto r : reads) {
        AddVariants(vars, rvs->GetVariants(r));
    }
    PrintVariants(vars);
    return reads;
}

std::list<BaseEdge*> ContigGenerator::Contig::ConstructPrimaryPath(uint8_t hap) {
    assert(hap == 0 || hap == 1);
    
    std::array<size_t,2> nodes = {0, 0};
    // calc cov
    std::vector<size_t> covs; 
    for (auto &s : segs_) {
        if (s.size() >= 2) {
            nodes[1] += s[1].size();
            for (auto p : s[0]) {
                const auto &rinfo = owner_->dataset_.GetReadInfo(p->OutNode()->ReadId());
                covs.push_back(rinfo.coverage[1]);
            }
            for (auto p : s[1]) {
                const auto &rinfo = owner_->dataset_.GetReadInfo(p->OutNode()->ReadId());
                covs.push_back(rinfo.coverage[1]);
            }
        } else if (s.size() >= 1) {
            nodes[0] += s[0].size();
        }
    }

    bool extend = hap == 0 || (hap == 1 && nodes[1]*10 >= nodes[0]);

    std::sort(covs.begin(), covs.end());
    size_t hapcov = std::accumulate(covs.begin()+covs.size()/4, covs.begin()+covs.size()*3/4, 0);
    if (hapcov > 0) {
        assert(covs.size() > 0);
        hapcov = hapcov*2 / covs.size();
    }


    std::list<BaseEdge*> path;
    size_t ih = 0;
    for (auto &s : segs_) {
        if (hap == 1 && s.size() == 1 && (s == segs_.front() || s == segs_.back()) && extend) {
            size_t cov = 0.0;
            for (auto p : s[0]) {
                const auto &rinfo = owner_->dataset_.GetReadInfo(p->OutNode()->ReadId());
                cov += rinfo.coverage[1];
            }
            if (cov > 0) cov /= s[0].size();
            if (cov < hapcov * 1.5 || s[0].size() >= covs.size() / 2) {
                DUMPER["asm"]("skip copying contig %d %s. (%zd < %zd*1.5)\n", id_, (s == segs_.front() ? "head" : "tail"), cov, hapcov);
                continue;
            }
        }

        std::list<BaseEdge*> *p = nullptr;
        if (s.size() == 1) {
            if (extend) p = &(s[0]);
        } else {
            size_t i = phasing_[ih++] == 0 ? hap : (hap+1) % 2; 
            p = &(s[i]);
        }
        if (p != nullptr) {
            if (path.size() > 0 && p->size() > 0 && path.back()->OutNode() != p->front()->InNode()) {
                auto gap = owner_->string_graph_.GetPath(path.back()->OutNode(), p->front()->InNode(), 100);
                // if (gap.size() == 0) // TODO
                path.insert(path.end(), gap.begin(), gap.end());
            }
            path.insert(path.end(), p->begin(), p->end());
        }
    }
    return path;
}

std::string ContigGenerator::Contig::MainName() const {
    std::ostringstream oss;
    oss << "ctg" << id_;
    return oss.str();
}

std::string ContigGenerator::Contig::SubName(size_t ibubble, size_t ipath) const {
    std::ostringstream oss;
    oss << MainName() << "-" << ibubble << "-" << ipath;
    return oss.str();
}

bool ContigGenerator::Contig::IsCircular() const {
    return pcontig.front()->InNode() == pcontig.back()->OutNode();
}

std::string ContigGenerator::Contig::Description(size_t len) const {
    std::ostringstream oss;
    oss << (IsPrimary() ? "" : pri_->MainName()) << ":"
        << (IsCircular() ? "circular" : "linear") << ":"
        << "length=" << len;
                    
    return oss.str();
}

bool ContigGenerator::Contig::IsCovered(const std::vector<std::list<BaseEdge*>>& paths) const {
    const int N = 3;
    std::vector<Seq::Id> lids, rids;
    for (auto &p : paths) {
        if (p.front()->InNode()->OutDegree() > 1) {
            BaseEdge* curr = p.front();
            for (size_t i = 0; i < N; ++i) {
                lids.push_back(curr->InRead()->id);
                if (curr->InNode()->InDegree() == 1) {
                    curr = curr->InNode()->InEdge<BaseEdge>(0);
                } else {
                    break;
                }
            }
        }

        if (p.back()->OutNode()->InDegree() > 1) {
            BaseEdge* curr = p.back();
            for (size_t i = 0; i < N; ++i) {
                rids.push_back(curr->OutRead()->id);
                if (curr->OutNode()->OutDegree() == 1) {
                    curr = curr->OutNode()->OutEdge<BaseEdge>(0);
                } else {
                    break;
                }
            }

        }
    }


    for (auto a : alts_) {
        size_t ir = 0;
        for (; ir < rids.size(); ++ir) {
            if (a->vreads.find(rids[ir]) != a->vreads.end()) {
                break;
            }
        }

        size_t il = 0;
        for (; il < lids.size(); ++il) {
            if (a->vreads.find(lids[il]) != a->vreads.end()) {
                break;
            }
        }

        DUMPER["asm"]("xcontig: %s -> %s: %zd/%zd, %zd/%zd\n", MainName().c_str(), a->MainName().c_str(), ir, rids.size(), il, lids.size());

        if (il < lids.size() && ir < rids.size()) {
            return a->pcontig.size() >= paths[0].size() / 2;
        } 
        
    }

    if (pri_ != nullptr) {
       
    size_t ir = 0;
    for (; ir < rids.size(); ++ir) {
        if (pri_->vreads.find(rids[ir]) != pri_->vreads.end()) {
            break;
        }
    }

    size_t il = 0;
    for (; il < lids.size(); ++il) {
        if (pri_->vreads.find(lids[il]) != pri_->vreads.end()) {
            break;
        }
    }
    DUMPER["asm"]("xcontig: %s -> %s: %zd/%zd, %zd/%zd\n", MainName().c_str(), pri_->MainName().c_str(), ir, rids.size(), il, lids.size());

    if (il < lids.size() && ir < rids.size()) {
        return pri_->pcontig.size() >= paths[0].size() / 2;
    } 
         
    }
    
    return false;
}
}