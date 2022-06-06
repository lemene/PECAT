#include "overlap_assemble.hpp"

#include <cassert>
#include <algorithm>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include <sstream>
#include <iostream>
#include <atomic>
#include <edlib.h>

#include "../file_io.hpp"
#include "../utility.hpp"


#include "../utils/logger.hpp"
#include "../simple_align.hpp"

namespace fsa {

OverlapAssemble::OverlapAssemble() : string_graph_(dataset_), path_graph_(dataset_) {
}

OverlapAssemble::~OverlapAssemble() {
}


ArgumentParser OverlapAssemble::GetArgumentParser() {
    ArgumentParser ap("fsa_ol_assemble", "Assemble overlaps", "1.0");
    opts_.SetArguments(ap);
    
    return ap;
}

void OverlapAssemble::CheckArguments() {
    opts_.CheckArguments();
    //DUMPER.SetLevel(opts_.dump);
    DUMPER.SetLevel(1);
    DUMPER.SetDirectory(opts_.output_directory);
}


void OverlapAssemble::Running() {
    // if read_file is provided, overlap-loading step can be accelerated

    if (!opts_.skip_purge) {
        dataset_.Load();

        dataset_.Purge();

    } else {
        dataset_.LoadPurged();

    }

    CreateStringGraph();

    CreatePathGraph();

    SaveGraph();

    if (dataset_.GetReadStore().Size() > 0) {
        SaveContigs();
    }
}

void OverlapAssemble::CreateStringGraph() {
    LOG(INFO)("Create string graph");
    
    string_graph_.AddOverlaps(dataset_.GetOverlapStore());

    LOG(INFO)("Simplify string graph");
    string_graph_.Simplify(opts_.reduction0, opts_.reducer0);
    
    LOG(INFO)("Identify paths from string graph");
    string_graph_.IdentifySimplePaths();
}

void OverlapAssemble::CreatePathGraph() {
    LOG(INFO)("Create path graph");
    path_graph_.BuildFrom(string_graph_);
    
    LOG(INFO)("Simplify path graph");
    path_graph_.Simplify(opts_.reduction1);
    
    LOG(INFO)("Identify paths from path graph");
    path_graph_.IdentifyPaths(opts_.select_branch);
}


void OverlapAssemble::SaveContigs() {
    LOG(INFO)("Save Contigs");

    std::vector<Contig> contigs;

    const auto& paths = path_graph_.GetPaths();

    for (size_t i = 0; i < paths.size(); i += 2) {
        contigs.push_back(Contig(i/2, paths[i], string_graph_));
        contigs.back().PhasedReads(dataset_.GetInconsistentOverlaps());
    }
    std::sort(contigs.begin(), contigs.end(), [](const Contig& a, const Contig &b) { return a.Weight() > b.Weight(); });
    LOG(INFO)("Contig size: %zd", contigs.size());

    // classify contigs
    for (size_t i=1; i<contigs.size(); ++i) {

        size_t threshold = opts_.max_trivial_length;
        auto &cluster = path_graph_.clusters_[path_graph_.edge2cluster_[contigs[i].path_->front()]];
        DUMPER["asm"]("%zd: Contig(%d,%d): Cluster(%zd): size %zd,%zd", i, contigs[i].id_, contigs[i].Weight(), path_graph_.edge2cluster_[contigs[i].path_->front()], cluster.Size(), cluster.Length());
        if (cluster.Length() < threshold) {
            contigs[i].SetTrivial();
            DUMPER["asm"]("Cluster SetTrivial0 %zd", contigs[i].id_);
        } else {
            if (contigs[i].Length() < threshold) {
                if (cluster.nontrivial.find(contigs[i].path_->front()) == cluster.nontrivial.end()) {
                    contigs[i].SetTrivial();
                    DUMPER["asm"]("Cluster SetTrivial1 %zd", contigs[i].id_);
                }
            }
        }
            

        if (!opts_.phased.empty()) {
            int pri = -1;
            for (size_t j = 0; j < i; ++j) {
                if (contigs[j].IsTrivial()) continue;
                if (contigs[i].IsDiploid(contigs[j], opts_)) {
                    if (contigs[j].IsPrimary()) {
                        if (pri < 0) pri = j;
                    } else {
                        pri = j;
                        break;
                    }
                }
            }
                    
            if (pri >= 0) {
                contigs[i].SetPrimary(contigs[pri]);
                contigs[pri].AddAlternate(contigs[i]);
            }
            
        }
    }

    std::mutex mutex;
    GzFileWriter fcontig_seqs(OutputPath("primary.fasta"));
    GzFileWriter fcontig_tiles(OutputPath("primary_tiles"));
    GzFileWriter fcontig_seqs1(OutputPath("alternate.fasta"));
    GzFileWriter fcontig_tiles1(OutputPath("alternate_tiles"));
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
    auto work_func = [&contigs, &index, this, combine_func](size_t id) {
        std::ostringstream oss_pri;
        std::ostringstream oss_pri_tiles;
        std::ostringstream oss_alt;
        std::ostringstream oss_alt_tiles;
        std::ostringstream oss_oth;
        std::ostringstream oss_oth_tiles;

        for (size_t curr = index.fetch_add(1); curr < contigs.size(); curr = index.fetch_add(1)) {
            if (contigs[curr].IsTrivial())  {
                if(this->ConstructContig(contigs[curr].pcontig).size() < opts_.max_trivial_length) {
                    contigs[curr].Save(oss_oth, *this, opts_.min_contig_length);
                    contigs[curr].SaveTiles(oss_oth_tiles, string_graph_);
                    contigs[curr].SaveBubbles(oss_oth, oss_oth_tiles, oss_oth, oss_oth_tiles, *this);
                    continue;
                }
            }
            if (contigs[curr].IsPrimary()) {
                contigs[curr].Save(oss_pri, *this, opts_.min_contig_length);
                contigs[curr].SaveTiles(oss_pri_tiles, string_graph_);
                contigs[curr].SaveBubbles(oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles, *this);
            } else {
                if (contigs[curr].GetPrimary()->IsPrimary()) {
                    contigs[curr].Save(oss_alt, *this, opts_.min_contig_length);
                    contigs[curr].SaveTiles(oss_alt_tiles, string_graph_);
                    contigs[curr].SaveBubbles(oss_pri, oss_pri_tiles, oss_oth, oss_oth_tiles, *this);
                } else {
                    if ( // !contigs[curr].GetPrimary()->GetPrimary()->IsPrimary() ||
                        contigs[curr].IsDiploid(*contigs[curr].GetPrimary()->GetPrimary(), opts_) ||
                        contigs[curr].IsOverlapped(*contigs[curr].GetPrimary()->GetPrimary(), opts_)) {

                        contigs[curr].Save(oss_oth, *this, opts_.min_contig_length);
                        contigs[curr].SaveTiles(oss_oth_tiles, string_graph_);
                        contigs[curr].SaveBubbles(oss_oth, oss_oth_tiles, oss_oth, oss_oth_tiles, *this);
                    } else {
                        contigs[curr].Save(oss_pri, *this, opts_.min_contig_length);
                        contigs[curr].SaveTiles(oss_pri_tiles, string_graph_);
                        contigs[curr].SaveBubbles(oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles, *this);
                    }
                }
            }
        }
        combine_func(oss_pri, oss_pri_tiles, oss_alt, oss_alt_tiles, oss_oth, oss_oth_tiles);
    };

    MultiThreadRun((size_t)opts_.thread_size, work_func);
    

}

std::string OverlapAssemble::ConstructContigS(const std::list<BaseEdge*> &contig) {
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

std::string OverlapAssemble::ConstructContigStraight(const std::list<BaseEdge*> &contig) {
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


std::string OverlapAssemble::ConstructContig(const std::list<BaseEdge*> &contig) {
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


std::vector<Seq::Tile> OverlapAssemble::EdgesToTiles(const std::vector<BaseEdge*> &contig) {
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

std::string OverlapAssemble::EdgeToSeq(const BaseEdge *e) {

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

void OverlapAssemble::SaveGraph() {
    
    LOG(INFO)("Save Graph");
    string_graph_.SaveEdges(OutputPath("graph_edges.gz"));
    path_graph_.SaveEdges(OutputPath("graph_paths.gz"));
    path_graph_.SavePaths(OutputPath("contig_paths.gz"));
}


double OverlapAssemble::ComputeSequenceSimilarity(const std::string &qseq, const std::string &tseq) {
    auto r = edlibAlign(qseq.c_str(), qseq.size(), tseq.c_str(), tseq.size(),edlibDefaultAlignConfig());
    double identity = r.status == EDLIB_STATUS_OK ? 1 - r.editDistance *1.0 / qseq.size() : 0;


    edlibFreeAlignResult(r);
    return identity;
}

OverlapAssemble::Contig::Contig(size_t id, const std::list<PathEdge*> &path, StringGraph &sg) : id_(id), path_(&path) {
    
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

        DUMPER["asm"]("ctgswap: %zd, %d\n", id_, p->Type().c_str());
        std::vector<std::list<BaseEdge*>> ctgs;
        for (size_t i = 0; i < p->SimplePathSize(); ++i) {
            ctgs.push_back(p->GetSimplePath(i));
        }

        if (ctgs.size() > 1) {
            if ( p->IsType("bubble") || 
                (p == path.front() && p->IsType("semi") && static_cast<SemiBubbleEdge*>(p)->Branchs() == 2) ||
                (p == path.back()  && p->IsType("semi") && static_cast<SemiBubbleEdge*>(p)->Branchs() == 1)) {

                DUMPER["asm"]("ctgswap -- : %zd, %d\n", id_, p->Type().c_str());
                auto mx = std::max_element(ctgs.begin(), ctgs.end(), [path_len](const std::list<BaseEdge*>& a, const std::list<BaseEdge*> &b) {
                    return path_len(a) < path_len(b);
                });
                if (ctgs.begin() != mx) {
                    auto t = ctgs[0];
                    ctgs[0] = *mx;
                    *mx = t;
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
        if (ctgs.size() > 1) {
            acontigs.push_back(std::make_pair(p, std::move(ctgs)));
        }
    }
    weight_ = path_len(pcontig);
}


static bool ValidateDipolid(const std::list<BaseEdge*> &path0, const std::list<BaseEdge*> &path1, PhaseInfoFile &ignored, ReadStore &rs) {
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
    return count*1.0 / path0.size() > 0.5;
}

// c is primary
bool OverlapAssemble::Contig::IsDiploid(const Contig& c, const AsmOptions& opts) const {
    size_t count = 0;
    for (auto r : reads) {
        if (c.vreads.find(r) != c.vreads.end()) {
            count ++;
        } else {
        }
    }
    auto rr = count*1.0 / reads.size() >= opts.diploid_rate && (int)count >= opts.diploid_count || count*1.0 / reads.size() >= 0.45;
    
    DUMPER["asm"]("contig: %zd -> %zd: %zd/%zd/%zd, %f, %d\n", id_, c.id_, count, reads.size(), c.reads.size(), count*1.0 / reads.size(), rr);
    return count*1.0 / reads.size() >= opts.diploid_rate && (int)count >= opts.diploid_count ||  count*1.0 / reads.size() >= 1.0;

}


// c is primary
bool OverlapAssemble::Contig::IsOverlapped(const Contig& c, const AsmOptions& opts) const {
    size_t count = 0;
    for (auto r : vreads) {
        if (c.vreads.find(r) != c.vreads.end()) {
            count ++;
        } 
    }
    
    DUMPER["asm"]("overlaped: %zd -> %zd: %zd/%zd\n", id_, c.id_, count, vreads.size());
    return count > 0 && count*1.0 / vreads.size() >= 0.5; // TODO Parameterization 

}

void OverlapAssemble::Contig::PhasedReads(PhaseInfoFile *phased) {
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

void OverlapAssemble::Contig::Save(std::ostream& os,  OverlapAssemble& ass, int min_contig_length) {
    std::vector<std::string> seqs = {ass.ConstructContig(pcontig)};
    assert(seqs.size() >= 1);
    if (seqs[0].length() >= (size_t)min_contig_length) {
        os << ">" << MainName() << " " << Description(seqs[0].size()) << "\n" << seqs[0] << "\n";
    }
}


void OverlapAssemble::Contig::SaveTiles(std::ostream& os, const StringGraph& sg) {
    for (auto e : pcontig) {
        os << MainName() << " edge=" <<  sg.IdToString(e->InNode()->Id()) << 
                "~" << sg.IdToString(e->OutNode()->Id()) << "\n";
    }
}


void OverlapAssemble::Contig::SaveBubbles(std::ostream &fctg,  std::ostream &ftile, OverlapAssemble& ass) {
    int ibubble = 1;

    for (const auto &bubble : acontigs) {
        const auto &paths = bubble.second;
        assert(paths.size() > 0);
        SaveBubbles(fctg, ftile, ass, ibubble, paths);

        ibubble++;
    }
}


void OverlapAssemble::Contig::SaveBubbles(std::ostream &fctg,  std::ostream &ftile, OverlapAssemble& ass, int id, const std::vector<std::list<BaseEdge*>>& paths) {
    assert(paths.size() > 0);

    for (size_t i = 1; i < paths.size(); ++i) {
        std::string seq = ass.ConstructContigS(paths[i]);
        fctg << ">" << SubName(id, i) << " length=" << seq.size() << "\n" 
             << seq << "\n";

        for (auto p : paths[i]) {
            ftile << SubName(id, i) <<
                    " edge=" <<  ass.string_graph_.IdToString(p->InNode()->Id()) << 
                    "~" << ass.string_graph_.IdToString(p->OutNode()->Id()) << "\n";

        }
    }
}

void OverlapAssemble::Contig::SaveBubbles(std::ostream& fctg0, std::ostream& ftile0, std::ostream& fctg1, std::ostream& ftile1, OverlapAssemble& ass) {
    int ibubble = 1;

    for (const auto &bubble : acontigs) {
        const auto &paths = bubble.second;
        assert(paths.size() > 0);
        
        bool is_covered = IsCovered(paths);
        std::ostream &fctg =  !is_covered ? fctg0  : fctg1;
        std::ostream &ftile = !is_covered ? ftile0 : ftile1;

        SaveBubbles(fctg, ftile, ass, ibubble, paths);

        ibubble++;
    }
}

std::string OverlapAssemble::Contig::MainName() const {
    std::ostringstream oss;
    oss << "ctg" << id_;
    return oss.str();
}

std::string OverlapAssemble::Contig::SubName(size_t ibubble, size_t ipath) const {
    std::ostringstream oss;
    oss << MainName() << "-" << ibubble << "-" << ipath;
    return oss.str();
}

bool OverlapAssemble::Contig::IsCircular() const {
    return pcontig.front()->InNode() == pcontig.back()->OutNode();
}

std::string OverlapAssemble::Contig::Description(size_t len) const {
    std::ostringstream oss;
    oss << (IsPrimary() ? "" : pri_->MainName()) << ":"
        << (IsCircular() ? "circular" : "linear") << ":"
        << "length=" << len;
                    
    return oss.str();
}

bool OverlapAssemble::Contig::IsCovered(const std::vector<std::list<BaseEdge*>>& paths) const {
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
    return false;
}

} // namespace fsa {
