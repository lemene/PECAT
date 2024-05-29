#include "alignment_graph.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdint.h>
#include <fstream>
#include <unordered_set>

#include <algorithm>
#include <numeric>
#include "aligner.hpp"
#include "./utils/logger.hpp"
#include "../read_store.hpp"
#include "crr_dataset.hpp"
#include "crr_options.hpp"

namespace fsa {

DnaSerialTable2 AlignmentGraph::Base2Num;

AlignmentGraph::AlignmentGraph(const CrrOptions& opts, const CrrDataset& ds) 
 : sopts_(opts), dataset_(ds) {
}

void AlignmentGraph::SetParameter(const std::string &name, const std::string &opts) {
    if (name == "score") {
        ParseScoreParamter(opts);
    } else {
        LOG(ERROR)("Not support parameter: %s", name.c_str());
    }
}


void AlignmentGraph::SetParameter(const std::string &name, double v) {
    LOG(ERROR)("Not support parameter: %s", name.c_str());
}

void AlignmentGraph::ParseScoreParamter(const std::string &opts) {
    // weight:lc=30
    // count:

    auto ss = SplitStringByChar(opts, ':');
    if (ss.size() >= 1) {
        if (ss[0] == "count") {
            LinkScore = &AlignmentGraph::LinkScoreCount;
        } else if (ss[0] == "weight") {
            LinkScore = &AlignmentGraph::LinkScoreWeight;
        } else {
            LOG(ERROR)("Not support parameter: score=%s", opts.c_str());
        }

        for (size_t i=1; i<ss.size(); ++i) {
            auto kv = SplitStringByChar(ss[i], '=');
            if (kv[0] == "lc") {
                opts_.range[0] = std::stoi(kv[1]);
            } else if (kv[0] == "cov") {
                auto sss = SplitStringByChar(kv[1], ',');
                if (sss.size() >= 1 && sss[0].size() > 0) {
                    opts_.range[0] = std::stoi(sss[0]);
                }
                if (sss.size() >= 2 && sss[1].size() > 0) {
                    opts_.range[1] = std::stoi(sss[1]);
                }
                if (sss.size() >= 3 && sss[2].size() > 0) {
                    opts_.rate[0] = std::stod(sss[2]);
                }
                if (sss.size() >= 4 && sss[3].size() > 0) {
                    opts_.rate[1] = std::stod(sss[3]);
                }
            } else if (kv[0] == "rd") {
                opts_.reduction_ = std::stod(kv[1]);
            } else if (kv[0] == "bs") {
                auto sss = SplitStringByChar(kv[1], ',');
                if (sss.size() >= 1 && sss[0].size() > 0) {
                    opts_.branch_score_[0] = std::stod(sss[0]);
                }
                if (sss.size() >= 2 && sss[1].size() > 0) {
                    opts_.branch_score_[1] = std::stod(sss[1]);
                }
                if (sss.size() >= 3 && sss[2].size() > 0) {
                    opts_.branch_score_[2] = std::stod(sss[2]);
                }
            } else if (kv[0] == "msel") {
                opts_.min_selected = std::stoi(kv[1]);
            } else if (kv[0] == "wr") {
                auto ws = SplitStringByChar(kv[1], ',');
                if (ws.size() != 2) LOG(ERROR)("paramter 'weight' format is 0.2,0.8");
                opts_.weight_range_[0] = std::stod(ws[0]);
                opts_.weight_range_[1] = std::stod(ws[1]);
            } else {
                LOG(ERROR)("Not support parameter: score=...%s...", ss[i].c_str());
            }
        }
    } else {
        LOG(ERROR)("Not support parameter: score=%s", opts.c_str());
    }
}

void AlignmentGraph::Build(const DnaSeq& target, const std::array<size_t,2> &range, const std::vector<Alignment> &aligned) {
    Clear();

    DEBUG_printf("target range: %zd %zd %zd\n", range[0], range[1], target.Size());
    // add Tags
    AddTarget(target, range);
    tid_ = aligned[0].tid;

    for (size_t i=0; i<aligned.size(); ++i) {
        AddQuery(i, aligned[i].query_start, aligned[i].aligned_query, aligned[i].target_start, aligned[i].aligned_target);

        query_infos_.scores_.push_back(Score(aligned[i].qid, aligned[i].target_start, aligned[i].target_end, cols.size()));
    }

    std::sort(tags_.begin(), tags_.end());

    size_t ics = 0;
    for (size_t ice=0; ice<tags_.size(); ++ice) {
        if (tags_[ics].curr.col == tags_[ice].curr.col) continue;
        BuildCol(ics, ice);
        ics = ice;
    }
    BuildCol(ics, tags_.size());

    tags_.clear();
    if (print_rubbish) SaveGraph("sss.txt", 0, cols.size());
}


void AlignmentGraph::BuildCol(size_t stag, size_t etag) {
    assert(stag < etag && etag <= tags_.size());
    size_t rowsize = tags_[etag-1].curr.row + 1;
    size_t icol = tags_[stag].curr.col;

    cols[icol].rows.assign(rowsize, NodeGroup());
    cols[icol].coverage = 0;
    cols[icol].weight = 0;
    cols[icol].queries.reset();

    size_t s = stag;
    for (size_t i = s; i < etag; ++i) {
        if (tags_[i].curr.row == 0) {
            cols[icol].queries.set(tags_[i].id, true);
            cols[icol].coverage ++;
        }
        
        if (tags_[i].curr == tags_[s].curr) continue;

        BuildNode(s, i);
        s = i;
    } 
    BuildNode(s, etag);

}

void AlignmentGraph::BuildNode(size_t stag, size_t etag) {
    assert(etag > stag);

    Node& node = (*this)[tags_[stag].curr];

    for (size_t i = stag; i < etag; ++i) {
        if (node.links.size() > 0 && tags_[i].prev == node.links.back().prev) {
            node.links.back().seqs.set(tags_[i].id, true);
            node.links.back().count++;
            node.seqs.set(tags_[i].id, true);
        } else {
            node.links.push_back(Link());
            node.links.back().prev = tags_[i].prev;
            node.links.back().seqs.set(tags_[i].id, true);
            node.links.back().count = 1;
            node.seqs.set(tags_[i].id, true);
        }
        node.count++;
    }
    node.links.shrink_to_fit();
}   

void AlignmentGraph::Clear() {
    cols.clear();
    target_ = nullptr;
    query_infos_.Clear();
    sequence_.clear();
    quality_.clear();
    tags_.clear();
}


AlignmentGraph::Segment AlignmentGraph::FindBestPathBasedOnCount() {
    Segment seg;

    double global_score = -1;

    for (size_t i = 0; i < cols.size(); i++) {
        
        for (size_t j = 0; j < cols[i].Size(); j++) {
            for (size_t k = 0; k < cols[i][j].Size(); k++) {

                Node &node = cols[i][j][k];
                
                node.SortLinks();

                for (const auto &link : node.links) {
                    if (link.count == 0) continue;
                    const Loc& prev = link.prev;
                    
                    double score = LinkScoreCount(i, j, link) +  (prev.col == -1 ? 0 : cols[prev.col][prev.row][prev.base].score);

                    if (score > node.score) {
                        node.score = score;
                        node.best_link = &link;
                    }
                }

                if (node.score > global_score) {
                    global_score = node.score;
                    seg.end = Loc(i, j, k);                       
                }
            }
        }
    }

    return seg;

}

AlignmentGraph::Segment AlignmentGraph::FindBestPathBasedOnWeight() {
    Segment seg ;

    double global_score = -1;
    
    ComputeSimilarity4();

    score_range_ = { 1.0, -1.0};
    for (size_t i=0; i<query_infos_.scores_.size(); ++i) {
        if (query_infos_.selected_.count(i) > 0) {
            auto s = query_infos_.scores_[i].WeightInGraph();
            score_range_[0] = std::min<double>(s, score_range_[0]);
            score_range_[1] = std::max<double>(s, score_range_[1]);
        }
    }
    // pre-compute
    for (size_t col = 0; col < cols.size(); col++) {
        cols[col].weight = 0;
        cols[col].selected = 0;
        if (cols[col].queries[0]) cols[col].weight += 0.5;  // TODO Target score

       // if (cols[col].queries[0]) cols[col].weight += 1;
        for (size_t i=0; i<query_infos_.scores_.size(); ++i) {
            if (cols[col].queries[i+1] && query_infos_.selected_.count(i) > 0) {
                cols[col].weight += query_infos_.scores_[i].WeightInGraph(score_range_, opts_.weight_range_);
                cols[col].selected += 1;
            }
        }
    }

    for (size_t i = 0; i < cols.size(); i++) {
        if (cols[i].selected < sopts_.min_coverage) continue;
        for (size_t j = 0; j < cols[i].Size(); j++) {
            for (size_t k = 0; k < cols[i][j].Size(); k++) {

                Node &node = cols[i][j][k];
                
                for (const auto &link : node.links) {
                    if (link.count == 0) continue;

                    Loc prev = link.prev;

                    double score = LinkScoreWeight(i, j, link) + (prev.col == -1 ? 0 : cols[prev.col][prev.row][prev.base].score);
                    DEBUG_printf("FFF (%zd, %zd, %zd) <- (%zd, %zd,%zd),%f\n", prev.col, prev.row, prev.base, i, j, k, score);

                    if (score > node.score) {
                        node.score = score;
                        node.best_link = &link;
                    }
                }

                if (node.score > global_score) {
                    global_score = node.score;
                    seg.end = Loc(i, j, k);                     
                }
            }
        }
    }

    return seg;

}


void AlignmentGraph::Reconstruct(const std::vector<Segment>& segs) {

    std::vector<std::string> cns;
    for (const auto& seg : segs) {
        std::string s= ReconstructSimple(seg);
        if (seg.type == 0 || seg.end.col - seg.begin.col > opts_.max_bubble_length_ ) {
            cns.push_back(ReconstructSimple(seg));
        } else {
            cns.push_back(ReconstructComplex(seg));
        }
    }

    sequence_ = cns[0];
    for (size_t i = 1; i < segs.size(); ++i) {
        size_t off = segs[i].begin.base == 4 ? 0 : 1;           // if the base is not '-', skip the base
        sequence_.insert(sequence_.end(), cns[i].begin()+off, cns[i].end());
    }
}

std::string AlignmentGraph::ReconstructSimple(const Segment& seg) {
    std::string cns;
    std::string qual;
    const Node *curr_node = Get(seg.end);
    Loc loc = seg.end;

    const std::vector<std::string> toBase = {"A", "C", "G", "T", ""};

    auto valid = [this](const Node *n, const Loc& l) {
        return cols[l.col].coverage >= (size_t)sopts_.min_coverage || n->best_link->count >= (size_t)sopts_.min_coverage / 2; 
    };
    // 找到有效区域
    std::vector<std::array<Loc,2>> range;
    int state = 0;  //
    Loc start = seg.end;
    Loc end;
    int bad_count = 0;
    
    while (loc.col >= 0 ) {
        if (curr_node->best_link != nullptr && loc != seg.begin) {
            DEBUG_printf("col: (%zd,%zd,%zd) %zd  %zd  %zd\n", loc.col, loc.row, loc.base, 
                cols[loc.col].coverage, sopts_.min_coverage, curr_node->best_link->count);
            if (state == 0) {
                if (valid(curr_node, loc)) {
                    start = loc;
                    state = 1;
                }
            } else if (state == 1) {
                if (!valid(curr_node, loc)) {
                    end = loc;
                    state = 2;
                }
            } else if (state == 2) {
                if (!valid(curr_node, loc)) {
                    bad_count ++;
                    if (bad_count > 4) {
                        range.push_back({start, end});
                    }
                } else {
                    state = 1;
                    bad_count = 0;
                }
            } else {
                assert(0 && "never come here");
            }
            loc = curr_node->best_link->prev;
            curr_node = Get(loc);
        } else {
            break;
        }
    }
    if (state == 1) {
        range.push_back({start, loc});
    }

    auto mx = std::max_element(range.begin(), range.end(), [](const std::array<Loc,2>& r0, const std::array<Loc,2> &r1) {
        return r0[0].col - r0[1].col < r1[0].col - r1[1].col;
    });

    if (range.size() > 0) {
        DEBUG_printf("vvv %zd %d %d %d %d\n", range.size(), (*mx)[0].col, (*mx)[1].col, seg.begin.col, seg.end.col);
        loc = (*mx)[0];
        true_range_[1] = loc.col;
        if (loc.col >=0 ) curr_node = Get(loc);
        while (loc.col >= 0 ) {
            cns += toBase[loc.base];
            if (curr_node->best_link != nullptr && loc != (*mx)[1]) {
                if (!valid(curr_node, loc) ) {
                    break;
                }
                loc = curr_node->best_link->prev;
                curr_node = Get(loc);
            } else {
                break;
            }
            if (loc.col != -1) true_range_[0] = loc.col; 
        }
    }
    DEBUG_printf("true_range: %zd %zd\n", true_range_[0], true_range_[1]);
    std::reverse(cns.begin(), cns.end());
    return cns;
}

std::string AlignmentGraph::ReconstructComplex(const Segment& seg) {
   
    assert(!"TODO");

    return "";
}


std::vector<AlignmentGraph::Segment> AlignmentGraph::SplitSegment(const Loc &end) {

    std::vector<Segment> segments;

    auto is_clear_node = [](const Node* n, size_t cov) {
        if (cov <= 10) return true;         // coverage偏少，则分为简单区域

        if (n->best_link == &n->links[0] && n->best_link->count >= cov / 2) {
            if (n->links.size() == 1 || (n->links.size() > 1 && n->links[0].count / 2 >= n->links[1].count)) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    };
    
    std::vector<Loc> traceback;
    std::vector<uint8_t> state; // 0 clear, 1 unclear


    const Node *curr_node = Get(end);
    Loc loc = end;
    while (loc.col >= 0) {
        traceback.push_back(loc);
        DEBUG_printf("seg is_complex: (%d %d %d), %d, %d\n", loc.col, loc.row, loc.base, cols[loc.col].coverage, is_clear_node(curr_node, cols[loc.col].coverage));
        state.push_back(is_clear_node(curr_node, cols[loc.col].coverage)? 0 : 1);
 
        if (curr_node->best_link != nullptr) {
            loc = curr_node->best_link->prev;
            curr_node = Get(loc);
        } else {
            break;
        }
    }


    const int SOLID_LEN = 10;

    Segment seg;
    seg.end = traceback[0];
    seg.begin = traceback.back();
    seg.type = state[0];
    size_t startIndex = 0;

    size_t index = 0;
    while (index < traceback.size()) {

        if (seg.type == 0) {
            for (; index < traceback.size(); ++index) {
                if (state[index] == 0) continue;
                else break;
            }

            if (index < traceback.size()) {
                assert(state[index] == 1);
                
                // 寻找分割点
                int count = 0;
                size_t endIndex = index;
                for (; endIndex> startIndex; --endIndex) {
                    if (traceback[endIndex-1].row == 0) count ++;
                    if (count >= SOLID_LEN/2) break;
                }
                if (count >= SOLID_LEN/2) {
                    seg.begin = traceback[endIndex-1];
                    segments.push_back(seg);
                    seg.end = traceback[endIndex-1];
                    seg.begin = traceback.back();
                    seg.type = 1;
                } else {
                    seg.type = 1;
                }

            } else {
                segments.push_back(seg);
            }
        } else {
            assert(seg.type == 1);
            int count = 0;
            for (; index < traceback.size(); ++index) {
                if (state[index] == 0) {
                    if (traceback[index].row == 0) {
                        count ++;
                        if (count == SOLID_LEN) break;
                    }
                } else {
                    count = 0;
                }
            }

            if (index < traceback.size()) {
                int count = 0;
                size_t endIndex = index;
                for (; endIndex > startIndex; --endIndex) {
                    if (traceback[endIndex-1].row == 0) count ++;
                    if (count >= SOLID_LEN/2) break;
                }
                assert(count >= SOLID_LEN/2);
                seg.begin = traceback[endIndex-1];
                segments.push_back(seg);
                seg.end = traceback[endIndex-1];
                seg.begin = traceback.back();
                seg.type = 0;
            
            } else {
                segments.push_back(seg);
            }
        }
    }

    std::reverse(segments.begin(), segments.end());

    return segments;

}


std::vector<AlignmentGraph::Segment> AlignmentGraph::SplitSegment2(const Loc &end) {

    std::vector<Segment> segments;

    std::vector<Loc> traceback;
    std::vector<uint8_t> state; // 0 clear, 1 unclear


    std::vector<bool> simple_cols(cols.size());
    for (size_t i=0; i<cols.size(); ++i) {
        simple_cols[i] = IsSimpleColumn(i);
    }

    const Node *curr_node = Get(end);
    Loc loc = end;
    while (loc.col >= 0) {
        traceback.push_back(loc);
        
        state.push_back(simple_cols[loc.col] ? 0 : 1);
 
        if (curr_node->best_link != nullptr) {
            loc = curr_node->best_link->prev;
            curr_node = Get(loc);
        } else {
            break;
        }
    }


    const int SOLID_LEN = 12;

    Segment seg;
    seg.end = traceback[0];
    seg.begin = traceback.back();
    seg.type = state[0];
    size_t startIndex = 0;

    size_t index = 0;
    while (index < traceback.size()) {

        if (seg.type == 0) {
            for (; index < traceback.size(); ++index) {
                if (state[index] == 0) continue;
                else break;
            }

            if (index < traceback.size()) {
                assert(state[index] == 1);
                
                // 寻找分割点
                int count = 0;
                size_t endIndex = index;
                for (; endIndex> startIndex; --endIndex) {
                    if (traceback[endIndex-1].row == 0) count ++;
                    if (count >= SOLID_LEN/2) break;
                }
                if (count >= SOLID_LEN/2) {
                    seg.begin = traceback[endIndex-1];
                    segments.push_back(seg);
                    seg.end = traceback[endIndex-1];
                    seg.begin = traceback.back();
                    seg.type = 1;
                } else {
                    seg.type = 1;
                }

            } else {
                segments.push_back(seg);
            }
        } else {
            assert(seg.type == 1);
            int count = 0;
            for (; index < traceback.size(); ++index) {
                if (state[index] == 0) {
                    if (traceback[index].row == 0) {
                        count ++;
                        if (count == SOLID_LEN) break;
                    }
                } else {
                    count = 0;
                }
            }

            if (index < traceback.size()) {
                int count = 0;
                size_t endIndex = index;
                for (; endIndex > startIndex; --endIndex) {
                    if (traceback[endIndex-1].row == 0) count ++;
                    if (count >= SOLID_LEN/2) break;
                }
                assert(count >= SOLID_LEN/2);
                seg.begin = traceback[endIndex-1];
                segments.push_back(seg);
                seg.end = traceback[endIndex-1];
                seg.begin = traceback.back();
                seg.type = 0;
            
            } else {
                segments.push_back(seg);
            }
        }
    }

    std::reverse(segments.begin(), segments.end());

    return segments;

}


AlignmentGraph::Loc AlignmentGraph::Locate(const Node& node) {
    for (size_t i = 0; i < cols.size(); ++i) {
        for (size_t j = 0; j < cols[i].Size(); ++j) {
            for (size_t k = 0; k < cols[i][j].Size(); ++k) {
                if (&cols[i][j][k] == &node)
                    return { (int)i, (int)j, (int)k };
            }
        }
    }
    return { -1, 0, 0 };
}

void AlignmentGraph::Consensus() {
    Segment seg = LinkScore == &AlignmentGraph::LinkScoreWeight 
        ? FindBestPathBasedOnWeight()
        : FindBestPathBasedOnCount();


    //auto segs = SplitSegment(seg.end);
    // DEBUG_printf("seg: size=%zd\n", segs.size());
    // for (auto &s : segs) {
    //     DEBUG_printf("seg: (%zd, %zd) %d\n", s.begin.col, s.end.col, s.type);
    // }
    if (seg.end.col > 0) {  // TODO should be replaced by assert(seg.end.col > 0 && "Must find one path");
        sequence_ = ReconstructSimple(seg);
    } else {
        sequence_ = "";
    }

}

void AlignmentGraph::AddTarget(const DnaSeq& target, const std::array<size_t, 2> &range) {
    assert(cols.size() == 0);   // It should be empty graph

    target_ = &target;
    range_ = range;

    cols.assign(target.Size(), Column());

    Loc prev_loc = { -1, 0, -1 };
    
    for (size_t i = range[0]; i < range[1]; ++i) {
        
        Loc curr_loc = { (int)i, 0, target[i]};
        tags_.push_back({curr_loc, prev_loc, 0});
        prev_loc = curr_loc;
    }    
    
}

void AlignmentGraph::AddQuery(size_t sid, size_t query_start, const std::string &aligned_query, size_t target_start, const std::string &aligned_target) {
    
    size_t row = 0;
    int index_t = (int)target_start - 1;   // move the index to -1 for the first operation ++.
    Loc prev_loc = Loc::Invalid();

    assert(aligned_target.size() == aligned_query.size());
    for (size_t i=0; i<aligned_target.size(); ++i) {

        char bt = aligned_target[i];
        char bq = aligned_query[i];
      
        if (bq != '-' && bt != '-') {           // match or mismatch
            ++index_t;
            row = 0;
        } else if (bq == '-' && bt != '-') {    // deletion
            ++index_t;
            row = 0;
        } else if (bq != '-' && bt == '-') {    // insert
            ++row;
        } else {
            assert(bq == '-' && bt == '-');     // trivial
            continue;
        }
        
        assert(bt == '-' || (*target_)[index_t] == Base2Num[bt]);

        if (index_t < (int)range_[0]) continue;
        if (index_t >= (int)range_[1]) break;

        Loc curr_loc(index_t, (int)(row), Base2Num[bq]);
        tags_.push_back({curr_loc, prev_loc, (int)sid+1});       // 0 for target
        prev_loc = curr_loc;
    }

}

int AlignmentGraph::Options::BranchThreshold(int cov) {
    if (cov > range[1]) {
        return rate[1] * cov;
    } else if (cov >= range[0]) {
        double slope = (rate[1]*range[1] - rate[0]*range[0]) / (range[1] - range[0]);
        return (cov - range[0]) * slope + rate[0]*range[0];
    } else {
        return INT_MAX;         // 很大的数目，同一比较
    }
}

std::vector<const AlignmentGraph::Link*> AlignmentGraph::CollectLinks(size_t i) {
    assert(i < cols.size());
    std::vector<const Link*> links;
    if (cols[i].Size() >= 1) {
        for (size_t j = 0; j < cols[i][0].Size(); ++j) {
            for (const auto& l : cols[i][0][j].links) {
                if (l.prev.col != -1) {
                    links.push_back(&l);
                }
            }
        }
    }
    return links;
}

void AlignmentGraph::ComputeSimilarity4() {
    assert(cols.size() > 0 );

    auto cands = CollectImportantBranches();
    DEBUG_printf("branches: %zd\n", cands.size());

    VerifyImportantBranches1(cands);

    auto find_link = [](const std::array<ImportantBranch::LinkCol,2>& links, size_t id) -> const Link* {
        for (size_t i=0; i<links.size(); ++i) {
            if (links[i].l->seqs[id]) { return links[i].l; }
        }
        return nullptr;
    };
    
    query_infos_.SplitWindows(range_);
    
    for (size_t il = 0; il < cands.size(); ++il) {
        if (!cands[il].valid) continue;

        size_t icol = cands[il].c;  
        DEBUG_printf("position: %zd\n", icol);

        auto& links = cands[il].links;
        auto tlink = find_link(cands[il].links, 0);
        
        for (size_t si = 0; si<query_infos_.scores_.size(); ++si) {
            int qid = si;
            auto &ss = query_infos_.scores_[si];

            ss.all_cross += 1;
            for (size_t iw = 0; iw < query_infos_.windows.size(); ++iw) {
                if (icol >= query_infos_.windows[iw][0] && icol <= query_infos_.windows[iw][1]) {
                    auto& bs = ss.block_scores[iw];
                    bs.all_cross += 1;

                }
            }
            if (cols[icol].queries[qid+1]) {
                query_infos_.scores_[si].cross += 1;

                auto qlink = find_link(links, qid+1);
                
                query_infos_.scores_[si].branches[icol] = qlink;
                if (tlink != nullptr) {
                    query_infos_.scores_[si].t_in_cross += 1;
                }
                if (qlink != nullptr) {
                    query_infos_.scores_[si].q_in_cross += 1;
                }

                if (tlink  != nullptr && qlink != nullptr) {
                    if (tlink == qlink) query_infos_.scores_[si].q_t_one += 1;
                    else                query_infos_.scores_[si].q_t_two += 1;
          
                }
                for (size_t iw = 0; iw < query_infos_.windows.size(); ++iw) {
                    if (icol >= query_infos_.windows[iw][0] && icol <= query_infos_.windows[iw][1]) {
                        auto& bs = ss.block_scores[iw];
                        bs.cross += 1;
                        if (tlink  != nullptr && qlink != nullptr) {
                            if (tlink == qlink) bs.q_t_one += 1;
                            else                bs.q_t_two += 1;
                        }

                    }
                }
            }
        }
    
    }
    query_infos_.SelectReads3(opts_.min_selected, range_);
}
   
// structure for analyzing important location and branches
auto AlignmentGraph::CollectImportantBranches() -> std::vector<ImportantBranch> {

    std::vector<ImportantBranch> cands;
    for (size_t i=range_[0]; i<range_[1]; ++i) {     // Skip the first and last bases
        std::vector<const Link*> links = CollectLinks(i);
        size_t link_count = std::accumulate(links.begin(), links.end(), 0, [](size_t a, const Link* b) {
            return a + b->count;
        });
        
        int min_coverage = opts_.range[0];
        if (links.size() >= 2 && (int)link_count >= min_coverage) {
            std::sort(links.begin(), links.end(), [](const Link* a, const Link *b) {
                return a->count > b->count || (a->count == b->count && a->prev < b->prev);
            });

            int branch_threshold = opts_.BranchThreshold(cols[i].coverage);
            DEBUG_printf("br_check: %zd, th=%d, cov=%d, count=%d  %d,%d\n", i, branch_threshold, cols[i].coverage, link_count,  links[0]->count, links[1]->count);
            if (links[1]->count >= (size_t)branch_threshold) {
                cands.push_back(ImportantBranch());
                
                cands.back().c = i;
                cands.back().links[0].l = links[0];
                cands.back().links[1].l = links[1];
                for (auto &lnk : cands.back().links) {
                    for (size_t r = 0; r < cols[i][0].Size(); ++r) {
                        for (const auto& l : cols[i][0][r].links) {
                            if (lnk.l == &l) {
                                lnk.r = r;
                            }
                        }
                    }
                }

            }
        }
    }
    return cands;
}


// verify location, locactions. They should look like this
//      ->->        |
//     /    \       |
//   ->      ->     |
//     \    /       |
//      ->->        |
void AlignmentGraph::VerifyImportantBranches1(std::vector<ImportantBranch>& cands) {
    // Merge adjacent branches
    std::vector<std::array<size_t, 2>> segs;
    for (size_t start = 0, end = 1; start < cands.size(); start = end, end = start+1) {

        for (; end < cands.size(); ++end) {
            if (cands[end-1].links[0].l->prev.col + 2 < cands[end].links[0].l->prev.col) {
                break;
            }
        }
        segs.push_back({start, end});
    } 
    DEBUG_printf("br_merged_size: %zd\n", segs.size());

    std::unordered_set<size_t> verified;
    for (const auto& sg : segs) {
        size_t start = sg[0];
        size_t end = sg[1];
        bool rz = VerifiyImportantBranch(cands, start, end);
        DEBUG_printf("br_check_merged %zd-%zd(%zd) result=%d\n", cands[start].c, cands[end-1].c, end-start, rz);
        if (rz) {
            assert(end > start);
            verified.insert(start);
            verified.insert(end-1);
            verified.insert((start+end-1) / 2);
        }
    }

    for (size_t i = 0; i < cands.size(); ++i) {
        if (verified.find(i) == verified.end()) {
            cands[i].valid = false;
        }
    }

    VerifyConsistent(cands);

}

void DEBUG_PrintBases(const std::string& msg, const std::vector<int> &bs) {
    DEBUG_printf("%s(%zd): ", msg.c_str(), bs.size());
    for (auto b : bs) {
        DEBUG_printf("%c", "ACGT-"[b]);
    }
    DEBUG_printf("\n");
}

bool HasHomopolymer(size_t n, const std::vector<int>& bs, size_t start, size_t end) {
    size_t count = 1;
    for (size_t i = start+1; i < end; ++i) {
        if (bs[i] == bs[i-1]) {
            count++;
        } else {
            count = 1;
        }
        if (count >= n) break;
    }
    return count >= n;
}


bool HasHomopolymer(size_t n, const std::vector<int>& head, const std::vector<int>& body, const std::vector<int> &tail) {
    std::vector<int> s(head.begin(), head.end());
    s.insert(s.end(), body.begin(), body.end());
    s.insert(s.end(), tail.begin(), tail.end());
    return HasHomopolymer(n, s, 0, s.size());
}

std::vector<int> AlignmentGraph::GetBranch(const Link* link, const Loc& right, size_t left, Loc& end) {
    std::vector<int> branch;
    if (right.base < 4) branch.push_back(right.base);
    Loc head = link->prev;
    const auto& stream = link->seqs;
    while (head.col > (int)left || (head.col == (int)left && head.row > 0)) {
        DEBUG_printf("pass extend running: %d, %d, %d\n", head.col, head.row, head.base);
        auto n = Get(head);
        if (n == nullptr) break;
        if (head.base < 4) branch.push_back(head.base);


        size_t th = opts_.BranchThreshold(cols[head.col].coverage);
        Link* mx = nullptr;
        for (auto &l : n->links) {
            //if (l.count >= th) {
                DEBUG_printf("pass extend count: %zd, %zd %zd, %f\n", l.count, th, (l.seqs & stream).count(), th*0.8);
                if ((l.seqs & stream).count() >= th*0.6) {
                    mx = &l;
                    break;
                }
            //}
        }
        if (mx == nullptr) break;
        head = mx->prev;  
    }

    if (head.col == (int)left && head.row == 0) {
        if (head.base < 4) branch.push_back(head.base);
        end = head;
    } else {
        DEBUG_printf("pass extend failed: %d == %zd, %d == 0, %d < 4\n", head.col, left, head.row, head.base);
        branch.clear();
    }
    std::reverse(branch.begin(), branch.end());
    return branch;
}


std::vector<int> AlignmentGraph::GetBranchLeft(const Loc& start, size_t len, const MyBitSet& seqs) {
    std::vector<int> ext;

    Loc curr = start;
    auto n = Get(curr);
    while (n != nullptr && ext.size() < len) {

        //size_t th = opts_.BranchThreshold(cols[curr.col].coverage);
        size_t th = seqs.count();
        Link* mx = nullptr;
        for (auto &l : n->links) {
            //if (l.count >= th) {
                if ((l.seqs & seqs).count() > th*0.5) {
                    mx = &l;
                    break;
                }
            //}
        }
        if (mx == nullptr) break;
        curr = mx->prev;
        if (curr.base < 4) ext.push_back(curr.base);
        n = Get(curr);  
    }

    std::reverse(ext.begin(), ext.end());
    return ext;
}

std::vector<int> AlignmentGraph::GetBranchRight(const Loc &start, size_t len, const MyBitSet& seqs) {
    std::vector<int> ext;

    Loc curr = start;
    size_t th = seqs.count();
    while (ext.size() < len) {
        DEBUG_printf("CURR: %d %d %d, %zd\n", curr.col, curr.row, curr.base, th);
        Link* mx = nullptr;
        for (size_t i = 0; i <= 4 && mx == nullptr; ++i) {
            Loc next = { curr.col+1, 0, (int)i};
            auto nn = Get(next);
            DEBUG_printf("next: %d %d %d, %lld\n", next.col, next.row, next.base, nn);
            if (nn != nullptr) {
                for (auto &l : nn->links) {
                    if (l.prev == curr) {
                        DEBUG_printf("==: %zd %zd %zd\n", l.seqs.count(), seqs.count(), (l.seqs & seqs).count());
                        if ((l.seqs & seqs).count() > th*0.5) {
                            mx = &l;
                            curr = next;
                            break;
                        }
                    }
                }
            }
        }

        for (size_t i = 0; i <= 4 && mx == nullptr; ++i) {
            Loc next = { curr.col, curr.row+1, (int)i};
            auto nn = Get(next);
            DEBUG_printf("next1: %d %d %d, %lld\n", next.col, next.row, next.base, nn);
            if (nn != nullptr) {
                for (auto &l : nn->links) {
                    assert(l.count == l.seqs.count());
                    if (l.prev == curr) {
                        DEBUG_printf("==1: %zd %zd %zd\n", l.seqs.count(), seqs.count(), (l.seqs & seqs).count());
                        if ((l.seqs & seqs).count() > th*0.5) {
                            mx = &l;
                            curr = next;
                            break;
                        }
                    }
                }
            }
        }
        if (mx == nullptr) break;
        
        if (curr.base < 4) ext.push_back(curr.base);
    }

    return ext;
}

bool AlignmentGraph::VerifiyImportantBranch(const std::vector<ImportantBranch>& cands, size_t start, size_t end) {

    const int N = 5;
    size_t tstart = cands[start].c;
    size_t tend =  cands[end-1].c + 1;

    // too close to ends
    if (tstart < range_[0] + N || tend >= range_[1] - N) {
        return false;
    }

    // get target bases
    std::vector<int> local(2*N + tend - tstart);
    for (size_t i = 0; i < 2*N + tend - tstart; ++i) {
        local[i] = (*target_)[tstart-N + i];
    }
    DEBUG_PrintBases("target", local);
    
    if (HasHomopolymer(5, local, N-3, N+tend-tstart+3)) {
        return false;
    }

    DEBUG_printf("pass homo\n");

    Loc left0;
    auto branch0 = GetBranch(cands[end-1].links[0].l, Loc(cands[end-1].c, 0, cands[end-1].links[0].r), cands[start].c-1, left0);
    DEBUG_PrintBases("branch_0", branch0);
    if (branch0.size() == 0) return false;

    Loc left1;
    auto branch1 = GetBranch(cands[end-1].links[1].l, Loc(cands[end-1].c, 0, cands[end-1].links[1].r), cands[start].c-1, left1);
    DEBUG_PrintBases("branch_1", branch1);
    if (branch1.size() == 0) return false;

    auto br_left_0 = GetBranchLeft(left0, 3, cands[start].links[0].l->seqs);
    auto br_left_1 = GetBranchLeft(left1, 3, cands[start].links[1].l->seqs);
    DEBUG_PrintBases("br_left_0", br_left_0);
    DEBUG_PrintBases("br_left_1", br_left_1);
    if (br_left_0.size() < 3 || br_left_1.size() < 3) return false;

    auto find_matched = [](const std::vector<int>& a, const std::vector<int>& b) -> std::array<size_t,2> {
        assert(a.size() == 3 && b.size() == 3);

        if (a[0] == b[0] && a[1] == b[1]) {
            return {0, 0};
        } else if (a[0] == b[1] && a[1] == b[2]) {
            return {0, 1};
        } else if (a[1] == b[0] && a[2] == b[1]) {
            return {1, 0};
        } else if (a[1] == b[1] && a[2] == b[2]) {
            return {1, 1};
        } else {
            return {a.size(), b.size()};
        }
    };
    auto ibrl = find_matched(br_left_0, br_left_1);
    if (ibrl[0] == br_left_0.size() || ibrl[1] == br_left_1.size() ) return false;
    
    auto br_right_0 = GetBranchRight({(int)cands[end-1].c, 0, cands[end-1].links[0].r}, 3, cands[end-1].links[0].l->seqs);
    DEBUG_PrintBases("br_right_0", br_right_0);
    auto br_right_1 = GetBranchRight({(int)cands[end-1].c, 0, cands[end-1].links[1].r}, 3, cands[end-1].links[1].l->seqs);
    DEBUG_PrintBases("br_right_1", br_right_1);
    if (br_right_0.size() < 3 || br_right_1.size() < 3) return false;

    auto ibrr = find_matched(br_right_0, br_right_1);
    if (ibrr[0] == br_right_0.size() || ibrr[1] == br_right_1.size() ) return false;
    
    DEBUG_printf("pass find_matched: %zd,%zd   %zd,%zd\n", ibrl[0], ibrl[1], ibrr[0], ibrr[1]);

    if (HasHomopolymer(5, br_left_0, branch0, br_right_0) || HasHomopolymer(5, br_left_1, branch1, br_right_1)) return false;

    std::vector<int> bubble0(br_left_0.begin()+ibrl[0], br_left_0.end());
    bubble0.insert(bubble0.end(), branch0.begin(), branch0.end());
    bubble0.insert(bubble0.end(), br_right_0.begin(), br_right_0.begin()+ibrr[0]+2);
    
    std::vector<int> bubble1(br_left_1.begin()+ibrl[1], br_left_1.end());
    bubble1.insert(bubble1.end(), branch1.begin(), branch1.end());
    bubble1.insert(bubble1.end(), br_right_1.begin(), br_right_1.begin()+ibrr[1]+2);
    DEBUG_PrintBases("bubble0", bubble0);
    DEBUG_PrintBases("bubble1", bubble1);
    
    auto has_insert = [](const std::vector<int> &branch0, const std::vector<int> &branch1) {
        size_t ihead = 0;
        for (; ihead <= std::min<size_t>(branch0.size(), branch1.size()); ++ihead) {
            if (branch0[ihead] != branch1[ihead]) break;
        }

        size_t itail = 0;
        for (; itail <= std::min<size_t>(branch0.size(), branch1.size()); ++itail) {
            if (branch0[branch0.size()-itail-1] != branch1[branch1.size()-itail-1]) break;
        }
        DEBUG_printf("head-tail %zd, %zd\n", ihead, itail);
        if (ihead + itail >= std::min<size_t>(branch0.size(), branch1.size())) return true;

        DEBUG_PrintBases("mid0", std::vector<int>(branch0.begin()+ihead, branch0.end()-itail));
        DEBUG_PrintBases("mid1", std::vector<int>(branch1.begin()+ihead, branch1.end()-itail));
        if (branch0.size() <= branch1.size()) {
            auto s = std::search(branch1.begin()+ihead, branch1.end()-itail, branch0.begin()+ihead, branch0.end()-itail);
            if (s != branch1.end()-itail) return true;
        } else {
            auto s = std::search(branch0.begin()+ihead, branch0.end()-itail, branch1.begin()+ihead, branch1.end()-itail);
            if (s != branch0.end()-itail) return true;
        }
        return false;
    };

    if (has_insert(bubble0, bubble1) && bubble0.size() <= bubble1.size() + 1 && bubble0.size() +1 >= bubble1.size())  return false;
    
    DEBUG_printf("pass has_insert\n");
    return true;
}


void AlignmentGraph::VerifyConsistent(std::vector<ImportantBranch>& brs) {

    const int MIN_COUNT = 6;
    const int MIN_INV = 100;
    const double MAX_CONSIST = 0.80;
    const double MIN_CONSIST = 0.75;
    const double MIN_CONSIST_COUNT = 2;
    const double MIN_CONSIST_RATE = 0.1;


    // calcuate consistent sorce between each branch pair.
    std::vector<std::array<uint16_t, 2>> scores(brs.size()*brs.size(), {0, 0});
    std::vector<std::array<uint16_t,2>>  v_count(brs.size(), {0, 0});

    for (size_t i = 0; i < brs.size(); ++i) {
        if (!brs[i].valid) continue;

        for (size_t j = 0; j < i; ++j) {
            if (!brs[j].valid) continue;

            std::array<uint16_t, 4> links = {0, 0, 0, 0};
            links[0] = (brs[i].links[0].l->seqs & brs[j].links[0].l->seqs).count();
            links[1] = (brs[i].links[0].l->seqs & brs[j].links[1].l->seqs).count();
            links[2] = (brs[i].links[1].l->seqs & brs[j].links[0].l->seqs).count();
            links[3] = (brs[i].links[1].l->seqs & brs[j].links[1].l->seqs).count();

            DEBUG_printf("br_calc_cons (%zd, %zd) %zd,%zd,%zd,%zd\n", brs[i].c, brs[j].c, links[0], links[1], links[2], links[3]);

            auto& ss = scores[i*brs.size() + j];
            ss[0] = links[0] + links[1] + links[2] + links[3];
            ss[1] = std::max(links[0] + links[3], links[1] + links[2]);

            if (brs[i].c < brs[j].c + MIN_INV && brs[i].c + MIN_INV > brs[j].c) continue;
            if (ss[0] < MIN_COUNT) continue;
            v_count[i][0] ++;
            v_count[j][0] ++;

            if (ss[1]*1.0 / ss[0] >= MAX_CONSIST) {
                v_count[i][1] ++;
                v_count[j][1] ++;
            }
        } 
    }


    for (size_t i = 0; i < brs.size(); ++i) {
        if (!brs[i].valid) continue;

        bool rz = v_count[i][1] < std::max<uint16_t>(MIN_CONSIST_COUNT, v_count[i][0]*MIN_CONSIST_RATE);
        DEBUG_printf("check_br_cons (%zd) %zd, %zd %d\n", brs[i].c, v_count[i][0], v_count[i][1], rz);
        if (rz) {
            brs[i].valid = false;
        } 
    }
}


void AlignmentGraph::QueryInfos::SelectReads3(size_t min_sel, const std::array<size_t,2> &range) {

    // select queries for correction
    size_t win_count = windows.size();
    assert(win_count >= 1);
    std::vector<std::vector<double>> segscores(windows.size());
    std::vector<double> score_all;

    auto is_enough_overlap = [](size_t r0, size_t r1, size_t w0, size_t w1) {
        size_t ol0 = std::max<size_t>(r0, w0);
        size_t ol1 = std::min<size_t>(r1, w1);

        return ol1 > ol0 && ((ol1 - ol0)*2 > (w1 - w0) || (ol1 - ol0)*2 > r1 - r0);
    };

    std::unordered_set<Seq::Id> added;
    for (size_t i = 0; i < scores_.size(); ++i) {
        const auto &s = scores_[i];
        score_all.push_back(s.Weight());
        DEBUG_printf("w-all: %zd %f %zd %zd\n", i, s.Weight(), s.tstart, s.tend);
    
        if (win_count > 1) {
            for (size_t j = 0; j < windows.size(); ++j) {
                DEBUG_printf("cross(%zd,%zd): %d / %d\n",i, j, s.block_scores[j].cross , s.block_scores[j].all_cross);
                if (is_enough_overlap(s.tstart, s.tend, windows[j][0], windows[j][1])) {
                    segscores[j].push_back(s.block_scores[j].Weight());
                    added.insert(s.qid);
                    DEBUG_printf("w-s%zd: %zd %f %d %d\n", j, i, segscores[j].back(), s.block_scores[j].q_t_one, s.block_scores[j].q_t_two);
                }
            }
        }
    }

    std::vector<double> segths(segscores.size());
    if (win_count > 1) {
        for (size_t i = 0; i < segths.size(); i++) {
            segths[i] = std::max(0.0, FindScoreThreshold3(segscores[i])-0.000001);
            DEBUG_printf("th(%zd): %f\n", i, segths[i]);
        }
    }
    double th_all = std::max(0.0, FindScoreThreshold3(score_all));
    DEBUG_printf("th(all): %f\n", th_all);
    
    std::unordered_set<int> excluded;
    std::unordered_set<int> all_selected;
    selected_.clear();
    for (size_t i = 0; i < scores_.size(); ++i) {
        const auto &s = scores_[i];
        if (s.Weight() >= th_all) {
            selected_.insert(i);
            all_selected.insert(i);
            DEBUG_printf("sel(all): %zd %f\n", i, s.Weight());
        }

        if (win_count > 1) {
            for (size_t j = 0; j < windows.size(); ++j) {
                if (is_enough_overlap(s.tstart, s.tend, windows[j][0], windows[j][1])) {
                    if (s.block_scores[j].Weight() >= segths[j]) {
                        selected_.insert(i);
                        DEBUG_printf("sel(%zd): %zd %f\n", j, i, s.block_scores[j].Weight());
                    } else if (s.block_scores[j].Weight() < 0 && s.all_cross >= 10) {
                        excluded.insert(i);
                    }
                }
            }
        }
        
    }
    DEBUG_printf("sel: %zd\n", selected_.size());
    for (auto e : excluded) {
        if (all_selected.find(e) == all_selected.end())
            selected_.erase(e);
    }
    DEBUG_printf("sel2: %zd\n", selected_.size());

    // add 
    for (size_t i = 0; i < scores_.size(); ++i) {
        const auto &s = scores_[i];

        if (s.q_t_two == 0) {
            selected_.insert(i);
        }
        
    }
    DEBUG_printf("sel: %zd\n", selected_.size());

    if (selected_.size() < min_sel) {
        std::vector<size_t> score_index(scores_.size());
        for (size_t i=0; i<score_index.size(); ++i) {
            score_index[i] = i;
        }
        
        std::sort(score_index.begin(), score_index.end(), [this](int a, int b) {
            return scores_[a].Weight() > scores_[b].Weight();
        });

        for (size_t i=0; i<score_index.size(); ++i) {
            selected_.insert(score_index[i]);
            if (selected_.size() >= min_sel) break;
        }
    }
}


double AlignmentGraph::QueryInfos::FindScoreThreshold3(const std::vector<double>& score) const {

    // parameters b w
    const double b = 0.01;
    const int w = 10;
    int n = int(2 / b) + 1;     // all scores are in range [-1,1]
    std::vector<int> hist(n, 0);

    for (auto s : score) {
        assert(s >= -1 && s <= 1);
        hist[int((s + 1) / b)] += 1;
    }

    std::vector<int> smoothed(hist.size()-w+1, 0);
    smoothed[0] = std::accumulate(hist.begin(), hist.begin()+w, 0);
    for (size_t i=1; i + w <hist.size(); ++i) {
        smoothed[i] = smoothed[i-1] - hist[i-1] + hist[i-1+w];
    }

    const int LOW = 4;
    const double DIFF = 0.25; 
    const int WIDTH = 40;
    int state = 0;
    int ps = 0;
    int accu = std::accumulate(hist.end()-w+1, hist.end(), 0);
    for (int i = smoothed.size() - 1 ; i > -1; --i) {
        //printf("line: %zd, %zd, %d, %d\n", i, smoothed[i], ps, smoothed[ps]);
        accu += hist[i];
        if (state == 0) {
            if (smoothed[i] >= LOW) { // SKIP smoothed[i] is small
                // found starting point
                ps = i;
                state = 1;
                //printf("start: %zd\n", i);
            }
        } else if (state == 1) {
            if (smoothed[ps] - smoothed[i] >= std::max<int>(smoothed[ps] * DIFF, LOW) && std::abs(i-ps) >= WIDTH / 2 ) {
                // Found peak
                //printf("peak: %zd, %zd\n", ps, i);
                state = 2;
                ps = i;
            } else {
                if (smoothed[i] > smoothed[ps]) {
                    ps = i;
                }
            }
        } else if (state == 2) {        // find troughs
            if (smoothed[i] - smoothed[ps] >= std::max<int>(smoothed[ps] * DIFF, LOW)) {
                //printf("troughs: %zd, %zd\n", ps, i);
                break;
            } else {
                if (smoothed[i] < smoothed[ps] || (smoothed[i] <= smoothed[ps] && i >= 100) ) {
                    ps = i;
                }
            }
        }
    }
    
    return (ps+w/2)  * b - 1;
}

size_t AlignmentGraph::QueryInfos::GetBlockSize() const {
    return std::accumulate(scores_.begin(), scores_.end(), 0, [](size_t a, const Score& b) {
        return a + b.tend - b.tstart;
    }) / scores_.size();
}

std::array<size_t, 3> AlignmentGraph::QueryInfos::GetWindowSize(const std::array<size_t, 2> &range) const {
    size_t step_size = std::accumulate(scores_.begin(), scores_.end(), 0, [](size_t a, const Score& b) {
        return a + b.tend - b.tstart;
    }) / scores_.size();

    step_size = std::max<size_t>(step_size/2, 10000);
    size_t count = (range[1] - range[0] + step_size / 2) / step_size;

    if (count <= 2) {
        size_t win_size = range[1] - range[0];
        return {win_size, (win_size + 1) / 2, 1};
    } else {
        step_size = (range[1] - range[0] + count - 1) / count;
        return {step_size*2, step_size, count - 1};
    }
}

void AlignmentGraph::QueryInfos::SplitWindows(const std::array<size_t, 2> &range) {
    auto bsizes = GetWindowSize(range);
    size_t win_size = bsizes[0];
    size_t step_size = bsizes[1];
    size_t win_count = bsizes[2];
    DEBUG_printf("win_size = %zd, step_size = %zd, win_count = %zd, range = (%zd, %zd)\n", 
        win_size, step_size, win_count, range[0], range[1]);

    assert(windows.size() == 0);
    size_t start = range[0];
    size_t end = start;
    while (end < range[1]) {
        //if (start + win_size + step_size / 2 > range[1]) {
        if (start + win_size > range[1]) {
            end = range[1];
        } else {
            end = start + win_size;
        }
        windows.push_back({start, end});
        start += step_size;
    }

    assert(windows.size() == win_count);
    
    for (auto &s : scores_) {
        s.block_scores.assign(win_count, Score::BlockScore());
    }
}

void AlignmentGraph::QueryInfos::SaveReadInfos(std::ostream& os, int tid, const ReadStore &rs) const {
    const std::string& tname = rs.QueryNameById(tid);
    for (size_t i=0; i<scores_.size(); ++i) {
        const auto & s = scores_[i];
        const std::string& qname = rs.QueryNameById(s.qid);
        bool sel = selected_.find(i) != selected_.end() ? 1 : 0;
        os << tname << " " << qname << " " << s.WeightInGraph() << " " << sel << "\n";
    }
}   

double mypow(double x, size_t n) { // GLIBC_2.29
    if (n > 0) {
        double y = x;
        for (int i = 1; i < n; ++i) {
            y *= x;
        }
        return y;
    } else {
        return 1;
    }
}

double AlignmentGraph::LinkScoreCount(size_t col, size_t row, const Link& link) { 
    
    double scale = std::max<double>(mypow(opts_.branch_score_[2], row)*opts_.branch_score_[0], opts_.branch_score_[1]);
    double compensate = std::max<double>(scale * cols[col].coverage, sopts_.min_coverage * opts_.branch_score_[0]);
    return link.count - compensate;
} 

double AlignmentGraph::LinkScoreWeight(size_t col, size_t row, const Link &link) {
    double s = 0;
    for (size_t i=0; i<query_infos_.scores_.size(); ++i) {
        if (link.seqs[i+1] && query_infos_.selected_.find(i) != query_infos_.selected_.end()) {
            s += query_infos_.scores_[i].WeightInGraph(score_range_, opts_.weight_range_);
        }
    }
    if (link.seqs[0]) {
        s += 0.5; // TODO Target score 
    }


    double scale = std::max<double>(mypow(opts_.branch_score_[2], row)*opts_.branch_score_[0], opts_.branch_score_[1]);
    double compensate = std::max<double>(scale * cols[col].weight, opts_.branch_score_[0] * cols[col].weight * sopts_.min_coverage / cols[col].coverage);
    
    return s - compensate;
}


bool AlignmentGraph::IsSimpleColumn(size_t col) {
    if (col + 1 < cols.size()) {
        auto links = CollectLinks(col+1);       // 下个col会收拢前一个的分支，方便统计

        std::sort(links.begin(), links.end(), [](const Link* a, const Link *b) {
            return a->count > b->count || (a->count == b->count && a->prev < b->prev);
        });

        return links.size() <= 1 || (links[0]->count * 3 >= cols[col].coverage && links[0]->count >= links[1]->count *4/2);
    } else {
        return true;
    }
}

void AlignmentGraph::SaveGraph(const std::string &fname, size_t s, size_t e) const {

    std::ofstream of(fname);

    of << "Target,Source,Weight,Seqs\n";
    for (size_t i=s; i<e; ++i) {
        const auto& c = cols[i];
        for (size_t ir = 0; ir < c.rows.size(); ++ir) {
            auto &r = c.rows[ir];

            for (size_t ib=0; ib< r.Size(); ++ib) {
                auto &b = r[ib];
            
                for (auto l: b.links) {
                    if (l.prev.col != -1) {
                        of << i << "_" << ir << "_" << "ACGT-"[ib] << ", "  
                           << l.prev.col << "_" << l.prev.row << "_" << "ACGT-"[l.prev.base] << "," 
                           << l.count << "," ;

                        if (l.seqs[0]) {
                            of << dataset_.QueryStringById(tid_) ;
                        }

                        for (size_t i = 0; i < query_infos_.scores_.size(); ++i) {
                            if (l.seqs[i+1]) {
                                of << '-' << dataset_.QueryStringById(query_infos_.scores_[i].qid);
                            }
                        }
                        of << '\n';
                    }
                }
            }
        }
    }

}


} // namespace fsa {
