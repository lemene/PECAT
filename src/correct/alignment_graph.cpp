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

namespace fsa {


DnaSerialTable2 AlignmentGraph::Base2Num;

AlignmentGraph::AlignmentGraph() {
}


void AlignmentGraph::SetParameter(const std::string &name, const std::string &opts) {
    if (name == "score") {
        ParseScoreParamter(opts);
    } else {
        LOG(ERROR)("Not support parameter: %s", name.c_str());
    }
}

void AlignmentGraph::SetParameter(const std::string &name, int v) {
    if (name == "min_coverage") {
        min_coverage_ = v;
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
            } else if (kv[0] == "rd") {
                reduction_ = std::stod(kv[1]);
            } else if (kv[0] == "bs") {
                auto sss = SplitStringByChar(kv[1], ',');
                if (sss.size() >= 1 && sss[0].size() > 0) {
                    branch_score_[0] = std::stod(sss[0]);
                }
                if (sss.size() >= 2 && sss[1].size() > 0) {
                    branch_score_[1] = std::stod(sss[1]);
                }
                if (sss.size() >= 3 && sss[2].size() > 0) {
                    branch_score_[2] = std::stod(sss[2]);
                }
            } else if (kv[0] == "msel") {
                min_selected = std::stoi(kv[1]);
            } else if (kv[0] == "wr") {
                auto ws = SplitStringByChar(kv[1], ',');
                if (ws.size() != 2) LOG(ERROR)("paramter 'weight' format is 0.2,0.8");
                weight_range_[0] = std::stod(ws[0]);
                weight_range_[1] = std::stod(ws[1]);
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

    // add Tags
    AddTarget(target, range);

    for (size_t i=0; i<aligned.size(); ++i) {
        AddQuery(i, aligned[i].query_start, aligned[i].aligned_query, aligned[i].target_start, aligned[i].aligned_target);

        query_infos_.scores_.push_back(Score(aligned[i].qid, cols.size()));
    }

    std::sort(tags_.begin(), tags_.end());

    size_t ics = 0;
    for (size_t ice=0; ice<tags_.size(); ++ice) {
        if (tags_[ics].curr.col == tags_[ice].curr.col) continue;
        BuildCol(ics, ice);
        ics = ice;
    }
    BuildCol(ics, tags_.size());

    //SaveGraph("sss.txt", 0, cols.size());
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
}   

void AlignmentGraph::Clear() {
    cols.clear();
    branches_.clear();
    branch_flags_.clear();
    target_ = nullptr;
    query_infos_.scores_.clear();
    sequence_.clear();
    quality_.clear();
    tags_.clear();
}


AlignmentGraph::Segment AlignmentGraph::FindBestPathBasedOnCount() {
    Segment seg;

    const Node* end = nullptr;
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
                    end = &node;
                    seg.end = Loc(i, j, k);                       
                }
            }
        }
    }

    return seg;

}

AlignmentGraph::Segment AlignmentGraph::FindBestPathBasedOnWeight() {
    Segment seg ;

    const Node* end = nullptr;
    double global_score = -1;

    ComputeSimilarity();

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
        if (cols[col].queries[0]) cols[col].weight += 0.5;

       // if (cols[col].queries[0]) cols[col].weight += 1;
        for (size_t i=0; i<query_infos_.scores_.size(); ++i) {
            if (cols[col].queries[i+1] && query_infos_.selected_.count(i) > 0) {
                cols[col].weight += query_infos_.scores_[i].WeightInGraph(score_range_, weight_range_);
                cols[col].selected += 1;
            }
        }
    }

    for (size_t i = 0; i < cols.size(); i++) {
        if (cols[i].selected < 4) continue;
        for (size_t j = 0; j < cols[i].Size(); j++) {
            for (size_t k = 0; k < cols[i][j].Size(); k++) {

                Node &node = cols[i][j][k];
                
                for (const auto &link : node.links) {
                    if (link.count == 0) continue;

                    Loc prev = link.prev;

                    double score = LinkScoreWeight(i, j, link) + (prev.col == -1 ? 0 : cols[prev.col][prev.row][prev.base].score);

                    if (score > node.score) {
                        node.score = score;
                        node.best_link = &link;
                    }
                }

                if (node.score > global_score) {
                    global_score = node.score;
                    end = &node;   
                    seg.end = {i, j, k};                     
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
        if (seg.type == 0 || seg.end.col - seg.begin.col > max_bubble_length_ ) {
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
    for (const auto &s : cns) {

    }
}

std::string AlignmentGraph::ReconstructSimple(const Segment& seg) {
    std::string cns;
    std::string qual;
    const Node *curr_node = Get(seg.end);
    Loc loc = seg.end;

    const std::vector<std::string> toBase = {"A", "C", "G", "T", ""};

    auto valid = [this](const Node *n, const Loc& l) {
        return cols[l.col].coverage >= min_coverage_ || n->best_link->count >= min_coverage_ / 2; 
    };
    // 找到有效区域
    std::vector<std::array<Loc,2>> range;
    int state = 0;  //
    Loc start = seg.end;
    Loc end;
    int bad_count = 0;
    
    while (loc.col >= 0 ) {
        if (curr_node->best_link != nullptr && loc != seg.begin) {
            DEBUG_printf("vvv %d cov=%d, count=%d, state=%d\n", loc.col, cols[loc.col].coverage, curr_node->best_link->count, state);
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
        }
    }
    std::reverse(cns.begin(), cns.end());
    return cns;
}

std::string AlignmentGraph::ReconstructComplex(const Segment& seg) {
   
   assert(!"TODO");

    return "";
}


std::vector<std::string> AlignmentGraph::RestoreSegment(const Segment &seg) {
    std::vector<std::string> seqs;

    assert(seg.begin.col != -1);

    std::vector<Loc> stack_loc;
    std::vector<size_t> stack_link;
    std::vector<std::bitset<MAX_COV>> stack_seqs;
    const std::vector<std::string> toBase = {"A", "C", "G", "T", ""};

    for (size_t i=0; i<5; ++i) {
        Loc loc = seg.end;
        loc.base = i;
        stack_loc.push_back(loc);  
        stack_link.push_back(0);
        stack_seqs.push_back(Get(loc)->seqs);

        while (stack_loc.size() > 0) {
            // 检查是否到达边界
            if (stack_loc.back().col == seg.begin.col && stack_loc.back().row == 0) {
                std::string s;
                for (auto loc : stack_loc) {
                    s += toBase[loc.base];
                }
                std::reverse(s.begin(), s.end());
                
                for (size_t i=0; i<stack_seqs.back().count(); ++i) seqs.push_back(s);
                //printf("find one:  %zd %s\n", stack_seqs.back().count(), s.c_str());
                stack_link.pop_back();
                stack_loc.pop_back();
                stack_seqs.pop_back();

            }

            auto n = Get(stack_loc.back());
            // 检查能否继续延长
            if (stack_link.back() < n->links.size()) {
                Link& next_link = n->links[stack_link.back()];
                auto next_seqs = (next_link.seqs & stack_seqs.back());
                if (next_seqs.count() > 0 && n->links[stack_link.back()].prev.col >= 0) {
                    stack_loc.push_back(n->links[stack_link.back()].prev);
                    stack_link.back()++;
                    stack_link.push_back(0);
                    stack_seqs.push_back(next_seqs);
                } else {
                    stack_link.back()++;
                }
                
            } else {
                stack_link.pop_back();
                stack_loc.pop_back();
                stack_seqs.pop_back();
            }
            
        }
    }

    return seqs;
}

std::vector<AlignmentGraph::Segment> AlignmentGraph::SplitSegment(const Loc &end) {

    std::vector<Segment> segments;

    auto is_clear_node = [](const Node* n, size_t cov) {
        if (cov <= 10) return true;         // coverage偏少，则分为简单区域

        if (n->best_link == &n->links[0] && n->best_link->count >= cov / 2) {
            if (n->links.size() == 1 || n->links.size() > 1 && n->links[0].count / 2 >= n->links[1].count) {
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
        //printf("is_complex: (%d %d %d), %d, %d\n", loc.col, loc.row, loc.base, cols[loc.col].coverage, is_clear_node(curr_node, cols[loc.col].coverage));
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
    for (auto seg : segments) {
        //printf("seg: %d %d %d, %s\n", seg.begin.col, seg.end.col, seg.type, ReconstructSimple(seg).c_str());
    }

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
    for (auto seg : segments) {
        //printf("seg: %d %d %d, %s\n", seg.begin.col, seg.end.col, seg.type, ReconstructSimple(seg).c_str());
    }

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
    branches_.assign(cols.size(), nullptr);
    branch_flags_.assign(cols.size(), false);

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

        if (index_t < range_[0]) continue;
        if (index_t >= range_[1]) break;

        Loc curr_loc(index_t, (int)(row), Base2Num[bq]);
        tags_.push_back({curr_loc, prev_loc, sid+1});       // 0 for target
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


std::vector<const AlignmentGraph::Link*> AlignmentGraph::CollectLinks1(size_t i) {
    std::vector<const Link*> links;


    assert (cols[i][0].Size() == 5);
    for (size_t j = 0; j < 4; ++j) {
        for (const auto& l : cols[i][0][j].links) {
            if (l.prev.col != -1) {
                links.push_back(&l);
            }
        }
    }

    for (size_t t=1; t < cols[i].Size(); ++t) {
        for (size_t j=0; j<cols[i][t].Size(); ++j) {
            for (const auto &l : cols[i][t][j].links) {
                if (l.prev == Loc {(int)i, 0, 4}) {
                    links.push_back(&l);
                }
            }
        }
    }

    for (const auto& l : cols[i+1][0][4].links) {
        if (l.prev == Loc {(int)i, 0, 4}) {
            links.push_back(&l);
        }
    }
    //for (size_t j = 0; j < cols[i+1][0].Size(); ++j) {
    //    for (const auto& l : cols[i+1][0][j].links) {
    //        if (l.prev == Loc {(int)i, 0, 4}) {
    //            links.push_back(&l);
    //        }
    //    }
    //}

    return links;
}


std::vector<const AlignmentGraph::Link*> AlignmentGraph::CollectLinks2(size_t i) {
    std::vector<const Link*> links;


    assert (cols[i][0].Size() == 5);
    for (size_t j = 0; j < 4; ++j) {
        for (const auto& l : cols[i][0][j].links) {
            if (l.prev.col != -1) {
                links.push_back(&l);
            }
        }
    }

    return links;
}

void AlignmentGraph::ComputeSimilarity() {
    assert(cols.size() > 0 );

    for (size_t i=range_[0]+1; i<range_[1]-1; ++i) {     // Skip the first and last bases
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

            size_t sz = 0;
            for (; sz < links.size(); ++sz) {
                if ((int)links[sz]->count < branch_threshold) break;    
            }

            if (sz >= 2) {
                auto find_link = [&links,sz](size_t id) -> const Link*{
                    for (size_t i=0; i<sz; ++i) {
                        if (links[i]->seqs[id]) { return links[i]; }
                    }
                    return nullptr;
                };

                const int qw = 1;

                auto tlink = find_link(0);
                branch_flags_[i] = true;
                branches_[i] = tlink;
                for (size_t si = 0; si<query_infos_.scores_.size(); ++si) {
                    int qid = si;
                    if (cols[i].queries[qid+1]) {
                        query_infos_.scores_[si].cross += qw;

                        auto qlink = find_link(qid+1);
                        
                        query_infos_.scores_[si].branches[i] = qlink;
                        if (tlink != nullptr) {
                            query_infos_.scores_[si].t_in_cross += qw;
                        }
                        if (qlink != nullptr) {
                            query_infos_.scores_[si].q_in_cross += qw;
                        }

                        if (tlink  != nullptr && qlink != nullptr) {
                            if (tlink == qlink) query_infos_.scores_[si].q_t_one+= qw;
                            else                query_infos_.scores_[si].q_t_two+= qw;
                        }
                    }
                }
            }
        }
    }

    // print point
    for (size_t i=range_[0]+1; i<range_[1]-1; ++i) {     // Skip the first and last bases
        if (branch_flags_[i]) {
            DEBUG_printf("branch point %zd\n", i);
        }
    }

    query_infos_.SelectReads(min_selected);
    
}

void AlignmentGraph::QueryInfos::SelectReads(int min_sel) {

    // select queries for correction
    std::vector<int> score_index(scores_.size());
    for (int i=0; i<score_index.size(); ++i) {
        score_index[i] = i;
    }
    
    std::sort(score_index.begin(), score_index.end(), [this](int a, int b) {
        return scores_[a].Weight() > scores_[b].Weight();
    });

    auto threshold = FindScoreThreshold2();
    DEBUG_printf("threshold %f\n", threshold);
    
    int count0 = 0;     
    int count1 = 0;
    int count2 = 0;
    for (auto s : scores_) {
        if (s.Weight() > threshold ) count0 ++;
        else if (s.Weight() == threshold ) count1 ++;
        else count2++;
    }

    selected_.clear();
    int ccc = count0;   
    size_t selsize = std::max(std::min<int>(ccc*1.0 + count1, (ccc+count1+count2)*3/2),   std::min(min_sel, ccc*2));
    selsize = std::max<size_t>(min_sel, selsize);

    for (size_t i=0; i<score_index.size(); ++i) {
        selected_.insert(score_index[i]);
        if (i >= selsize) break;
    }
}

double AlignmentGraph::QueryInfos::FindScoreThreshold() {

    // parameters b w
    const double b = 0.01;
    const int w = 10;

    int n = int(2 / b) + 1;
    std::vector<int> hist(n, 0);

    for (auto s : scores_) {
        hist[int((s.Weight() + 1) / b)] += 1;
    }

    std::vector<int> smoothed(hist.size()-w+1, 0);
    smoothed[0] = std::accumulate(hist.begin(), hist.begin()+w, 0);
    for (size_t i=1; i + w <hist.size(); ++i) {
        DEBUG_printf("hist: %zd %d\n", i, hist[i]);
        smoothed[i] = smoothed[i-1] - hist[i-1] + hist[i-1+w];
        DEBUG_printf("smooth: %zd %d\n", i, smoothed[i]);
    }

    std::vector<std::array<int, 3>> peaks;  // { left, right, peak}
    std::vector<int> troughs;

    int state = 0;
    int ps = 0;
    const double rate = 0.2;
    for (size_t i = 0 ; i < smoothed.size(); ++i) {
        if (state == 0) {
            if (smoothed[i] >= 6) {
                ps = i;
                state = 1;
            }
        } else if (state == 1) {
            if (smoothed[i] < std::max<int>(smoothed[ps] * 0.8, 6)) {
                int end = i - 1;
                int start = ps;
                for (int ii = ps; ii > 0; ii--) {
                    if (smoothed[ii] >= smoothed[ps] * 0.8 && smoothed[ii] >= 6) {
                        start = ii;
                    } else {
                        break;
                    }
                }
                peaks.push_back({start, end, ps});
                DEBUG_printf("pushpeak %d %d %d\n", start, end, ps);
                state = 2;
                ps = i;
            } else {
                if (smoothed[i] > smoothed[ps]) {
                    ps = i;
                }
            }
        } else if (state == 2) {        // find troughs
            DEBUG_printf("check low %d %d %d\n", i, smoothed[i] , smoothed[ps]);
            if (smoothed[i] * 0.8 > smoothed[ps] ) {
                state = 0;
            } else {
                if (smoothed[i] < smoothed[ps]) {
                    ps = i;
                }
            }
        }
    }

    DEBUG_printf("peak: %zd\n", peaks.size());
    

    if (peaks.size() >= 2) {
        const auto &last0 = peaks.back();
        const auto &last1 = peaks[peaks.size()-2];

        return (last1[1] + w - 1 + last0[0])/2 * b - 1;
        
        
    } else if (peaks.size() == 1) {
        const auto &last = peaks.back();
        return (last[0] + last[1] + w - 1)/2 * b - 1;
    } else {
        return 0.0;
    }
}


double AlignmentGraph::QueryInfos::FindScoreThreshold1() {

    // select queries for correction
    std::vector<size_t> score_index(scores_.size());
    for (size_t i=0; i<score_index.size(); ++i) {
        score_index[i] = i;
    }

    std::sort(score_index.begin(), score_index.end(), [this](int a, int b) {
        return scores_[a].Weight() > scores_[b].Weight();
    });

    // TODO two parameters: ratio and minimum number used for calculating MAD 
    size_t count = std::max<size_t>(std::min<size_t>(0, score_index.size()), score_index.size() / 4);
    std::vector<double> data(count, 0);
    for (size_t i=0; i < count; ++i) {
        data[i] = scores_[score_index[i]].Weight();
    }

    double median = 0;
    double mad = 0;
    ComputeMedianAbsoluteDeviation(data, median, mad);

    return median-3*1.4826*mad;
}


double AlignmentGraph::QueryInfos::FindScoreThreshold2() {

    // parameters b w
    const double b = 0.01;
    const int w = 10;

    int n = int(2 / b) + 1;
    std::vector<int> hist(n, 0);

    for (auto s : scores_) {
        hist[int((s.Weight() + 1) / b)] += 1;
    }

    std::vector<int> smoothed(hist.size()-w+1, 0);
    smoothed[0] = std::accumulate(hist.begin(), hist.begin()+w, 0);
    for (size_t i=1; i + w <hist.size(); ++i) {
        smoothed[i] = smoothed[i-1] - hist[i-1] + hist[i-1+w];
    }

    int state = 0;
    int ps = 0;
    int start = 0;
    int end = 0;
    int peak = 0;
    const double rate = 0.2;
    int accu = std::accumulate(hist.end()-w+1, hist.end(), 0);
    for (int i = smoothed.size() - 1 ; i > -1; --i) {
        accu += hist[i];
        if (state == 0) {
            if (smoothed[i] >= 6) {
                ps = i;
                state = 1;
                start = i;
            }
        } else if (state == 1) {
            if (smoothed[i] < std::max<int>(smoothed[ps] * 0.8, 6) && std::abs(i-ps) > 3) {
                DEBUG_printf("check peak %d %d %d\n", i, smoothed[i] , smoothed[ps]);
                
                peak = ps;
                state = 2;
                ps = i;
            } else {
                DEBUG_printf("check up %d %d %d\n", i, smoothed[i] , smoothed[ps]);
                if (smoothed[i] > smoothed[ps]) {
                    ps = i;
                }
            }
        } else if (state == 2) {        // find troughs
            DEBUG_printf("check low %d %d %d %d %zd\n", i, smoothed[i] , smoothed[ps], accu, scores_.size());
            if (smoothed[i] * 0.8 > smoothed[ps] ) {
                state = 0;
                if (accu > std::max<int>(10, scores_.size() / 4)) break;    // TODO 可以设置更小的值，但需要检测峰之间的距离，并且取少量数据
            } else {
                if (smoothed[i] < smoothed[ps]) {
                    ps = i;
                }
            }
        }
    }
    
    DEBUG_printf("troughs: %d, %d, %d, %d\n", start, end, ps, peak);
    int all = std::accumulate(hist.begin(), hist.end(), 0);
    int sel = std::accumulate(hist.begin()+ps, hist.end(), 0);
    DEBUG_printf("troughs: %d %d, %d\n", accu, sel, all);
    if (sel > all * 3 / 4) {
        size_t accu = 0;
        for (int i = hist.size() - 1; i > 0; --i) {
            accu += hist[i];
            if (accu > scores_.size() / 2) {
                return i*b - 1;
            }
        }
        return -1;

        //return ((start + ps + w)/2 * b - 1)/2;
    } else {
        return (ps + w/2) * b - 1;
    }
}


void AlignmentGraph::QueryInfos::SaveReadInfos(std::ostream& os, int tid, const ReadStore &rs) const {
    const std::string& tname = rs.QueryNameById(tid);
    for (size_t i=0; i<scores_.size(); ++i) {
        const auto & s = scores_[i];
        const std::string& qname = rs.QueryNameById(s.qid);
        bool sel = selected_.find(i) != selected_.end() ? 1 : 0;
        os << tname << " " << qname << " " << s.Weight() << " " << sel << "\n";
    }
}   


double AlignmentGraph::LinkScoreCount(size_t col, size_t row, const Link& link) { 
    //return link.count - cols[col].coverage*std::max(branch_score_*4/(4+row), 0.3); 
    
    double scale = std::max<double>(std::pow<double>(branch_score_[2], row)*branch_score_[0], branch_score_[1]);
    double compensate = std::max<double>(scale * cols[col].coverage, min_coverage_ * branch_score_[0]);
    return link.count - compensate;
} 

double AlignmentGraph::LinkScoreWeight(size_t col, size_t row, const Link &link) {
    double s = 0;
    for (size_t i=0; i<query_infos_.scores_.size(); ++i) {
        if (link.seqs[i+1] && query_infos_.selected_.find(i) != query_infos_.selected_.end()) {
            s += query_infos_.scores_[i].WeightInGraph(score_range_, weight_range_);
        }
    }

    //double compensate = std::pow<double>(reduction_, row) * branch_score_;
    //printf("score: %d %d %zd %f %f %f\n",row, col, cols[col].coverage, compensate, min_coverage_ * branch_score_ / cols[col].coverage, std::max(branch_score_*4/(4+row), 0.3));
    //return s - std::max<double>(compensate*cols[col].weight, min_coverage_ * branch_score_ * cols[col].weight / cols[col].coverage );

    double scale = std::max<double>(std::pow<double>(branch_score_[2], row)*branch_score_[0], branch_score_[1]);
    double compensate = std::max<double>(scale * cols[col].weight, branch_score_[0] * cols[col].weight * min_coverage_ / cols[col].coverage);

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
                    if (l.prev.col != -1)
                        of << i << "_" << ir << "_" << "ACGT-"[ib] << ", "  
                           << l.prev.col << "_" << l.prev.row << "_" << "ACGT-"[l.prev.base] << "," 
                           << l.count << "," 
                           << l.seqs.to_string() << "\n";
                }
            }
        }
    }

}


} // namespace fsa {
