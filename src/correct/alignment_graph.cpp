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

    VerifyImportantBranches0 = dataset_.variants == nullptr ? 
        &AlignmentGraph::VerifyImportantBranches1 : 
        &AlignmentGraph::VerifyImportantBranchesByVariants;
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
    
    if (sopts_.debug_flag == 0) {
        ComputeSimilarity4();
    } else {
        ComputeSimilarity();
    }

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
                cols[col].weight += query_infos_.scores_[i].WeightInGraph(score_range_, weight_range_);
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
            DEBUG_printf("col: (%zd,%zd,%zd) %zd, %zd, %zd\n", loc.col, loc.row, loc.base, 
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

    auto segs = SplitSegment(seg.end);
    DEBUG_printf("seg: size=%zd\n", segs.size());
    for (auto &s : segs) {
        DEBUG_printf("seg: (%zd, %zd) %d\n", s.begin.col, s.end.col, s.type);
    }
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

        if ( !(bt == '-' || (*target_)[index_t] == Base2Num[bt])) {
            LOG(INFO)("Error %c(%d) %d(%zd) == %d", bt, bt, (*target_)[index_t], index_t, Base2Num[bt]);
            LOG(INFO)("SSS: %zd, %zd, %zd", sid, query_start,target_start);
            printf("q:%s\nt:%s\n", aligned_query.c_str(), aligned_target.c_str());
            printf("%s\n", (*target_).ToString()->c_str());
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
                DEBUG_printf("position: %zd\n", i);

                auto tlink = find_link(0);
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

    query_infos_.SelectReads(min_selected);
    
}


void AlignmentGraph::ComputeSimilarity4() {
    assert(cols.size() > 0 );

    auto cands = CollectImportantBranches();
    DEBUG_printf("candidate location: %zd\n", cands.size());

    (this->*VerifyImportantBranches0)(cands);

    auto find_link = [](const std::array<ImportantBranch::LinkCol,2>& links, size_t id) -> const Link*{
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

    query_infos_.SelectReads3(min_selected, range_);
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
            DEBUG_printf("cand: %zd, th=%d, cov=%d, count=%d  %d,%d\n", i, branch_threshold, cols[i].coverage, link_count,  links[0]->count, links[1]->count);
            if (links[1]->count >= (size_t)branch_threshold) {
                cands.push_back(ImportantBranch());
                
                //printf("position0: %zd\n", i);
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
    std::unordered_set<size_t> verified;

    size_t start = 0;
    while (start < cands.size()) {
        // find consecutive locations
        size_t end = start+1;
        for (; end < cands.size(); ++end) {
            if (cands[end-1].links[0].l->prev.col + 2 < cands[end].links[0].l->prev.col) {
                break;
            }
        }

        DEBUG_printf("Verify range(0) %zd-%zd(%zd)\n", cands[start].c, cands[end-1].c, end-start);
        bool check_ok = !sopts_.skip_branch_check ? VerifiyImportantBranch(cands, start, end) : true;
        DEBUG_printf("Verify range(1) %zd-%zd(%zd) %d\n", cands[start].c, cands[end-1].c, end-start, check_ok);
        if (check_ok) {
            assert(end > start);
            
            DEBUG_printf("AddPosition: %zd %zd %zd\n", start, end-1, (start+end-1) / 2);
            verified.insert(start);
            verified.insert(end-1);
            verified.insert((start+end-1) / 2);
            //for (size_t i = start ; i < end; ++i) {
            //    verified.insert(i);
            //}
        }
        start = end;
    }
    
    for (size_t i = 0; i < cands.size(); ++i) {
        if (verified.find(i) == verified.end()) {
            cands[i].valid = false;
        }
    }
    std::unordered_set<size_t> removed;
    VerifyImportantSitesByConsistent(cands, 0.6, removed);
    VerifyImportantSitesByConsistent(cands, 0.7, removed);
    ReactivateImportantSitesByConsistent(cands, 0.7, removed);
    //VerifyImportantBranchesByDensity(cands);
}

void AlignmentGraph::VerifyImportantBranchesByVariants(std::vector<ImportantBranch>& cands) {

    const int N = 3;
    std::unordered_set<int> vars;
    auto stdvars = dataset_.variants->QuerySnps(tid_);
    for (auto s : stdvars) {
        for (int i = s - N; i < s + N + 1; ++i) {
            vars.insert(i);
        }
    }

    for (size_t i = 0; i < cands.size(); ++i) {
        if (vars.find(cands[i].c) == vars.end()) {
            cands[i].valid = false;
        }
    }
}

void PrintBases(const std::string& msg, const std::vector<int> &bs) {
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

std::vector<int> AlignmentGraph::ExtendLeft(const Loc& start, size_t len, const MyBitSet& seqs) {
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

std::vector<int> AlignmentGraph::ExtendRight(const Loc &start, size_t len, const MyBitSet& seqs) {
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
    DEBUG_printf("Verange: %zd %zd - %zd %zd\n", tstart, tend, range_[0], range_[1]);
    if (tstart < range_[0] + N || tend >= range_[1] - N) {
        return false;
    }

    std::vector<int> local(2*N + tend - tstart);
    for (size_t i = 0; i < 2*N + tend - tstart; ++i) {
        local[i] = (*target_)[tstart-N + i];
    }
    PrintBases("target", local);
    
    if (HasHomopolymer(4, local, N-3, N+tend-tstart+3)) {
        return false;
    }

    DEBUG_printf("pass homo\n");
    // 终点合并
    // if (!(cands[end-1].links[0].l->prev != cands[end-1].links[1].l->prev && cands[end-1].links[0].r == cands[end-1].links[1].r && cands[end-1].links[0].r != 4)) {
    //     return false;
    // }
    // DEBUG_printf("pass end\n");

    auto extend_left = [this](const Link* link, const Loc& right, size_t left, Loc& end) {
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
    };

    Loc left0;
    auto branch0 = extend_left(cands[end-1].links[0].l, Loc(cands[end-1].c, 0, cands[end-1].links[0].r), cands[start].c-1, left0);
    if (branch0.size() == 0) return false;
    PrintBases("branch0", branch0);

    Loc left1;
    auto branch1 = extend_left(cands[end-1].links[1].l, Loc(cands[end-1].c, 0, cands[end-1].links[1].r), cands[start].c-1, left1);
    if (branch1.size() == 0) return false;
    PrintBases("branch1", branch1);
    
    DEBUG_printf("pass extend\n");

    auto brl0 = ExtendLeft(left0, 3, cands[start].links[0].l->seqs);
    auto brl1 = ExtendLeft(left1, 3, cands[start].links[1].l->seqs);
    PrintBases("brl0", brl0);
    PrintBases("brl1", brl1);
    if (brl0.size() < 3 || brl1.size() < 3) return false;

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
    auto ibrl = find_matched(brl0, brl1);
    if (ibrl[0] == brl0.size() || ibrl[1] == brl1.size() ) return false;
    
    auto brr0 = ExtendRight({(int)cands[end-1].c, 0, cands[end-1].links[0].r}, 3, cands[end-1].links[0].l->seqs);
    PrintBases("brr0", brr0);
    auto brr1 = ExtendRight({(int)cands[end-1].c, 0, cands[end-1].links[1].r}, 3, cands[end-1].links[1].l->seqs);
    PrintBases("brr1", brr1);
    if (brr0.size() < 3 || brr1.size() < 3) return false;

    auto ibrr = find_matched(brr0, brr1);
    if (ibrr[0] == brr0.size() || ibrr[1] == brr1.size() ) return false;
    
    DEBUG_printf("pass find_matched: %zd,%zd   %zd,%zd\n", ibrl[0], ibrl[1], ibrr[0], ibrr[1]);

    auto has_homo = [](const std::vector<int>& head, const std::vector<int>& body, const std::vector<int> &tail) {
        std::vector<int> s(head.begin(), head.end());
        s.insert(s.end(), body.begin(), body.end());
        s.insert(s.end(), tail.begin(), tail.end());
        return HasHomopolymer(4, s, 0, s.size());
    };
    if (has_homo(brl0, branch0, brr0) || has_homo(brl1, branch1, brr1)) return false;

    DEBUG_printf("pass has_homo\n");

    std::vector<int> bubble0(brl0.begin()+ibrl[0], brl0.end());
    bubble0.insert(bubble0.end(), branch0.begin(), branch0.end());
    bubble0.insert(bubble0.end(), brr0.begin(), brr0.begin()+ibrr[0]+2);
    
    std::vector<int> bubble1(brl1.begin()+ibrl[1], brl1.end());
    bubble1.insert(bubble1.end(), branch1.begin(), branch1.end());
    bubble1.insert(bubble1.end(), brr1.begin(), brr1.begin()+ibrr[1]+2);
    PrintBases("bubble0", bubble0);
    PrintBases("bubble1", bubble1);
    // //assert(branch0.size() >= 2 && branch1.size() >= 2);
    // if (branch0.size() < 2 && branch1.size() < 2) return false;
    // if (branch0.front() != branch1.front() || branch0.back() != branch1.back())  return false;
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

        PrintBases("mid0", std::vector<int>(branch0.begin()+ihead, branch0.end()-itail));
        PrintBases("mid1", std::vector<int>(branch1.begin()+ihead, branch1.end()-itail));
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
  
    DEBUG_printf("pass insert0\n");
    if (branch0.size() <= branch1.size()) {
        auto s = std::search(branch1.begin(), branch1.end(), branch0.begin(), branch0.end());
        if (s != branch1.end()) return false;
    } else {
        auto s = std::search(branch0.begin(), branch0.end(), branch1.begin(), branch1.end());
        if (s != branch0.end()) return false;
    }
    if (branch0.size() == 2 || branch1.size() == 2) return false;

    DEBUG_printf("pass insert\n");
    return true;
}

void AlignmentGraph::VerifyImportantSitesByConsistent(std::vector<ImportantBranch> &cands, double th, std::unordered_set<size_t> &removed) {
    //std::unordered_set<size_t> removed;

    for (size_t i = 0; i < cands.size(); ++i) {
        if (!cands[i].valid) continue;

        size_t icol = cands[i].c; 
        std::array<size_t, 2> scores = {0, 0};

        for (size_t j = 0; j < cands.size(); ++j) {
            if (!cands[j].valid) continue;

            size_t jcol = cands[j].c; 
            if (icol < jcol + 5 && jcol < icol + 5) continue;
            std::array<size_t, 4> ss = {0, 0, 0, 0};
            ss[0] = (cands[i].links[0].l->seqs & cands[j].links[0].l->seqs).count();
            ss[1] = (cands[i].links[0].l->seqs & cands[j].links[1].l->seqs).count();
            ss[2] = (cands[i].links[1].l->seqs & cands[j].links[0].l->seqs).count();
            ss[3] = (cands[i].links[1].l->seqs & cands[j].links[1].l->seqs).count();
            DEBUG_printf("consistent(%zd,%zd) %zd,%zd,%zd,%zd\n", icol, jcol, ss[0], ss[1],ss[2],ss[3]);

            if (ss[0] + ss[1] + ss[2] + ss[3] < 10) continue;
            if (ss[0] + ss[3] >= ss[1] + ss[2]) {
                scores[0] += ss[0] + ss[3];
                scores[1] += ss[1] + ss[2];
            } else {
                scores[0] += ss[1] + ss[2];
                scores[1] += ss[0] + ss[3];
            }

        } 
        DEBUG_printf("consistent pos=%zd, %zd / %zd = %0.2f < %0.2f\n", icol, scores[0], scores[1], scores[0]*1.0/(scores[0]+scores[1]), th);
        if (scores[0]*1.0/(scores[0]+scores[1]) < th) {
            removed.insert(i);
        }
    }

    for (auto i : removed) {
        cands[i].valid = false;
    }
}



void AlignmentGraph::ReactivateImportantSitesByConsistent(std::vector<ImportantBranch> &cands, double th, const std::unordered_set<size_t> &removed) {
    std::unordered_set<size_t> restored;

    for (auto i : removed) {

        size_t icol = cands[i].c; 
        std::array<size_t, 2> scores = {0, 0};

        for (size_t j = 0; j < cands.size(); ++j) {
            if (!cands[j].valid) continue;

            size_t jcol = cands[j].c; 
            if (icol < jcol + 5 && jcol < icol + 5) continue;
            std::array<size_t, 4> ss = {0, 0, 0, 0};
            ss[0] = (cands[i].links[0].l->seqs & cands[j].links[0].l->seqs).count();
            ss[1] = (cands[i].links[0].l->seqs & cands[j].links[1].l->seqs).count();
            ss[2] = (cands[i].links[1].l->seqs & cands[j].links[0].l->seqs).count();
            ss[3] = (cands[i].links[1].l->seqs & cands[j].links[1].l->seqs).count();
            DEBUG_printf("consistent(%zd,%zd) %zd,%zd,%zd,%zd\n", icol, jcol, ss[0], ss[1],ss[2],ss[3]);

            if (ss[0] + ss[1] + ss[2] + ss[3] < 10) continue;
            if (ss[0] + ss[3] >= ss[1] + ss[2]) {
                scores[0] += ss[0] + ss[3];
                scores[1] += ss[1] + ss[2];
            } else {
                scores[0] += ss[1] + ss[2];
                scores[1] += ss[0] + ss[3];
            }

        } 
        DEBUG_printf("consistent pos=%zd, %zd / %zd = %0.2f < %0.2f\n", icol, scores[0], scores[1], scores[0]*1.0/(scores[0]+scores[1]), th);
        if (scores[0]*1.0/(scores[0]+scores[1]) >= th) {
            restored.insert(i);
        }
    }

    for (auto i : restored) {
        cands[i].valid = true;
    }
}

void AlignmentGraph::VerifyImportantBranchesByDensity(std::vector<ImportantBranch> &cands) {
    std::vector<size_t> density(cands.size(), 0);

    const size_t D = 10000;
    std::list<size_t> position;

    for (size_t i = 0; i < cands.size(); ++i) {
        if (!cands[i].valid) continue;

        if (position.size() > 0 && cands[i].c - cands[position.front()].c > D) {
            for (auto p : position) {
                if (position.size() > density[p] ) {
                    density[p] = position.size();
                }
            }
            while(position.size() > 0 && cands[i].c - cands[position.front()].c > D) {
                position.pop_front();
            }
        }
        assert (position.size() == 0 || cands[i].c - cands[position.front()].c <= D);
        position.push_back(i);
    }

    for (auto p : position) {
        if (position.size() > density[p] ) {
            density[p] = position.size();
        }
    }

    // for (size_t i = 0; i < cands.size(); ++i) {
    //     if (!cands[i].valid) continue;
    //     density[i]++;
    //     for (size_t j = i+1; j < cands.size(); ++j) {
    //         if (cands[j].c - cands[i].c <= D) {
    //             if (cands[j].valid) {
    //                 density[i]++;
    //                 density[j]++;
    //             }
    //         } else {
    //             break;
    //         }
    //     }
    // }

    const size_t C = 6;

    for (size_t i = 0; i < cands.size(); ++i) {
        if (cands[i].valid) {
            DEBUG_printf("density: %zd, %zd, %zd\n", i, cands[i].c, density[i]);
            if (density[i] < C) {
                cands[i].valid = false;
            }
        }
    }
}
    

void AlignmentGraph::QueryInfos::SelectReads(int min_sel) {

    // select queries for correction
    std::vector<int> score_index(scores_.size());
    for (size_t i = 0; i < score_index.size(); ++i) {
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

void AlignmentGraph::QueryInfos::SelectReads2(size_t min_sel, const std::array<size_t,2> &range) {

    // select queries for correction
    std::vector<double> score_left;
    std::vector<double> score_right;
    std::vector<double> score_all;

    size_t mid = (range[0] + range[1]) / 2;
    size_t len = (range[1] - range[0]) / 2;
    size_t left_range = mid - len / 2;
    size_t right_range = mid + len / 2;

    size_t comm = 0;
    for (const auto &s : scores_) {
        auto m = (s.tstart + s.tend) / 2;
        auto w = s.Weight();
        int c = 0;
        if (m <= mid || s.tstart <= left_range) {
            score_left.push_back(w);
            //printf("w-left: %f\n", w);
            c++;
        } 
        if (m >= mid || s.tend >= right_range) {
            score_right.push_back(w);
            //printf("w-right: %f\n", w);
            c++;
        }
        if (c == 2) comm ++;
        score_all.push_back(w);
        //printf("w-all: %f\n", w);
    }

    //printf("comm: %zd / %zd\n", comm, scores_.size());
    double th_left = 2.0;
    double th_right = 2.0;
    if (comm < scores_.size() /2 ) {
        th_left = std::max(0.0, FindScoreThreshold3(score_left)-0.000001);
        th_right = std::max(0.0, FindScoreThreshold3(score_right)-0.000001);
    } 
    double th_all = std::max(0.0, FindScoreThreshold3(score_all));


    //printf("th: %f %f %f\n", th_left, th_right, th_all);

    //assert(selected_.size() == 0);
    selected_.clear();
    for (size_t i = 0; i < scores_.size(); ++i) {
        const auto &s = scores_[i];
        auto m = (s.tstart + s.tend) / 2;
        auto w = s.Weight();
        if ((m <= mid || s.tstart <= left_range) && w > th_left) {
            selected_.insert(i);
        } 
        if ((m >= mid || s.tend >= right_range) && w > th_right) {
            selected_.insert(i);
        }
        if (w >= th_all) {
            selected_.insert(i);
        }
    }

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


double AlignmentGraph::LinkScoreCount(size_t col, size_t row, const Link& link) { 
    //return link.count - cols[col].coverage*std::max(branch_score_*4/(4+row), 0.3); 
    
    double scale = std::max<double>(std::pow<double>(branch_score_[2], row)*branch_score_[0], branch_score_[1]);
    double compensate = std::max<double>(scale * cols[col].coverage, sopts_.min_coverage * branch_score_[0]);
    return link.count - compensate;
} 

double AlignmentGraph::LinkScoreWeight(size_t col, size_t row, const Link &link) {
    double s = 0;
    for (size_t i=0; i<query_infos_.scores_.size(); ++i) {
        if (link.seqs[i+1] && query_infos_.selected_.find(i) != query_infos_.selected_.end()) {
            s += query_infos_.scores_[i].WeightInGraph(score_range_, weight_range_);
        }
    }
    if (link.seqs[0]) {
        s += 0.5; // TODO Target score 
    }

    //double compensate = std::pow<double>(reduction_, row) * branch_score_;
    //printf("score: %d %d %zd %f %f %f\n",row, col, cols[col].coverage, compensate, min_coverage_ * branch_score_ / cols[col].coverage, std::max(branch_score_*4/(4+row), 0.3));
    //return s - std::max<double>(compensate*cols[col].weight, min_coverage_ * branch_score_ * cols[col].weight / cols[col].coverage );

    double scale = std::max<double>(std::pow<double>(branch_score_[2], row)*branch_score_[0], branch_score_[1]);
    double compensate = std::max<double>(scale * cols[col].weight, branch_score_[0] * cols[col].weight * sopts_.min_coverage / cols[col].coverage);
    
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
