#include "alignment_array.hpp"

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
#include "logger.hpp"

namespace fsa {

DnaSerialTable2 AlignmentArray::Base2Num;

AlignmentArray::AlignmentArray() {
}


void AlignmentArray::SetParameter(const std::string &name, const std::string &v) {
    if (name == "score") {
    } else {
        LOG(FATAL)("Not support parameter: %s", name.c_str());
    }
}

void AlignmentArray::SetParameter(const std::string &name, int v) {
    LOG(FATAL)("Not support parameter: %s", name.c_str());
}

void AlignmentArray::SetParameter(const std::string &name, double v) {
    LOG(FATAL)("Not support parameter: %s", name.c_str());
}

TimeCounter nodes_tc("InitializeNodes");
void AlignmentArray::InitializeNodes(const DnaSeq &target, const std::vector<Alignment> &aligned) {
    TimeCounter::Mark mark(nodes_tc);
    cols_.assign(target.Size()+1, Column());

    for (size_t i=0; i<cols_.size() - 1; ++i) {
        cols_[i].index = 1;     // for target
    }

    for (const auto &al : aligned) {
        int offt = al.target_start - 1;
        int row = 0;
        for (size_t i=0; i<al.AlignSize(); ++i) {
            auto bs = al.GetAlign(i); // {q,t}

            if (bs[1] == bs[0] && bs[1] != '-') {
                ++offt;
                row = 0;
            } else if (bs[1] != '-' && bs[0] == '-') {
                ++offt;
                row = 0;
            } else if (bs[1] == '-' && bs[0] != '-') {
                ++row;
            } else {
                if (!(bs[1] == '-' && bs[0] == '-')) {
                    printf(" %c %c\n", bs[1], bs[0]);
                    fflush(stdout);
                }
                assert(bs[1] == '-' && bs[0] == '-');
                continue;
            }
            cols_[offt].index = std::max(cols_[offt].index, row+1);
        }
    }

    for (size_t i=1; i<cols_.size(); ++i) {
        cols_[i].index += cols_[i-1].index;
    }

    for (size_t i=cols_.size(); i>=2; --i) {
        cols_[i-1].index = cols_[i-2].index;
    }
    cols_[0].index = 0;
    nodegroups_.assign(cols_.back().index, NodeGroup());
    nodes_tc.Add(cols_.back().index);

    //for (size_t i=0; i<cols_.size()-1; ++i) {
    //    cols_[i].rows = &nodegroups_[cols_[i].index];
    //    cols_[i].size = cols_[i+1].index - cols_[i].index;
    //}
}


TimeCounter build_tc("Build");
void AlignmentArray::Build(const DnaSeq& target, const std::array<size_t,2> &range, const std::vector<Alignment> &aligned, const std::vector<int> &qids) {
    TimeCounter::Mark mark(build_tc);
    Clear();
    assert(aligned.size() == qids.size());

    InitializeNodes(target, aligned);

    AddTarget(target);
    range_ = range;

    for (size_t i=0; i<aligned.size(); ++i) {
        AddQuery(i, aligned[i]);
    }
}


void AlignmentArray::Clear() {
    cols_.clear();
    target_ = nullptr;
    ranges_.clear();
    scores_.clear();
    sequence_.clear();
}

void AlignmentArray::Save(const std::string &fname) const {
    std::ofstream out(fname);

    // for (size_t i = 0; i < cols_.size(); ++i) {
    //     for (size_t j = 0; j < cols_[i].rows.size(); ++j) {
    //         for (size_t k = 0; k < 5; ++k) {
    //             //const Node &node = cols[i].rows[j].base[k];

    //             //const char* acgt = ".ACGT-";
    //         }
    //     }
    // }
}

TimeCounter path_tc("FindBestPath");

AlignmentArray::Loc AlignmentArray::FindBestPath() {
    TimeCounter::Mark mark(path_tc);
    Loc end = {-1,-1,-1};
    double global_score = -1;

    std::vector<int> prev_links;
    for (size_t i = 0; i < ColSize(); i++) {
        for (size_t j = 0; j < RowSize(i); j++) {
            for (size_t k = 0; k < NodeSize(); ++k) {

                Node& node = GetNode(i, j, k);
                if (node.seqs.count() > 0) {
                    if (j == 0 && i > 0) {
                        for (size_t jj = 0; jj < RowSize(i-1); ++jj) {
                            for (size_t kk = 0; kk < NodeGroup::Size(); ++kk) {
                                Loc ploc = {i-1, (int)jj, (int)kk};
                                
                                const Node& pnode = GetNode(i-1, jj, kk);
                                if (pnode.seqs.count() > 0) {
                                    size_t c = (pnode.seqs & node.seqs).count();

                                    c = std::max(0, (int)(c - prev_links[jj*NodeGroup::Size()+kk]));

                                    double w = row_weight_*4/(4+(jj>0?(jj-1):0));
                                    double score = c - cols_[i].coverage*w + pnode.score;

                                    if (score > node.score) {
                                        node.score = score;
                                        node.best = ploc;
                                    }
                                }
                            }
                        }
                    } else if ( j > 0 ){
                        for (size_t kk = 0; kk < NodeGroup::Size(); ++kk) {
                            Loc ploc = {i, (int)(j-1), (int)kk};

                            const Node& pnode = GetNode(i, j-1, kk);
                            if (pnode.seqs.count() > 0) {
                                size_t c = (pnode.seqs & node.seqs).count();
                                prev_links[(j-1)*NodeGroup::Size()+kk] += c;

                                double w = row_weight_*4/(4+j-1);

                                double score = c - cols_[i].coverage*w + pnode.score;
                                if (score > node.score) {
                                    node.score = score;
                                    node.best = ploc;

                                }
                            }
                        }
                    }
                
                    if (node.score > global_score) {
                        global_score = node.score;
                        end = {i, j, k};
                    }
                }
            }

            if (j == 0) {
                prev_links.assign(RowSize(i)*NodeGroup::Size(), 0);
            }
        }
    }
    
    return end;

}

AlignmentArray::Loc AlignmentArray::FindBestPath1() {
    TimeCounter::Mark mark(path_tc);
    Loc end = {-1,-1,-1};
    double global_score = -1;

    // 

    std::vector<std::array<std::array<int,5>, 5>> wts_0;
    std::vector<std::array<std::array<int,5>, 5>> wts_1;

    //LOG(INFO)("IN");

    for (size_t i = 1; i < ColSize(); i++) {
        Column& col = GetCol(i-1);
        if (col.coverage <= 1) continue;

        const size_t rowsize = RowSize(i-1);

        wts_0.assign(RowSize(i-1)-1+1, {0});
        wts_1.assign(RowSize(i-1)+1, {0});

        for (size_t j=1; j<RowSize(i-1); ++j) {
            for (size_t k=0; k<5; ++k) {
                Node& node = GetNode(i-1, j, k);
                if (node.seqs.count() > 0) {
                    for (size_t kk = 0; kk < 5; ++kk) {
                        Node& pnode = GetNode(i-1, j-1, kk);
                        wts_0[j-1][k][kk] = (node.seqs & pnode.seqs).count();
                    }
                }

            }
        }

        for (size_t k=0; k<5; ++k) {
            Node& node = GetNode(i, 0, k);
            if (node.seqs.count() > 0) {
                for (size_t j=0; j<RowSize(i-1); ++j) {
                    for (size_t kk = 0; kk < 5; ++kk) {
                        Node& pnode = GetNode(i-1, j, kk);
                        if (pnode.seqs.count() > 0) {
                             wts_1[j][k][kk] = (node.seqs & pnode.seqs).count();
                        }
                    }
                }
            }
        }

        // modification
        for (size_t j=0; j < RowSize(i-1)-1; ++j) {
            for (size_t k=0; k<5; ++k) {
                for (size_t kk=0; kk<5; ++kk) {
                    auto w = (wts_0[j][0][kk] + wts_0[j][1][kk] + wts_0[j][2][kk] + wts_0[j][3][kk] + wts_0[j][4][kk]);
                    wts_1[j][k][kk] = std::max(0, wts_1[j][k][kk] - w);
                }
            }
        }

        // merge nodes
        for (size_t k = 0; k < 5; ++k) {
            Node& node = GetNode(i, 0, k);
            if (node.seqs.count() > 0) {
                for (size_t kk = 0; kk < 5; ++kk) {
                    std::vector<size_t> merged;
                    int wa = 0;
                    for (size_t j = 0; j < RowSize(i-1); ++j) {
                        if (wts_1[j][k][kk] > 0) {
                            merged.push_back(j);
                            wa += wts_1[j][k][kk];
                        }
                    }

                    for (auto m : merged) {
                        wts_1[m][k][kk] = wa;
                    }

                    
                }

            }
            
        }

        // find best
        for (size_t j = 1; j<RowSize(i-1); ++j) {
            for (size_t k = 0; k < 5; ++k) {
               Node& node = GetNode(i-1, j, k);
               if (node.seqs.count() > 0) {
                   for (size_t kk=0; kk<5; ++kk) {

                       auto c = wts_0[j-1][k][kk];
                        if (c > 0) {
                            Node& pnode = GetNode(i-1, j-1, kk);
                            
                            double w = row_weight_*4/(4+j-1);

                            double score = c - cols_[i-1].coverage*w + pnode.score;
                            if (score > node.score) {
                                node.score = score;
                                node.best = {i-1, j-1, kk};
                            }
                        }
                   }
               }
                
                if (node.score > global_score) {
                    global_score = node.score;
                    end = {i, j, k};
                }
            }
        } 

        for (size_t k = 0; k < 5; ++k) {
            Node& node = GetNode(i, 0, k);

            if (node.seqs.count() > 0) {
                for (size_t j=0; j <RowSize(i-1); ++j) {
                    for (size_t kk = 0; kk < 5; ++kk) {
                        auto c = wts_1[j][k][kk];
                        if (c > 0) {
                            Node& pnode = GetNode(i-1, j, kk);
                            
                            double w = row_weight_*4/(4+j-1);

                            double score = c - cols_[i-1].coverage*w + pnode.score;
                            if (score > node.score) {
                                node.score = score;
                                node.best = {i-1, j, kk};
                            }

                        }
                    }
                }
                if (node.score > global_score) {
                    global_score = node.score;
                    end = {i, 0, k};
                }
            }

        }

    }
    
    return end;

}

AlignmentArray::Loc AlignmentArray::FindBestPathWithCount() {
    Loc end = {-1,-1,-1};
    double global_score = -1;

    Loc ploc = {-1, -1, -1};
    for (size_t i = 0; i < ColSize(); i++) {
        Column& col = GetCol(i);
        if (col.coverage <= 1) continue;

        for (size_t j=0; j<RowSize(i); ++j) {
            for (size_t k=0; k<5; ++k) {
                Node &node = GetNode(i,j,k);
                if (node.seqs.count() >= col.coverage / 2) {
                    node.best = ploc;
                    ploc = {i, j, k};
                    break;
                }
            }
            if (ploc.row != j) break;
        }

       
    }
    
    return ploc;

}

TimeCounter recon_tc("Reconstruct");
void AlignmentArray::Reconstruct(const Loc &end) {
    TimeCounter::Mark mark(recon_tc);
    assert(sequence_.size() == 0);
    sequence_.reserve(int(cols_.size()*1.2));

    for (Loc loc = end; loc.col >= 0; loc = Get(loc)->best) {
        char b = "ACGT-"[loc.base];

        if (b != '-') {
            sequence_.push_back(b);
        }
    }

    std::reverse(sequence_.begin(), sequence_.end());
}

TimeCounter cns_tc("Consensus");
void AlignmentArray::Consensus() {
    TimeCounter::Mark mark(cns_tc);
    Reconstruct(FindBestPath());
}


void AlignmentArray::AddTarget(const DnaSeq& target) {

    target_ = &target;

    for (size_t i = 0; i < target.Size(); ++i) {
        GetCol(i).AddSeq(0);
    }
}

TimeCounter query_tc("AddQuery");
void AlignmentArray::AddQuery(size_t sid, const Alignment &al) {
    TimeCounter::Mark mark(query_tc);
    size_t row = 0;
    int offset_t = (int)al.target_start - 1;   // move the index to -1 for the first operation ++.
    int offset_q = (int)al.query_start - 1;
    

    for (size_t i=0; i<al.AlignSize(); ++i) {
    ///printf("q1:%s-\nt1:%s-\n", al.aligned_query.c_str(),  al.aligned_target.c_str());fflush(stdout);

        auto bs = al.GetAlign(i); // {bq, bt}
        if (bs[1] == bs[0] && bs[1]  != '-') {
            ++offset_t;
            ++offset_q;
            row = 0;
        } else if (bs[0] == '-' && bs[1]  != '-') {
            ++offset_t;
            row = 0;
        } else if (bs[1]  == '-' && bs[0] != '-') {
            ++offset_q;
            ++row;
        } else {
                if (!(bs[1] == '-' && bs[0] == '-')) {
                    printf(" %c %c\n", bs[1], bs[0]);
                    fflush(stdout);
                }
            assert(bs[1]  == '-' && bs[0] == '-');
            continue;
        }

        // assert(bt == '-' || (*target_)[offset_t] == Base2Num[bt]);

        if (offset_t < (int)range_[0]) continue;
        if (offset_t >= (int)range_[1]) break;

        // assert(offset_t >= 0);
        Loc curr_loc = { offset_t, (int)row, Base2Num[bs[0]] };
        if (row == 0) {
            GetCol(offset_t).AddSeq(sid+1);
        }
        Get(curr_loc)->AddSeq(sid+1);

    }
}



double AlignmentArray::LinkScoreWeight(size_t col, size_t row) {
    double s = 0;
    double c = 0;
    //if (link.seqs[0]) s += 1;

    assert(0);
    return 0;
    //return s-(branch_score_*4/(4+row))*cols_[col].weight;

}


} // namespace fsa {
