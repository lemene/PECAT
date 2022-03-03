#include "edlib_aligner.hpp"

#include "edlib.h"
#include "ksw2.h"
#include "./utils/logger.hpp"
#include "utility.hpp"

namespace fsa {

EdlibAligner::EdlibAligner(const std::vector<std::string> &opts) {
    // edlib:bs=1000

    assert(opts[0] == "edlib");
    for (size_t i=1; i<opts.size(); ++i) {
        auto kv = SplitStringByChar(opts[i], '=');
        if (kv[0] == "bs") {
            block_size_ = std::stoi(kv[1]);
        } else if (kv[0] == "mc") {
            match_count_ = std::stoi(kv[1]);
        } else if (kv[0] == "be") {
            block_error_ = std::stod(kv[1]);
        } else {
            LOG(ERROR)("Not support parameter: aligner=edlib:..:%s:..", opts[i].c_str());
        }
    }
}

EdlibAligner::~EdlibAligner() {

}

void Print(const char* qseq, const char* tseq, size_t N) {
    printf("q:");
    for (size_t i=0; i<N; ++i) {
        printf("%c", "ACGT"[qseq[i]]);
    }
    printf("\nt:");
    for (size_t i=0; i<N; ++i) {
        printf("%c", "ACGT"[tseq[i]]);
    }
    printf("\n\n");
}

void PrintSeq(const char* seq, size_t N, const char* msg) {
    printf("%s:", msg);
    for (size_t i=0; i<N; ++i) {
        printf("%c", "ACGT"[seq[i]]);
    }
    printf("\n\n");
}

bool CheckAlign(const char* qseq, size_t qsize, const char* tseq, size_t tsize, size_t qstart, size_t tstart, std::vector<unsigned char> alignment) {
    size_t qi = qstart;
    size_t ti = tstart;
    for (size_t i=0; i < alignment.size(); ++i) {
        //printf("%zd:%c %zd:%c  %zd:%zd:%d \n", qi, "ACGT"[qseq[qi]] , ti, "ACGT"[tseq[ti]], alignment.size(), i, alignment[i] );
        if (alignment[i] == EDLIB_EDOP_MATCH) {
            assert(qseq[qi] == tseq[ti]);
            qi ++;
            ti ++;
        } else if (alignment[i] == EDLIB_EDOP_INSERT) {
            qi ++;
        } else if (alignment[i] == EDLIB_EDOP_DELETE) {
            ti ++;
        } else {
            assert(alignment[i] == EDLIB_EDOP_MISMATCH);
            assert(qseq[qi] != tseq[ti]);
            qi ++;
            ti ++;
        }
    }
}

bool EdlibAligner::ExtendRight(const char* qseq, size_t qsize, const char* tseq, size_t tsize, size_t qstart, size_t tstart, AlignResult& result) {
    size_t qindex = qstart;
    size_t tindex = tstart;
    result.qstart = qstart;
    result.qend = qstart;
    result.tstart = tstart;
    result.tend = tstart;

    int try_count = 1;
    int max_try_count = 1;
    while (true) {
        if (try_count % 2 == 1) {

            size_t qbsize = std::min(block_size_, qsize - qindex);
            size_t tbsize = std::min((size_t)(qbsize*1.5), tsize - tindex);
            if (qbsize*1 > tbsize) {
                qbsize = tbsize / 1;
            }

            DEBUG_printf("edlib tqsize: %zd %zd, %zd, %zd\n", qindex, tindex, qbsize, tbsize);

            if (qbsize == 0 || tbsize == 0) break;

            auto r = edlibAlign(qseq+qindex, qbsize, tseq+tindex, tbsize,
                edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
            
            if (r.status == EDLIB_STATUS_OK) {
                assert(r.numLocations >= 1);
                assert(r.startLocations[0] == 0);

                size_t qbend = qbsize - 1;
                size_t tbend = r.endLocations[0];   // 不需要 -1

                DEBUG_printf("edlib end %zd %zd %d\n", qbend, tbend, r.editDistance);
                int distance = r.editDistance;
                int match = 0; 
                int alend = r.alignmentLength - 1;
                for (; alend >= 0; --alend) {        
                    if (r.alignment[alend] == EDLIB_EDOP_MATCH) {
                        match += 1;
                        qbend --;
                        tbend --;
                    } else if (r.alignment[alend] == EDLIB_EDOP_INSERT) {
                        match = 0;
                        qbend --;
                        distance --; 
                    } else if (r.alignment[alend] == EDLIB_EDOP_DELETE) {
                        match = 0;
                        tbend --;
                        distance --;
                    } else {
                        assert(r.alignment[alend] == EDLIB_EDOP_MISMATCH);
                        match = 0;
                        qbend --;
                        tbend --;
                        distance --;
                    }

                    if (match >= match_count_ ) {
                        DEBUG_printf("edlib check distance %d %d %d %f\n", match, distance, tbend + 1 + match, distance*1.0/(tbend + 1 + match));
                    }

                    if (match >= match_count_ && distance*1.0/(tbend + 1 + match) <= block_error_) {
                        break;
                    }
                }
                        
                DEBUG_printf("edlib check end %d %d %d\n", match, distance, tbend );
                if ((int)tbend <= 0) {
                    std::string subq, subt;
                    AlignedString(r.alignment, r.alignmentLength, qseq+qindex, tseq+tindex, subq, subt); 
                    DEBUG_printf("qsub:%s\ntsub:%s\n", subq.c_str(), subt.c_str());
                }

                if (alend > 0 && match >= match_count_) {
                    result.alignment.insert(result.alignment.end(), r.alignment, r.alignment + alend +  match);
                    result.distance += distance;
                    tindex += tbend + 1 + match;
                    qindex += qbend + 1 + match;
                    DEBUG_printf("edlib distance: %d, %d %d, %f\n", distance, qbend + 1 + match, tbend + 1 + match, distance*1.0/(tbend + 1 + match));
                } else {
                    if (try_count >= 3) break;
                    try_count ++;
                }
            } else {
                break;
            }
            
            edlibFreeAlignResult(r);
        } else {

            size_t tbsize = std::min(block_size_, tsize - tindex);
            size_t qbsize = std::min((size_t)(tbsize*1.5), qsize - qindex);
            if (tbsize*1 > qbsize) {
                tbsize = qbsize / 1;
            }

            DEBUG_printf("edlib tqsize: %zd %zd, %zd, %zd\n", qindex, tindex, qbsize, tbsize);

            if (qbsize == 0 || tbsize == 0) break;

            auto r = edlibAlign(tseq+tindex, tbsize, qseq+qindex, qbsize, 
                edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
            
            if (r.status == EDLIB_STATUS_OK) {
                assert(r.numLocations >= 1);
                assert(r.startLocations[0] == 0);

                size_t tbend = tbsize - 1;
                size_t qbend = r.endLocations[0];   // 不需要 -1

                
                //DEBUG_printf("edlib check %d %d %d %d\n", aaa[0], aaa[1], aaa[2], aaa[3]);
                DEBUG_printf("edlib end %zd %zd %d\n", qbend, tbend, r.editDistance);
                int distance = r.editDistance;
                int match = 0; 
                int alend = r.alignmentLength - 1;
                for (; alend >= 0; --alend) {        
                    if (r.alignment[alend] == EDLIB_EDOP_MATCH) {
                        match += 1;
                        qbend --;
                        tbend --;
                    } else if (r.alignment[alend] == EDLIB_EDOP_INSERT) {
                        match = 0;
                        tbend --;
                        distance --; 
                    } else if (r.alignment[alend] == EDLIB_EDOP_DELETE) {
                        match = 0;
                        qbend --;
                        distance --;
                    } else {
                        assert(r.alignment[alend] == EDLIB_EDOP_MISMATCH);
                        match = 0;
                        qbend --;
                        tbend --;
                        distance --;
                    }

                    if (match >= match_count_ ) {
                        DEBUG_printf("edlib check distance %d %d %d %f\n", match, distance, tbend + 1 + match, distance*1.0/(tbend + 1 + match));
                    }

                    if (match >= match_count_ && distance*1.0/(tbend + 1 + match) <= block_error_) {
                        break;
                    }
                }
                        
                DEBUG_printf("edlib check end %d %d %d\n", match, distance, tbend );
                if ((int)tbend <= 0) {
                    std::string subq, subt;
                    AlignedString(r.alignment, r.alignmentLength, qseq+qindex, tseq+tindex, subq, subt); 
                    DEBUG_printf("qsub:%s\ntsub:%s\n", subq.c_str(), subt.c_str());
                }

                if (alend > 0 && match >= match_count_) {
                    for (int ii = 0; ii < alend+match; ++ii ) {
                        if (r.alignment[ii] == EDLIB_EDOP_INSERT) {
                            r.alignment[ii] = EDLIB_EDOP_DELETE;
                        } else if (r.alignment[ii] == EDLIB_EDOP_DELETE) {
                            
                            r.alignment[ii] = EDLIB_EDOP_INSERT;
                        }
                    }
                    result.alignment.insert(result.alignment.end(), r.alignment, r.alignment + alend +  match);
                    result.distance += distance;
                    tindex += tbend + 1 + match;
                    qindex += qbend + 1 + match;
                    DEBUG_printf("edlib distance: %d, %d %d, %f\n", distance, qbend + 1 + match, tbend + 1 + match, distance*1.0/(tbend + 1 + match));
                } else {
                    if (try_count >= 3) break;
                    try_count ++;
                }
            } else {
                break;
            }
            
            edlibFreeAlignResult(r);
        }
        result.qstart = qstart;
        result.qend = qindex;
        result.tstart = tstart;
        result.tend = tindex;

    }
    DEBUG_printf("edlib extend: %zd, %d %d - %zd, %d %d\n", qsize, result.qstart, result.qend, tsize, result.tstart, result.tend);
    //CheckAlign(qseq, qsize, tseq, tsize, result.qstart, result.tstart, result.alignment);
}

bool EdlibAligner::Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t, 2> qrange, std::array<size_t, 2> trange, Alignment& al) {
    //LOG(INFO)("O EdlibAligner::AlignBlock");
    //Alignment al0;
    //AlignBlock(qseq, qsize, tseq, tsize, qrange, trange, al0);
    AlignResult result_right;
    ExtendRight(qseq, qsize, tseq, tsize, qrange[0], trange[0], result_right);

    AlignResult result_left;
    std::string qseq_left(qseq, qseq+qrange[0]);
    std::reverse(qseq_left.begin(), qseq_left.end());

    std::string tseq_left(tseq, tseq+trange[0]);
    std::reverse(tseq_left.begin(), tseq_left.end());

    ExtendRight(qseq_left.c_str(), qseq_left.size(), tseq_left.c_str(), tseq_left.size(), 0, 0, result_left);
    std::reverse(result_left.alignment.begin(), result_left.alignment.end());
    result_left.alignment.insert(result_left.alignment.end(), result_right.alignment.begin(), result_right.alignment.end());
    result_left.distance += result_right.distance;
    result_left.qstart = qseq_left.size() - result_left.qend;
    result_left.qend = result_right.qend;
    result_left.tstart = tseq_left.size() - result_left.tend;
    result_left.tend = result_right.tend;


    //printf("Align-new: (%zd %zd) - (%zd %zd) - %d, %f\n",result_left.qstart, result_left.qend, result_left.tstart, result_left.tend, result_left.distance, (double)result_left.distance/ (result_left.tend - result_left.tstart));
    al.target_start = result_left.tstart;
    al.target_end = result_left.tend;
    al.query_start = result_left.qstart;
    al.query_end = result_left.qend;
    al.distance = result_left.distance;
    AlignedString(&result_left.alignment[0], result_left.alignment.size(), qseq+al.query_start, tseq+al.target_start, 
            al.aligned_query, al.aligned_target);
    //LOG(INFO)("X EdlibAligner::AlignBlock");
    return true;
}

bool EdlibAligner::Align1(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t, 2> qrange, std::array<size_t, 2> trange, Alignment& al) {

    size_t tstart = trange[0];
    size_t tend = trange[1];
    size_t qstart = qrange[0];
    size_t qend = qrange[1];
    
    auto r = edlibAlign(qseq+qstart, qend-qstart, tseq+tstart, tend-tstart,
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (r.status == EDLIB_STATUS_OK) {
            assert(r.numLocations >= 1);
            assert(r.startLocations[0] == 0);

            size_t qbend = qsize - 1;
            size_t tbend = r.endLocations[0];   // 不需要 -1
            al.target_start = tstart;
            al.target_end = tend;
            al.query_start = qstart;
            al.query_end = qend;
            al.distance = r.editDistance;
            AlignedString(r.alignment, r.alignmentLength, qseq+al.query_start, tseq+al.target_start, 
                    al.aligned_query, al.aligned_target);
            DEBUG_printf("Distance: %d\n", al.distance);
    }

    return true;
}

void EdlibAligner::AlignedString(const unsigned char* alignment, int alignmentLenght, const char* query, const char* target, std::string& aligned_query, std::string& aligned_target) {
    size_t index_target  = 0;
    size_t index_query = 0; 
    for (int i=0; i<alignmentLenght; ++i) {

        if (alignment[i] == EDLIB_EDOP_MATCH) {
            aligned_target.push_back("ACGT"[target[index_target++]]);
            aligned_query.push_back("ACGT"[query[index_query++]]);

        } else if (alignment[i] == EDLIB_EDOP_INSERT) {
            aligned_target.push_back('-');
            aligned_query.push_back("ACGT"[query[index_query++]]);

        } else if (alignment[i] == EDLIB_EDOP_DELETE) {
            aligned_target.push_back("ACGT"[target[index_target++]]);
            aligned_query.push_back('-');

        } else {
            assert(alignment[i] == EDLIB_EDOP_MISMATCH);
            aligned_target.push_back("ACGT"[target[index_target++]]);
            //aligned_query.push_back('-');
            //aligned_target.push_back('-');
            aligned_query.push_back("ACGT"[query[index_query++]]);
        }
    }
}






// std::array<int, 2> Aligner::ExtendAlignment(const uint8_t *t, int tlen, const uint8_t *q, int qlen) {
    
//     // parameter from minimap2
//     int8_t a = 2;      // score of match
//     int8_t b = -4;     // score of mismatch
//     int8_t ambi = 1;
//     int gapo = 4;
//     int gape = 2;
//     int zdrop = 200;
//     int bw = int(500*1.5+1);
//     int8_t mat[25] = {
//         a, b, b, b, ambi,
//         b, a, b, b, ambi,
//         b, b, a, b, ambi,
//         b, b, b, a, ambi, 
//         ambi, ambi, ambi, ambi, ambi, 
//     };

//     ksw_extz_t ez;
//     memset(&ez, 0, sizeof(ksw_extz_t));
//     //  ksw_extz(0, qs, &qseq[0], ts, &tseq[0], 5, mat, gapo, gape, bw, zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, &ez);
//     ksw_extz(0, qlen, q, tlen, t, 5, mat, gapo, gape, bw, zdrop, KSW_EZ_EXTZ_ONLY, &ez);
//     kfree(0, ez.cigar);
        
//     int qend = ez.reach_end ? qlen : (0 + ez.max_q + 1);
//     int tend = ez.reach_end ? (0 + ez.mqe_t + 1) : (0 + ez.max_t + 1);
//     return {tend, qend};
// }

// void Aligner::ExtendLeft(uint8_t *t, int tlen, uint8_t *q, int qlen, int &tstart, int &qstart) {

//         std::reverse(t, t+tstart);
//         std::reverse(q, q+qstart);
//         auto ends = Aligner::ExtendAlignment(t, tstart, q, qstart);
//         std::reverse(t, t+tstart);
//         std::reverse(q, q+qstart);

//         tstart -= ends[0];
//         qstart -= ends[1];
// }

// void Aligner::ExtendRight(const uint8_t *t, int tlen, const uint8_t *q, int qlen, int &tend, int &qend) {
    
//     assert(tlen > tend && qlen > qend);

//     auto tostr = [](const uint8_t *t, int tlen) {
//         std::string str(tlen, '-');
//         for (int i=0; i<tlen; ++i) {
//             str[i] = "ACGT"[t[i]];
//         }
//         return str;
//     };
//     auto ends = ExtendAlignment(t+tend, tlen-tend, q+qend, qlen-qend);
//     tend += ends[0];
//     qend += ends[1];
// }

// Alignment Aligner::AlignEdlib1(const std::string &query, const std::array<int, 4> &range) {
//     Alignment result(target_, query);
//     int qs = range[0];
//     int qe = range[1];
//     int ql = query.size();
//     int ts = range[2];
//     int te = range[3];
//     int tl = target_.size();


//     std::vector<uint8_t> &tseq = target_0123_;(target_.size(), 0);

//     std::vector<uint8_t> qseq(query.size(), 0);
//     std::transform(query.begin(), query.end(), qseq.begin(), [](uint8_t a){ return DnaSeq::Serial(a);});

//     if (qs > 0 && ts > 0) {
//         ExtendLeft(&tseq[0], tl, &qseq[0], ql, ts, qs);
//     }
    
//     if (ql - qe > 0 && tl - te > 0) {
//         ExtendRight(&tseq[0], tl, &qseq[0], ql, te, qe);
//     } 
     
//     EdlibAlignResult r = edlibAlign(query.c_str()+qs, qe-qs, target_.c_str()+ts, te-ts,
//         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

//     if (r.status == EDLIB_STATUS_OK) {
//         auto ident = ComputeIdentity(r.alignment, r.alignmentLength, local_window_size_);
//         if (ident[0] >= min_global_identity_ && ident[1] >= min_local_identity_ ) {
//             assert(r.numLocations >= 1);
//             result.target_start = ts + r.startLocations[0];
//             result.target_end = ts + r.endLocations[0] + 1;
//             result.query_start = qs;
//             result.query_end = qe;
//             result.distance = r.editDistance;

//             AlignedString(r.alignment, r.alignmentLength, query.c_str()+result.query_start, target_.c_str()+result.target_start, 
//                             result.aligned_query, result.aligned_target);
//             TrimResult(result);
//         }
//     }

//     edlibFreeAlignResult(r);
//     return result;
// }


}