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

bool CheckAlign(const char* qseq, size_t qsize, const char* tseq, size_t tsize, size_t qstart, size_t tstart, std::vector<unsigned char> alignment) {
    size_t qi = qstart;
    size_t ti = tstart;
    for (size_t i=0; i < alignment.size(); ++i) {
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

    return true;
}

bool EdlibAligner::ExtendRight(const Seq &q, const Seq &t, size_t qstart, size_t tstart, AlignResult &result) {
    size_t qindex = qstart;
    size_t tindex = tstart;
    result.qstart = qstart;
    result.qend = qstart;
    result.tstart = tstart;
    result.tend = tstart;

    if (q.size == 16994)
    DEBUG_printf("edlib_ext_r: %zd -> %zd, %zd -> %zd\n", qstart, q.size, tstart, t.size);

    size_t qleft = q.size - qindex;
    size_t tleft = t.size - tindex;

    while (qleft >= block_size_*1.5 && tleft >= block_size_*1.5) {
        size_t qbsize = block_size_;
        size_t tbsize = block_size_*1.5;

        if (q.size == 16994) {
            DEBUG_printf("edlib_forward0_s: q = %zd -> %zd, t = %zd -> %zd\n", qindex, qbsize, tindex, tbsize);
        }
        AlignResult bresult;
        if (WrapEdlibPrefix(q, t, {qindex, qindex+qbsize}, {tindex, tindex+tbsize}, bresult)) {
            CheckAndTrimRightEnd(bresult, match_count_, block_error_);

            if (bresult.Valid()) {
                assert(bresult.alstart == 0 && bresult.alend <= bresult.alignment.size());
                result.alignment.insert(result.alignment.end(), bresult.alignment.begin()+bresult.alstart, bresult.alignment.begin()+bresult.alend);
                result.distance += bresult.distance;
                tindex = bresult.tend;
                qindex = bresult.qend;
                if (q.size == 16994) {
                    DEBUG_printf("edlib_forward0_e: q = %zd, t = %zd distance=%d\n", qindex, tindex, bresult.distance);
                }

            } else {
                break;
            }
        } else {
            break;
        }

        qleft = q.size - qindex;
        tleft = t.size - tindex;
    }

    if (qleft < block_size_*1.5 || tleft < block_size_*1.5) {
        if (qleft <= tleft) {
            size_t qbsize = qleft;
            size_t tbsize = std::min<size_t>(qleft*1.5, t.size-tindex);

            
            if (q.size == 16994) {
                DEBUG_printf("edlib_forward1_s: q = %zd -> %zd, t = %zd -> %zd\n", qindex, qbsize, tindex, tbsize);
            }

            AlignResult bresult;
            if (WrapEdlibPrefix(q, t, {qindex, qindex+qbsize}, {tindex, tindex+tbsize}, bresult)) {
                CheckAndTrimRightEnd(bresult, match_count_, block_error_);

                if (bresult.Valid()) {
                    assert(bresult.alstart == 0 && bresult.alend <= bresult.alignment.size());
                    result.alignment.insert(result.alignment.end(), bresult.alignment.begin()+bresult.alstart, bresult.alignment.begin()+bresult.alend);
                    result.distance += bresult.distance;
                    tindex = bresult.tend;
                    qindex = bresult.qend;
                                    
                    if (q.size == 16994) {
                        DEBUG_printf("edlib_forward1_e: q = %zd, t = %zd distance=%d\n", qindex, tindex, bresult.distance);
                    }
                } 
            }
        } else {
            size_t qbsize = std::min<size_t>(tleft*1.5, q.size-qindex);
            size_t tbsize = tleft;

            if (q.size == 16994) {
                DEBUG_printf("edlib_forward2_s: q = %zd -> %zd, t = %zd -> %zd\n", qindex, qbsize, tindex, tbsize);
                DEBUG_printf("edlib_forward2_s:q:");
                for (size_t i = qindex; i < qindex+qbsize; ++i) {
                    DEBUG_printf("%c", "ACGT"[q.seq[i]]);

                }
                DEBUG_printf("\n");
                DEBUG_printf("edlib_forward2_s:t:");
                for (size_t i = tindex; i < tindex+tbsize; ++i) {
                    DEBUG_printf("%c", "ACGT"[t.seq[i]]);

                }
                DEBUG_printf("\n");
            }

            AlignResult bresult;
            if (WrapEdlibPrefix(t, q, {tindex, tindex+tbsize}, {qindex, qindex+qbsize}, bresult)) {
            if (q.size == 16994) {
                DEBUG_printf("edlib_forward2_s: bq %zd %zd, bt %zd %zd\n", bresult.tstart, bresult.tend, bresult.qstart, bresult.qend);
            }
                CheckAndTrimRightEnd(bresult, match_count_, block_error_);

                if (bresult.Valid()) {
                    assert(bresult.alstart == 0 && bresult.alend <= bresult.alignment.size());

                    for (size_t i = bresult.alstart; i < bresult.alend; ++i) {
                        if (bresult.alignment[i] == EDLIB_EDOP_INSERT) {
                            result.alignment.push_back(EDLIB_EDOP_DELETE);
                        } else if (bresult.alignment[i] == EDLIB_EDOP_DELETE) {
                            result.alignment.push_back(EDLIB_EDOP_INSERT);
                        } else {
                            result.alignment.push_back(bresult.alignment[i]);
                        }
                    }

                    result.distance += bresult.distance;
                    qindex = bresult.tend;
                    tindex = bresult.qend;   

                    if (q.size == 16994) {
                        DEBUG_printf("edlib_forward2_e: q = %zd, t = %zd distance=%d\n", qindex, tindex, bresult.distance);
                    }
                } 
            }
        }
    }
    result.qstart = qstart;
    result.qend = qindex;
    result.tstart = tstart;
    result.tend = tindex;
    
    result.alend = result.alignment.size();
    result.alstart = 0;
        
    CheckAlign(q.seq, q.size, t.seq, t.size, result.qstart, result.tstart, result.alignment);
    
    return true;
}


bool EdlibAligner::Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t, 2> qrange, std::array<size_t, 2> trange, Alignment& al) {

    AlignResult result_right;
    ExtendRight({qseq, qsize}, {tseq, tsize}, qrange[0], trange[0], result_right);

    AlignResult result_left;
    std::string qseq_left(qseq, qseq+qrange[0]);
    std::reverse(qseq_left.begin(), qseq_left.end());

    std::string tseq_left(tseq, tseq+trange[0]);
    std::reverse(tseq_left.begin(), tseq_left.end());

    ExtendRight({qseq_left.c_str(), qseq_left.size()}, {tseq_left.c_str(), tseq_left.size()}, 0, 0, result_left);
    assert(result_left.alstart == 0);
    result_left.alignment.erase(result_left.alignment.begin()+result_left.alend, result_left.alignment.end());
    std::reverse(result_left.alignment.begin(), result_left.alignment.end());
    result_left.alignment.insert(result_left.alignment.end(), result_right.alignment.begin()+result_right.alstart, 
         result_right.alignment.begin()+result_right.alend);
    result_left.distance += result_right.distance;
    result_left.qstart = qseq_left.size() - result_left.qend;
    result_left.qend = result_right.qend;
    result_left.tstart = tseq_left.size() - result_left.tend;
    result_left.tend = result_right.tend;

    al.target_start = result_left.tstart;
    al.target_end = result_left.tend;
    al.query_start = result_left.qstart;
    al.query_end = result_left.qend;
    al.distance = result_left.distance;
    AlignedString(&result_left.alignment[0], result_left.alignment.size(), qseq+al.query_start, tseq+al.target_start, 
            al.aligned_query, al.aligned_target);
    return true;
}


void EdlibAligner::CheckAndTrimRightEnd(AlignResult &result, size_t match_count, double block_error) {
    DEBUG_printf("TRIMxxxx %zd >= %zd, %zd >= %zd, %zd >= %zd\n",result.alend ,result.alstart, result.qend , result.qstart, result.tend , result.tstart );
    assert(result.alend >= result.alstart && result.qend >= result.qstart && result.tend >= result.tstart);
    int match = 0; 
    int alend = result.alend - 1;
    int qbend = result.qend - 1;
    int tbend = result.tend - 1;
    int distance = result.distance;
    for (; alend >= (int)result.alstart; --alend) {        
        DEBUG_printf("TRIMxxxx %d, %d\n", alend, result.alignment[alend]);
        if (result.alignment[alend] == EDLIB_EDOP_MATCH) {
            match += 1;
            qbend --;
            tbend --;
        } else if (result.alignment[alend] == EDLIB_EDOP_INSERT) {
            match = 0;
            qbend --;
            distance --; 
        } else if (result.alignment[alend] == EDLIB_EDOP_DELETE) {
            match = 0;
            tbend --;
            distance --;
        } else {
            assert(result.alignment[alend] == EDLIB_EDOP_MISMATCH);
            match = 0;
            qbend --;
            tbend --;
            distance --;
        }

        if ((size_t)match >= match_count && distance*1.0/(tbend-result.tstart + 1 + match) <= block_error) {
            break;
        }
    }
    DEBUG_printf("Trim Check %zd -  %zd - %zd - %zd\n", alend , qbend, tbend, distance);

    if ((size_t)match >= match_count && distance*1.0/(tbend-result.tstart + 1 + match) <= block_error) {
        result.alend = alend + match;
        result.tend = tbend + 1 + match;
        result.qend = qbend + 1 + match;
        result.distance = distance;
        assert(result.alend <= result.alignment.size());

    } else {
        result.alend = result.alstart;
        result.tend = result.tstart;
        result.qend = result.qstart;
        result.distance = 0;
    }
    DEBUG_printf("Trim End %zd - %zd, %zd - %zd, %zd - %zd\n",result.alend ,result.alstart, result.qend , result.qstart, result.tend , result.tstart );

}

bool EdlibAligner::WrapEdlibPrefix(const Seq &q, const Seq &t, std::array<size_t,2> qblock, std::array<size_t, 2> tblock, AlignResult &result) {

    size_t tstart = tblock[0];
    size_t tend = tblock[1];
    size_t qstart = qblock[0];
    size_t qend = qblock[1];
    DEBUG_printf("WrapEdlibPrefix 0: %zd - %zd, %zd - %zd\n", qstart, qend, tstart, tend);

    auto r = edlibAlign(q.seq+qstart, qend-qstart, t.seq+tstart, tend-tstart, 
        edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
            
    if (r.status == EDLIB_STATUS_OK) {
        result.qstart = qstart;
        result.qend = qend;
        result.tstart = tstart;
        result.tend = tstart + r.endLocations[0] + 1;    // +1 符合字符串终点表示习惯
        result.distance = r.editDistance;
        result.alstart = 0;
        result.alend = r.alignmentLength;
        result.alignment.assign(r.alignment, r.alignment+r.alignmentLength);
    }
    edlibFreeAlignResult(r);
    
    if (t.size == 16994)
    DEBUG_printf("WrapEdlibPrefix 1: %d - %d, %d - %d, %d\n", result.qstart, result.qend, result.tstart, result.tend, result.distance);
    return r.status == EDLIB_STATUS_OK;
}

bool EdlibAligner::WrapEdlibGlobal(const Seq &q, const Seq &t, std::array<size_t,2> qblock, std::array<size_t, 2> tblock, AlignResult &result) {

    size_t tstart = tblock[0];
    size_t tend = tblock[1];
    size_t qstart = qblock[0];
    size_t qend = qblock[1];
    
    auto r = edlibAlign(q.seq+qstart, qend-qstart, t.seq+tstart, tend-tstart,
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (r.status == EDLIB_STATUS_OK) {
        assert(r.numLocations >= 1);
        assert(r.startLocations[0] == 0);

        result.qstart = qstart;
        result.qend = qend;
        result.tstart = tstart;
        result.tend = tend;
        result.distance = r.editDistance;
        result.alstart = 0;
        result.alend = r.alignmentLength;
        result.alignment.assign(r.alignment, r.alignment+r.alignmentLength);
    }
    edlibFreeAlignResult(r);
    return r.status == EDLIB_STATUS_OK;
}

void EdlibAligner::AlignedString(const unsigned char* alignment, int alignmentLenght, const char* query, const char* target, std::string& aligned_query, std::string& aligned_target) {
    size_t index_target  = 0;
    size_t index_query = 0; 
    for (int i=0; i<alignmentLenght; ++i) {

        if (alignment[i] == EDLIB_EDOP_MATCH) {
            aligned_target.push_back("ACGT"[(int)target[index_target++]]);
            aligned_query.push_back("ACGT"[(int)query[index_query++]]);

        } else if (alignment[i] == EDLIB_EDOP_INSERT) {
            aligned_target.push_back('-');
            aligned_query.push_back("ACGT"[(int)query[index_query++]]);

        } else if (alignment[i] == EDLIB_EDOP_DELETE) {
            aligned_target.push_back("ACGT"[(int)target[index_target++]]);
            aligned_query.push_back('-');

        } else {
            assert(alignment[i] == EDLIB_EDOP_MISMATCH);
            aligned_target.push_back("ACGT"[(int)target[index_target++]]);
            //aligned_query.push_back('-');     
            //aligned_target.push_back('-');
            aligned_query.push_back("ACGT"[(int)query[index_query++]]);
        }
    }
}

}