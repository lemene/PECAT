#include "ksw2_aligner.hpp"

#include "ksw2.h"
#include "./utils/logger.hpp"
#include "utility.hpp"
#include <cstring>

namespace fsa {

Ksw2Aligner::Ksw2Aligner(const std::vector<std::string> &opts) {
    // edlib:bs=1000

    assert(opts[0] == "ksw2");
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

Ksw2Aligner::~Ksw2Aligner() {

}


void AppendAlignedString(const uint32_t * cigar, size_t cigarLen, const char* query, const char* target, std::string& aligned_query, std::string& aligned_target) {
    size_t index_target  = 0;
    size_t index_query = 0; 
    for (size_t i=0; i<cigarLen; ++i) {
        uint32_t count = cigar[i] >> 4;
        char type = "MID"[cigar[i]&0xF];
        if (type == 'M') {
            for (uint32_t c =0; c <count; ++c) {
                if (target[index_target] == query[index_query]) {
                    aligned_target.push_back("ACGT"[(int)target[index_target++]]);
                    aligned_query.push_back("ACGT"[(int)query[index_query++]]);
                } else {
                    aligned_target.push_back("ACGT"[(int)target[index_target++]]);
                    //aligned_query.push_back('-');
                    //aligned_target.push_back('-');
                    aligned_query.push_back("ACGT"[(int)query[index_query++]]);
                }
            }
        } else if (type == 'I') {
            for (uint32_t c =0; c <count; ++c) {
                aligned_target.push_back('-');
                aligned_query.push_back("ACGT"[(int)query[index_query++]]);
            }

        } else {// if (type == 'D') {
            assert (type == 'D');
            for (uint32_t c =0; c <count; ++c) {
                aligned_target.push_back("ACGT"[(int)target[index_target++]]);
                aligned_query.push_back('-');
            }

        } 
    }
    
}



bool Ksw2Aligner::Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t, 2> qrange, std::array<size_t, 2> trange, Alignment& al) {

    // parameter from minimap2
    int8_t a = 2;      // score of match
    int8_t b = -2;     // score of mismatch
    int8_t ambi = 1;
    int gapo = 4;
    int gape = 2;
    int zdrop = 200;
    int bw = int(1000*1.5+1);
    int8_t mat[25] = {
        a, b, b, b, ambi,
        b, a, b, b, ambi,
        b, b, a, b, ambi,
        b, b, b, a, ambi, 
        ambi, ambi, ambi, ambi, ambi, 
    };

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    size_t tstart = trange[0];
    size_t tend = trange[1];
    size_t qstart = qrange[0];
    size_t qend = qrange[1];


    auto rs = ExtendRight((const uint8_t*)(qseq + qend), qsize - qend, (const uint8_t*)(tseq + tend), tsize - tend);
    DEBUG_printf("ksw2: right %zd, %zd\n", rs[0], rs[1]);

    qend += rs[0];
    tend += rs[1];

    rs = ExtendLeft((const uint8_t*)(qseq), qstart, (const uint8_t*)(tseq), tstart);
    DEBUG_printf("ksw2: left %zd, %zd\n", rs[0], rs[1]);
    qstart -= rs[0];
    tstart -= rs[1];

    const uint8_t* qs = (const uint8_t*)(qseq + qstart);
    const uint8_t* ts = (const uint8_t*)(tseq + tstart);

    //  ksw_extz(0, qs, &qseq[0], ts, &tseq[0], 5, mat, gapo, gape, bw, zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, &ez);
    ksw_extz(0, qend-qstart, &qs[0], tend-tstart, &ts[0], 5, mat, gapo, gape, bw, zdrop, 0, &ez);
        
    qend = ez.reach_end ? qend-qstart : (0 + ez.max_q + 1) ;
    tend = ez.reach_end ? (0 + ez.mqe_t + 1) : (0 + ez.max_t + 1) ;

    qend += qstart;
    tend += tstart;
    
    AppendAlignedString(ez.cigar, ez.n_cigar, qseq+qstart, tseq+tstart, al.aligned_query, al.aligned_target);
    kfree(0, ez.cigar);
    al.target_start = tstart;
    al.target_end = tend;
    al.query_start = qstart;
    al.query_end = qend;
    al.distance = 0;
    for (size_t i=0; i<al.aligned_target.size(); ++i) {
        if (al.aligned_target[i] != al.aligned_query[i]) al.distance += 1;
    }
    return true;
}




std::array<int, 2> Ksw2Aligner::ExtendRight(const uint8_t *q, int qlen, const uint8_t *t, int tlen) {
    
    // parameter from minimap2
    int8_t a = 2;      // score of match
    int8_t b = -4;     // score of mismatch
    int8_t ambi = 1;
    int gapo = 4;
    int gape = 2;
    int zdrop = 200;
    int bw = int(500*1.5+1);
    int8_t mat[25] = {
        a, b, b, b, ambi,
        b, a, b, b, ambi,
        b, b, a, b, ambi,
        b, b, b, a, ambi, 
        ambi, ambi, ambi, ambi, ambi, 
    };

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    //  ksw_extz(0, qs, &qseq[0], ts, &tseq[0], 5, mat, gapo, gape, bw, zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, &ez);
    ksw_extz(0, qlen, q, tlen, t, 5, mat, gapo, gape, bw, zdrop, KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT, &ez);
    kfree(0, ez.cigar);
        
    int qend = ez.reach_end ? qlen : (0 + ez.max_q + 1);
    int tend = ez.reach_end ? (0 + ez.mqe_t + 1) : (0 + ez.max_t + 1);
    return {qend, tend};
}

std::array<int, 2> Ksw2Aligner::ExtendLeft( const uint8_t *q, int qlen, const uint8_t *t, int tlen) {
    
    std::vector<uint8_t> vq(q, q+qlen);
    std::reverse(vq.begin(), vq.end());
    std::vector<uint8_t> vt(t, t+tlen);
    std::reverse(vt.begin(), vt.end());
    auto ends = ExtendRight(&vq[0], vq.size(), &vt[0], vt.size());
    return ends;
}


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