#include "aligner.hpp"

#include <cassert>
#include <cstring>

#include <algorithm>
#include <atomic>
#include <iostream>


#include "./utils/logger.hpp"
#include "sequence.hpp"
#include "simple_align.hpp"
#include "utility.hpp"

#include "./align/diff_aligner.hpp"
#include "./align/edlib_aligner.hpp"
#include "./align/ksw2_aligner.hpp"


namespace fsa {

Aligner::Aligner()  {    
}

Aligner::~Aligner() {
}


void Aligner::SetTarget(const DnaSeq& tseq) {
    target_ = &tseq; 
    target_0123_ = target_->ToUInt8();
}

void Aligner::SetTarget(const DnaSeq& tseq, const std::array<size_t, 2> trange) {
    SetTarget(tseq);
    trange_[0] = std::max((size_t)0, trange[0]);
    trange_[1] = std::min(tseq.Size(), trange[1]);

}
void Aligner::SetParameter(const std::string& name, double v) {
    if (name == "min_identity") {
        min_identity_ = v;
    } else if (name == "min_local_identity") {
        min_local_identity_ = v;
    } else {
        LOG(ERROR)("Not support parameter: %s", name.c_str());
    }
}

void Aligner::SetParameter(const std::string &name, const std::string &v) {
    if (name == "aligner") {
        SetAligner(v);
    } else  {
        LOG(ERROR)("Not support parameter: %s", name.c_str());
    }
}

void Aligner::SetAligner(const std::string& opts) {
    auto ss = SplitStringByChar(opts, ':');
    if (ss.size() >= 1) {
        if (ss[0] == "diff") {
            worker.reset(new DiffAligner(ss));
        } else if (ss[0] == "edlib") {
            worker.reset(new EdlibAligner(ss));
        } else if (ss[0] == "ksw2") {
            worker.reset(new Ksw2Aligner(ss));
        } else {
            LOG(ERROR)("Not support parameter: aligner=%s", opts.c_str());
        }
    } else {
        LOG(ERROR)("Not support parameter: aligner=%s", opts.c_str());

    }
}

std::array<double,2> Aligner::ComputeIdentity(const std::string& alq, const std::string& alt, size_t window_size) {
    assert(alq.size() == alt.size() && window_size <= alq.size());

    std::vector<int> score(alq.size(), 0);
    for (size_t i=0; i < alq.size(); ++i) {
        if (alq[i] != alt[i]) {
            score[i] = 1;
        }
    }
    
    std::vector<int> local_identity(alq.size() - window_size + 1, 0);
            
    local_identity[0] = std::accumulate(score.begin(), score.begin()+window_size, 0);
    for (size_t i=1; i<local_identity.size(); ++i) {
        local_identity[i] = local_identity[i-1] - score[i-1] + score[i+window_size-1];
    }

    return {100.0 - std::accumulate(score.begin(), score.end(), 0)*1.0 / alq.size() * 100, 
        100 - *std::max_element(local_identity.begin(), local_identity.end()) * 1.0 / window_size * 100};
}

bool Aligner::CheckAlignedString(const std::string &q, const std::string &t) {
    if (q.size() != t.size()) return false;

    for (size_t i=0; i< q.size(); i++) {
        if (q[i] != t[i] && q[i] != '-' && t[i] != '-') {
            printf("%zd, %c, %c\n", i, q[i], t[i]);
            return false;
        }
    }

    return true;
}


void Aligner::AppendAlignedString(const uint32_t * cigar, size_t cigarLen, const char* query, const char* target, std::string& aligned_query, std::string& aligned_target) {
    size_t index_target  = 0;
    size_t index_query = 0; 
    for (size_t i=0; i<cigarLen; ++i) {
        uint32_t count = cigar[i] >> 4;
        char type = "MID"[cigar[i]&0xF];

        if (type == 'M') {
            for (uint32_t c =0; c <count; ++c) {
                if (target[index_target] == query[index_query]) {
                    aligned_target.push_back(target[index_target++]);
                    aligned_query.push_back(query[index_query++]);
                } else {
                    aligned_target.push_back(target[index_target++]);
                    aligned_query.push_back('-');
                    aligned_target.push_back('-');
                    aligned_query.push_back(query[index_query++]);
                }
            }
        } else if (type == 'I') {
            for (uint32_t c =0; c <count; ++c) {
                aligned_target.push_back('-');
                aligned_query.push_back(query[index_query++]);
            }

        } else {// if (type == 'D') {
            assert (type == 'D');
            for (uint32_t c =0; c <count; ++c) {
                aligned_target.push_back(target[index_target++]);
                aligned_query.push_back('-');
            }

        } 
    }
    
}





 bool Aligner::Align(const DnaSeq& query, bool rc, const std::array<int, 4> &range, Alignment &al) {
    al.Reset(target_, &query);
    std::vector<uint8_t>& tseq = target_0123_;
    std::vector<uint8_t> qseq = query.ToUInt8(0, -1, rc);

    int ts = range[2];
    int qs = !rc ? range[0] : qseq.size() - range[1];
    int te = range[3];
    int qe = !rc ? range[1] : qseq.size() - range[0];

    // TODO 
    auto new_s = FindExactMatch(tseq, qseq, {ts, qs});
    DEBUG_printf("pos0: %d %d %d %d\n", qs, qe, ts,te);
    ts = new_s[0];
    qs = new_s[1];
    DEBUG_printf("pos1: %d %d %d %d\n", qs, qe, ts,te);


    auto r = worker->Align((const char*)&qseq[0], qseq.size(),
                       (const char*)&tseq[0], tseq.size(), {(size_t)qs, (size_t)qe}, {(size_t)ts, (size_t)te}, al); 
    if (r) {
        
        bool valid = false;
        DEBUG_printf("q:%s\nt:%s\n", query.ToString()->c_str(), target_->ToString()->c_str());
        DEBUG_printf("alq: %s\nalt: %s\n", al.aligned_query.c_str(), al.aligned_target.c_str());
        DEBUG_printf("global idents: %f > %f\n", al.Identity(), min_identity_);
        
        if (al.Identity() >= min_identity_ ) {
            valid = al.TrimEnds();
            if (valid && min_local_identity_ > 0 && al.aligned_query.size() >= 1.5*local_window_size_) {
                auto idents = ComputeIdentity(al.aligned_query,  al.aligned_target, local_window_size_);
                valid = idents[1] >= min_local_identity_;
                
                DEBUG_printf("local idents: %f > %f\n", idents[1], min_local_identity_);
            }
    
        }
        if (!valid) {
            al.target_start = 0;
            al.target_end = 0;
            al.query_start = 0;
            al.query_end = 0;
        } else {
            al.ComputeLocalDistance(local_window_size_);
        }
    }
    
    return al.Valid(); 
}

std::array<int, 2> Aligner::FindExactMatch(const std::vector<uint8_t>& tseq, const std::vector<uint8_t>& qseq, const std::array<int , 2> &s) {
    const int k = 4;
    const int w = 100;

    std::vector<uint64_t> tks(w, 0);
    std::vector<uint64_t> qks(w, 0);

    assert((int)tseq.size() > s[0] + w && (int)qseq.size() > s[1] + w);

    auto calc_kmer = [](const uint8_t* s, int k) {
        uint64_t kmer = 0;
        for (int i=0; i<k; ++i) {
            kmer = (kmer << 2) + s[i];
        }
        return kmer;
    };
    
    auto move_kmer = [](uint64_t kmer, int k, uint8_t b) {
        uint64_t mask = ~(3L << (2*(k-1)));

        kmer = kmer & mask;
        kmer = (kmer << 2)  + b;
        return kmer;
    };

    tks[0] = calc_kmer(&tseq[s[0]], k);
    for (int i=1; i<w; ++i) {
        tks[i] = move_kmer(tks[i-1], k, tseq[s[0]+k+i-1]);
    }
    
    qks[0] = calc_kmer(&qseq[s[1]], k);
    for (int i=1; i<w; ++i) {
        qks[i] = move_kmer(qks[i-1], k, qseq[s[1]+k+i-1]);
    }

    for (int it = 0; it < w; ++it) {
        for (int iq = 0; iq < w; ++iq) {
            if (tks[it] == qks[iq]) {
                return { s[0]+it, s[1]+iq};
            }
        }
    }

    return s;
}




} // namespace fsa {
