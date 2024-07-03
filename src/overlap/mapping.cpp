#include "mapping.hpp"

#include "edlib.h"
namespace fsa {


void Mapping::BuildIndex() {
    assert(sorted_by_start_.size() == 0 && "Rebuild");
    LOG(INFO)("sort overlaps by start position");
    sorted_by_start_.reserve(ol_store_.Size());
    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        sorted_by_start_.push_back(&ol_store_.Get(i));
    }

    std::sort(sorted_by_start_.begin(), sorted_by_start_.end(), [](const Overlap* a, const Overlap *b) {
        return a->b_.id < b->b_.id || (a->b_.id == b->b_.id && a->b_.start < b->b_.start);
    });

    if (sorted_by_start_.size() > 0) {
        BuildTargetIndex();
        BuildQueryRange();
        BuildQueryIndex();
    }
}

void Mapping::BuildTargetIndex() {
    LOG(INFO)("BuildTargetIndex");
    assert(sorted_by_start_.size() > 0);

    std::pair<Seq::Id, std::array<size_t, 2>> curr = {sorted_by_start_[0]->b_.id, {0, 1}};
    for (size_t i = 0; i < sorted_by_start_.size(); ++i) {
        if (sorted_by_start_[i]->b_.id != curr.first) {
            targets_.insert(curr);
            curr.first = sorted_by_start_[i]->b_.id;
            curr.second[0] = i;
            curr.second[1] = i+1;
        } else {
            curr.second[1] = i+1;
        }
    }
    targets_.insert(curr);
}

void Mapping::BuildQueryRange() {
    LOG(INFO)("BuildQueryRange");
    assert(sorted_by_start_.size() > 0);

    query_ranges_.reserve(sorted_by_start_.size());

    size_t curr = 0;
    for (size_t i = 0; i < sorted_by_start_.size(); ++i) {
        if (sorted_by_start_[curr]->b_.id != sorted_by_start_[i]->b_.id) {
            curr = i;
        }

        if (sorted_by_start_[curr]->b_.end <= sorted_by_start_[i]->b_.start) {
            for (; curr <= i; ++curr) {
                if (sorted_by_start_[curr]->b_.end > sorted_by_start_[i]->b_.start) {
                    break;
                }
            }
        }
        query_ranges_.push_back({i, curr});
    }
    LOG(INFO)("BuildQueryRange1");

    std::sort(query_ranges_.begin(), query_ranges_.end(), [this](const std::array<size_t, 2> &a, const std::array<size_t, 2> &b) {
        const Overlap* aol = sorted_by_start_[a[0]];
        const Overlap* bol = sorted_by_start_[b[0]];
        return (aol->a_.id < bol->a_.id) || 
               (aol->a_.id == bol->a_.id && aol->b_.id < bol->b_.id) ||
               (aol->a_.id == bol->a_.id && aol->b_.id == bol->b_.id  && aol->b_.start < bol->b_.start);
    });
}

void Mapping::BuildQueryIndex() {
    LOG(INFO)("BuildQueryIndex");
    
    assert(query_ranges_.size() > 0);
    std::pair<Seq::Id, std::array<size_t, 2>> curr = {sorted_by_start_[query_ranges_[0][0]]->a_.id, {0, 1}};
    for (size_t i = 0; i < sorted_by_start_.size(); ++i) {
        if (sorted_by_start_[query_ranges_[i][0]]->a_.id != curr.first) {
            queries_.insert(curr);
            curr.first = sorted_by_start_[query_ranges_[i][0]]->a_.id;
            curr.second[0] = i;
            curr.second[1] = i+1;
        } else {
            curr.second[1] = i+1;
        }
    }
    assert(queries_.find(curr.first) == queries_.end());
    queries_.insert(curr);
}

std::unordered_set<Seq::Id> Mapping::GetMappedReads() const {
    std::unordered_set<Seq::Id> mapped;

    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        const auto& ol = ol_store_.Get(i);

        if (ol.AlignedSize() >= 0.9 * ol.a_.len) {
            mapped.insert(ol.a_.id);
        }
    }
    return mapped;
}

std::vector<Mapping::Pair> Mapping::QueryOverlaps(const std::string &name) const {
    auto id = ol_store_.GetStringPool().QueryIdByString(name);
    return id != StringPool::NID ? QueryOverlaps(id) : std::vector<Pair>();
}

std::vector<Mapping::Pair> Mapping::QueryOverlaps(Seq::Id id) const {
    assert (id != StringPool::NID) ;

    std::vector<Mapping::Pair> ols;
    auto tgt = queries_.find(id);
    if (tgt != queries_.end()) {
        
        for (size_t i = tgt->second[0]; i < tgt->second[1]; ++i) {
            auto range = query_ranges_[i];
            const Overlap* tol = sorted_by_start_[range[0]];

            for (size_t ii = range[1]; ii < sorted_by_start_.size(); ++ii) {
                if (ii == range[0]) continue;

                auto qol = sorted_by_start_[ii];
                const int offset = 3000;
                if (qol->b_.end >= tol->b_.start + offset && qol->b_.start + offset <= tol->b_.end ) {
                    ols.push_back({qol, tol});
                }

                if (qol->b_.start > tol->b_.end) {
                    break;
                }
            }
        }

    }

    return ols;
}

void Mapping::Pair::ToOverlap() {
    a_.id = query->a_.id;
    a_.len = query->a_.len;
    a_.strand = query->SameDirect() == target->SameDirect() ? 0 : 1;
    
    b_.id = target->a_.id;
    b_.len = target->a_.len;
    b_.strand = 0;
    identity_ = std::min(query->identity_, target->identity_);

    size_t als_len = std::max<size_t>(query->b_.end, target->b_.end) - 
                     std::min<size_t>(query->b_.start, target->b_.start);
    aligned_query.assign(als_len, -1);
    aligned_target.assign(als_len, -1);

    size_t als_start = std::min<size_t>(query->b_.start, target->b_.start);
    auto align_cigar = [](const Overlap* query, decltype(aligned_query)& aligned, size_t als_start) {
        size_t qcurr = 0;
        size_t tcurr = 0;
        for (const auto& d : query->detail_) {
            switch (d.type) {
            //case 'M':
            case '=':
            //case 'X':
                
                for (size_t i = 0; i < d.len; ++i) {
                    auto p = query->a_.strand == 0 ? (query->a_.start + qcurr + i) : (query->a_.end - qcurr - i-1);
                    
                    aligned[query->b_.start + tcurr + i-als_start] = p;
                }
                qcurr += d.len;
                tcurr += d.len;
                break;
            case 'X':
                qcurr += d.len;
                tcurr += d.len;
                break;
            case 'D':
                // for (size_t i = 0; i < d.len; ++i) {
                //     auto p = query->a_.strand == 0 ? query->a_.start + qcurr : query->a_.end - qcurr ;
                // }
                tcurr += d.len;
                break;

            case 'I':
                qcurr += d.len;
                break;

            default:
                LOG(ERROR)("Not support cigar type '%c'.", d.type);
            }
        }
    };

    align_cigar(query, aligned_query, als_start);
    align_cigar(target, aligned_target, als_start);

    const size_t N = 4;
    start = std::max<size_t>(query->b_.start, target->b_.start) - als_start;
    end = std::min<size_t>(query->b_.end, target->b_.end) - als_start;

    for (; start + N < end; ++start) {
        if (aligned_query[start] == -1 || aligned_target[start]== -1) continue;
        if (std::abs(aligned_query[start] - aligned_query[start+N]) == N &&
            std::abs(aligned_target[start] - aligned_target[start+N]) == N) {      
            break;
        }
    }
    assert(aligned_query[start] != -1 && aligned_target[start]!= -1);

    for (; start + N < end; end--) {
        if (aligned_query[end-1] == -1 || aligned_target[end-1]== -1) continue;
        if (std::abs(aligned_query[end-1] - aligned_query[end-1-N]) == N &&
            std::abs(aligned_target[end-1] - aligned_target[end-1-N]) == N) {
            break;
        }
    }
    assert(aligned_query[end-1] != -1 && aligned_target[end-1]!= -1);
    // if (end - start > 10000) {    
    //     printf("start-end: %zd, %zd\n", start, end);
    //     for (size_t i = start; i < end ; ++i) {
    //         printf("%d xx %d\n", aligned[i][0], aligned[i][1]);
    //     }
    //     fflush(stdout);
    //     assert(0);
    // }
    if (query->a_.strand == 0) {
        a_.start = aligned_query[start];
        a_.end = aligned_query[end-1] + 1;
    } else {
        a_.start = aligned_query[end-1];
        a_.end = aligned_query[start] + 1;
    }
    
    if (target->a_.strand == 0) {
        b_.start = aligned_target[start];
        b_.end = aligned_target[end-1] + 1;
    } else {
        b_.start = aligned_target[end-1];
        b_.end = aligned_target[start] + 1;
    }    

    //printf("%d %d\n", target->SameDirect(), target->SameDirect()); 
    //printf("%d %d %d <-> %d %d %d\n", b_.start , b_.end , b_.len , a_.start , a_.end , a_.len );
    //fflush(stdout);
    assert(0 <= b_.start  && b_.start <= b_.end && b_.end <= b_.len);
    assert(0 <= a_.start  && a_.start <= a_.end && a_.end <= a_.len);
}

std::vector<uint8_t> Mapping::Pair::AlignBases(Seq::Id tid, const DnaSeq &qseq, const DnaSeq &tseq) {
    const DnaSeq* proxy_qseq = &qseq;
    DnaSeq qseq_rv;
    std::vector<int> alt;
    std::vector<int> alq;
    if (tid == b_.id) { // consistent  
        if (query->SameDirect() && target->SameDirect()) {
            alq.assign(aligned_query.begin()+start, aligned_query.begin()+end);
            alt.assign(aligned_target.begin()+start, aligned_target.begin()+end);
        } else if (query->SameDirect() && !target->SameDirect()) {
            //alq.assign(aligned_query.begin()+start, aligned_query.begin()+end);
            alq.assign(aligned_query.rbegin()+aligned_query.size()-end, aligned_query.rbegin()+aligned_query.size()-start);
            for (size_t i = 0; i < alq.size(); ++i) {
                if (alq[i] >= 0) {
                    alq[i] = qseq.Size() - 1 - alq[i];
                }
            }
            qseq_rv = DnaSeq::ReverseComplement(qseq);
            proxy_qseq = &qseq_rv; 
            assert(Seq::ReverseComplement(*(qseq.ToString())) == *(qseq_rv.ToString()));
            alt.assign(aligned_target.rbegin()+aligned_target.size()-end, aligned_target.rbegin()+aligned_target.size()-start);
        } else if (!query->SameDirect() && target->SameDirect()) {
            alq.assign(aligned_query.begin()+start, aligned_query.begin()+end);
            for (size_t i = 0; i < alq.size(); ++i) {
                if (alq[i] >= 0) {
                    alq[i] = qseq.Size() - 1 - alq[i];
                }
            }
            qseq_rv = DnaSeq::ReverseComplement(qseq);
            proxy_qseq = &qseq_rv;
            alt.assign(aligned_target.begin()+start, aligned_target.begin()+end);
        } else {
            // reverse
            alq.assign(aligned_query.rbegin()+aligned_query.size()-end, aligned_query.rbegin()+aligned_query.size()-start);
            alt.assign(aligned_target.rbegin()+aligned_target.size()-end, aligned_target.rbegin()+aligned_target.size()-start);
        }
    } else {
        assert(tid == a_.id); // inconsistent

        if (query->SameDirect() && target->SameDirect()) {
            alq.assign(aligned_target.begin()+start, aligned_target.begin()+end);
            alt.assign(aligned_query.begin()+start, aligned_query.begin()+end);
        } else if (query->SameDirect() && !target->SameDirect()) {
            alq.assign(aligned_target.rbegin()+aligned_target.size()-end, aligned_target.rbegin()+aligned_target.size()-start);
            for (size_t i = 0; i < alq.size(); ++i) {
                if (alq[i] >= 0) {
                    alq[i] = qseq.Size() - 1 - alq[i];
                }
            }
            qseq_rv = DnaSeq::ReverseComplement(qseq);
            proxy_qseq = &qseq_rv;
            alt.assign(aligned_query.rbegin()+aligned_query.size()-end, aligned_query.rbegin()+aligned_query.size()-start);
        } else if (!query->SameDirect() && target->SameDirect()) {
            alq.assign(aligned_target.begin()+start, aligned_target.begin()+end);
            for (size_t i = 0; i < alq.size(); ++i) {
                if (alq[i] >= 0) {
                    alq[i] = qseq.Size() - 1 - alq[i];
                }
            }
            qseq_rv = DnaSeq::ReverseComplement(qseq);
            proxy_qseq = &qseq_rv;
            alt.assign(aligned_query.begin()+start, aligned_query.begin()+end);
        } else {
            // reverse
            alq.assign(aligned_target.rbegin()+aligned_target.size()-end, aligned_target.rbegin()+aligned_target.size()-start);
            alt.assign(aligned_query.rbegin()+aligned_query.size()-end, aligned_query.rbegin()+aligned_query.size()-start);
        }
    }
    // for (size_t i = start; i < end; ++i) {
    //     if (aligned_query[i] != -1 && aligned_target[i] != -1) {
    //         printf("old-al: %zd %d %d - %d %d\n", i, aligned_query[i], aligned_target[i], qseq[aligned_query[i]], tseq[aligned_target[i]]);
    //     } else {
    //         printf("old-al: %zd %d %d \n", i, aligned_query[i], aligned_target[i] );
    //     }
    // }
    // for (size_t i = 0; i < alt.size(); ++i) {
    //     if (alq[i] != -1 && alt[i] != -1) {
    //         printf("new-al: %zd %d %d - %d %d\n", i, alq[i], alt[i], (*proxy_qseq)[alq[i]], tseq[alt[i]]);
    //     } else {
    //         printf("new-al: %zd %d %d\n", i, alq[i], alt[i]);
    //     }
    // }
    //assert(alq.size() == alt.size());
    //printf("direct %d %d %d\n", tid == b_.id, query->SameDirect(), target->SameDirect());
    //assert(0);
    return AlignBases00(*proxy_qseq, alq, tseq, alt);
}

std::vector<uint8_t> Mapping::Pair::AlignBases00(
    const DnaSeq &qseq, const std::vector<int> &al_q_2_ref, 
    const DnaSeq &tseq, const std::vector<int> &al_t_2_ref) {

    std::vector<uint8_t> alignment;
    alignment.reserve(AlignedLength()*2);

    size_t al_ref_start = 0;
    size_t al_ref_end = al_t_2_ref.size();
    assert(al_ref_end > al_ref_start);
    assert(al_q_2_ref[al_ref_start] >= 0);
    assert(al_t_2_ref[al_ref_start] >= 0);

    if (tseq[al_t_2_ref[al_ref_start]] == qseq[al_q_2_ref[al_ref_start]]) {
        alignment.push_back(EDLIB_EDOP_MATCH);
    } else {
        alignment.push_back(EDLIB_EDOP_MISMATCH);
    }

    int status = 0; // 0 MATCH, 1 MISMATCH
    size_t status_start = al_ref_start;
    for (size_t i = al_ref_start+1; i < al_ref_end; ++i) {
        // printf("running: %zd in (%zd %zd) %d %d %zd\n", i, al_ref_start, al_ref_end, al_q_2_ref[i], al_t_2_ref[i], alignment.size());
        if (al_q_2_ref[i] != -1 && al_t_2_ref[i] != -1) { // MATCH
            if (status == 1) {
                // mismatch block: (status_start, i)
                // realign the block
                // printf("mismatch1: %zd %zd\n", status_start, i);
                
                auto ss = RealignBlock(qseq, al_q_2_ref, tseq, al_t_2_ref, status_start, i);
                alignment.insert(alignment.end(), ss.begin(), ss.end());

                status_start = i;
            }

            status = 0;
            if (al_q_2_ref[i] - al_q_2_ref[status_start] == i - status_start &&
                al_t_2_ref[i] - al_t_2_ref[status_start] == i - status_start) {
                // no gap
                alignment.push_back(qseq[al_q_2_ref[i]] == tseq[al_t_2_ref[i]] ?
                    EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH);
            } else {
                // exist gaps
                // match block: [status_start, i)
                // printf("match: %zd %zd\n", status_start, i);

                // mismatch block: (i-1, i)
                // printf("mismatch: %zd %zd\n", i-1, i);
                auto ss = RealignBlock(qseq, al_q_2_ref, tseq, al_t_2_ref, i-1, i);
                assert(ss.size() > 0);
                alignment.insert(alignment.end(), ss.begin(), ss.end());
                
                alignment.push_back(qseq[al_q_2_ref[i]] == tseq[al_t_2_ref[i]] ?
                    EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH);
                status_start = i;
            }

        } else {    // MISMATCH
            if (status == 0) {
                // match block: [status_start, i)
                // printf("match1: %zd %zd\n", status_start, i);
                status = 1;
                status_start = i - 1; // 开区间
            }
        }
    }
    return alignment;
}

std::vector<uint8_t> Mapping::Pair::RealignBlock(
    const DnaSeq &qseq, const std::vector<int> &al_q_2_ref, 
    const DnaSeq &tseq, const std::vector<int> &al_t_2_ref, 
    size_t start, size_t end) {
    // 
    assert(al_q_2_ref[start] >= 0 && al_q_2_ref[end] >= 0 && al_q_2_ref[end] > al_q_2_ref[start]);
    assert(al_t_2_ref[start] >= 0 && al_t_2_ref[end] >= 0 && al_t_2_ref[end] > al_t_2_ref[start]);

    std::vector<uint8_t> ts;
    std::vector<uint8_t> qs;
    for (int ii = al_q_2_ref[start]+1; ii < al_q_2_ref[end]; ++ii) {
        qs.push_back(qseq[ii]);
    }
    for (int ii = al_t_2_ref[start]+1; ii < al_t_2_ref[end]; ++ii) {
        ts.push_back(tseq[ii]);
    }

    if (qs.size() == 0) {
        return std::vector<uint8_t>(ts.size(), EDLIB_EDOP_DELETE);
    }

    if (qs.size() == 1) {
        std::vector<uint8_t> ss;
        auto p = std::find(ts.begin(), ts.end(), qs[0]);
        for (auto i = ts.begin(); i < p; ++i) {
            ss.push_back(EDLIB_EDOP_DELETE);
        }
        if (p < ts.end()) {
            ss.push_back(EDLIB_EDOP_MATCH);
            for (auto i =p+1; i < ts.end(); ++i) {
                ss.push_back(EDLIB_EDOP_DELETE);
            }
        } else {
            ss.push_back(EDLIB_EDOP_INSERT);
        }
        return ss;
    }
    
    

    if (ts.size() == 0) {
        return std::vector<uint8_t>(qs.size(), EDLIB_EDOP_INSERT);
    }

    if (ts.size() == 1) {
        std::vector<uint8_t> ss;
        auto p = std::find(qs.begin(), qs.end(), ts[0]);
        for (auto i = qs.begin(); i < p; ++i) {
            ss.push_back(EDLIB_EDOP_INSERT);
        }
        if (p < qs.end()) {
            ss.push_back(EDLIB_EDOP_MATCH);
            for (auto i =p+1; i < qs.end(); ++i) {
                ss.push_back(EDLIB_EDOP_INSERT);
            }
        } else {
            ss.push_back(EDLIB_EDOP_DELETE);
        }
        return ss;
    }

    std::vector<uint8_t> ss;
    auto r = edlibAlign((const char*)&qs[0], qs.size(), (const char*)&ts[0], ts.size(), 
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
            
    if (r.status == EDLIB_STATUS_OK) {
        ss.assign(r.alignment, r.alignment+r.alignmentLength);
    } else {
        ss.assign(qs.size() + ts.size(), EDLIB_EDOP_INSERT);
        std::fill(ss.begin()+qs.size(), ss.end(), EDLIB_EDOP_DELETE);
    }
    edlibFreeAlignResult(r);
    return ss;
}

} // namespace fsa