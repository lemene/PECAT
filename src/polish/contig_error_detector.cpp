#include "contig_error_detector.hpp"

namespace fsa {

ContigErrorDetector::ContigErrorDetector(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds) {
    size_t tlen = ds.seq_store_.GetSeqLength(tid);
    ctg_cov_.assign(tlen, BaseCoverage());
}

void ContigErrorDetector::Detect() {
    ComputeCoverage();
    EvaluateQuality();
    CollectCandidates();

}

void ContigErrorDetector::ComputeCoverage() {
    // short name
    const ReadStore& seq_store = dataset_.seq_store_;
    auto  ol_group = dataset_.grouper_.Get(tid_);
    // 
    
    const DnaSeq& target = seq_store.GetSeq(tid_);
    for (size_t i = 0; i < target.Size(); ++i) {
        ctg_cov_[i].ref = target[i];
    }

    for (size_t i = 0; i < ol_group.Size(); i++) {
        for (size_t j = 0; j < ol_group.Size(i); j++) {
            auto ol = ol_group.Get(i, j);
            const DnaSeq& qseq = seq_store.GetSeq(ol->a_.id);
            const DnaSeq& tseq = seq_store.GetSeq(ol->b_.id);
            assert(ol->b_.strand == 0);
        
            std::vector<uint8_t> tal;   
            tal.reserve(ol->AlignedLength()*2);
            std::vector<uint8_t> qal;   
            qal.reserve(ol->AlignedLength()*2);

            if (!ol->IsProper(500)) {
                if (ol->SameDirect()) {
                    if (ol->a_.start >= 500) {
                        ctg_cov_[ol->b_.start].clips ++;
                    }
                    if (ol->a_.len - ol->a_.end >= 500) {
                        ctg_cov_[ol->b_.end-1].clips ++;
                    }
                } else {
                    if (ol->a_.start >= 500) {
                        ctg_cov_[ol->b_.end-1].clips ++;
                    }
                    if (ol->a_.len - ol->a_.end >= 500) {
                        ctg_cov_[ol->b_.start].clips ++;
                    }

                }
            }
        
            auto get_base = [](const Overlap::Read &r, const DnaSeq& seq, size_t idx) {
                return r.strand == 0 ? seq[r.start+idx] : (3 - seq[r.end - idx - 1]);
            };
        
            size_t qidx = 0;        // not from ol.b_.start;
            size_t tidx = ol->b_.start;
            size_t distance = 0;

            for (const auto &d : ol->detail_) {
                switch (d.type){
                case 'M':
                case '=':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        uint8_t cq = get_base(ol->a_, qseq, qidx+i);
                        uint8_t ct = tseq[tidx+i];
                        qal.push_back(cq+1);
                        tal.push_back(ct+1);
                        ctg_cov_[tidx+i].bases[ct]++;
                        if (cq != ct) {
                            distance++;
                        }
                    }
                    qidx += d.len;
                    tidx += d.len;
                    break;
                case 'D':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        char ct = tseq[tidx+i];
                        qal.push_back(0);
                        tal.push_back(ct+1);
                        ctg_cov_[tidx].bases[4]++;
                    }
                    distance += d.len;
                    tidx += d.len;
                    break; 
                case 'I':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        char cq = get_base(ol->a_, qseq, qidx+i);
                        qal.push_back(cq+1);
                        tal.push_back(0);
                    }
                    ctg_cov_[tidx].bases[5]++;
                    ctg_cov_[tidx].inssize += d.len;
                    qidx += d.len;
                    distance += d.len;
                    break;
                default:
                    LOG(ERROR)("never come here");
                }
            } 
        }

    }

    // for (size_t i = 0; i < ctg_cov_.size(); ++i) {
    //     LOG(INFO)("ctgcov(%d): %s", i, ctg_cov_[i].ToString().c_str());
    // }

}


void ContigErrorDetector::ComputeCoverage1() {
    // short name
    const ReadStore& seq_store = dataset_.seq_store_;
    auto  ol_group = dataset_.grouper_.Get(tid_);
    // 
    
    const DnaSeq& target = seq_store.GetSeq(tid_);
    for (size_t i = 0; i < target.Size(); ++i) {
        ctg_cov_[i].ref = target[i];
    }

    for (size_t i = 0; i < ol_group.Size(); i++) {
        for (size_t j = 0; j < ol_group.Size(i); j++) {
            auto ol = ol_group.Get(i, j);
            const DnaSeq& qseq = seq_store.GetSeq(ol->a_.id);
            const DnaSeq& tseq = seq_store.GetSeq(ol->b_.id);
            assert(ol->b_.strand == 0);
        
            std::vector<uint8_t> tal;   
            tal.reserve(ol->AlignedLength()*2);
            std::vector<uint8_t> qal;   
            qal.reserve(ol->AlignedLength()*2);

            if (!ol->IsProper(500)) {
                if (ol->SameDirect()) {
                    if (ol->a_.start >= 500) {
                        ctg_cov_[ol->b_.start].clips ++;
                    }
                    if (ol->a_.len - ol->a_.end >= 500) {
                        ctg_cov_[ol->b_.end-1].clips ++;
                    }
                } else {
                    if (ol->a_.start >= 500) {
                        ctg_cov_[ol->b_.end-1].clips ++;
                    }
                    if (ol->a_.len - ol->a_.end >= 500) {
                        ctg_cov_[ol->b_.start].clips ++;
                    }

                }
                continue;
            }
        
            auto get_base = [](const Overlap::Read &r, const DnaSeq& seq, size_t idx) {
                return r.strand == 0 ? seq[r.start+idx] : (3 - seq[r.end - idx - 1]);
            };
        
            size_t qidx = 0;        // not from ol.b_.start;
            size_t tidx = ol->b_.start;
            size_t distance = 0;

            for (const auto &d : ol->detail_) {
                switch (d.type){
                case 'M':
                case '=':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        uint8_t cq = get_base(ol->a_, qseq, qidx+i);
                        uint8_t ct = tseq[tidx+i];
                        qal.push_back(cq+1);
                        tal.push_back(ct+1);
                        ctg_cov_[tidx+i].bases[ct]++;
                        if (cq != ct) {
                            distance++;
                        }
                    }
                    qidx += d.len;
                    tidx += d.len;
                    break;
                case 'D':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        char ct = tseq[tidx+i];
                        qal.push_back(0);
                        tal.push_back(ct+1);
                        ctg_cov_[tidx].bases[4]++;
                    }
                    distance += d.len;
                    tidx += d.len;
                    break; 
                case 'I':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        char cq = get_base(ol->a_, qseq, qidx+i);
                        qal.push_back(cq+1);
                        tal.push_back(0);
                    }
                    ctg_cov_[tidx].bases[5]++;
                    ctg_cov_[tidx].inssize += d.len;
                    qidx += d.len;
                    distance += d.len;
                    break;
                default:
                    LOG(ERROR)("never come here");
                }
            } 
        }

    }

    // for (size_t i = 0; i < ctg_cov_.size(); ++i) {
    //     LOG(INFO)("ctgcov(%d): %s", i, ctg_cov_[i].ToString().c_str());
    // }

}

class WinIterator {
public:
    WinIterator(size_t len, size_t win_size, size_t stride) : length_(len), win_size_(win_size), stride_(stride) {
        count_ = (length_ - win_size_ + stride_ - 1) / stride_ + 1;
        LOG(INFO)("len: %zd %zd", len, count_);
    }

    std::array<size_t, 2> Get() const { return {stride_*index_, std::min<size_t>(length_, stride_*(index_+1))}; }
    void Start() {
        index_ = 0;
    }
    void Next() {
        index_ ++;
    }
    bool IsEnd() {
        return index_ > count_;
    }

protected:
    size_t index_ = 0;
    size_t length_;
    size_t win_size_;
    size_t stride_;
    size_t count_;
};

void ContigErrorDetector::CollectCandidates() {
    std::vector<std::array<size_t, 2>> cands;
    WinIterator witr = WinIterator(ctg_cov_.size(), 400, 200);
    for (witr.Start(); !witr.IsEnd(); witr.Next()) {
        auto win = witr.Get();

        std::vector<int> ss(win[1] - win[0], 0);
        std::vector<int> bs(win[1] - win[0], 0);
        std::vector<int> ms(win[1] - win[0], 0);
        size_t clips = 0;

        for (size_t i = win[0]; i < win[1]; ++i) {
            const auto& c = ctg_cov_[i];
            ss[i - win[0]] = std::accumulate(ctg_cov_[i].bases, ctg_cov_[i].bases+5, 0);
            bs[i - win[0]] += ctg_cov_[i].inssize;

            auto m = std::max_element(ctg_cov_[i].bases, ctg_cov_[i].bases+6);

            if ((m-c.bases) != ctg_cov_[i].ref) {
                ms[i - win[0]] = *m;
            }
            clips += c.clips;
        }

        auto min = std::min_element(ss.begin(), ss.end());
        if (*min < 3) {
            cands.push_back(win);
        }

        // if (std::accumulate(ms.begin(), ms.end(), 0) * 1.0 / std::accumulate(bs.begin(), bs.end(), 0) < 0.8) {
        //     cands.push_back(win);
        // }

        // if (clips > 5) {
        //     cands.push_back(win);            
        // }
    }

    // merge
    std::vector<std::array<size_t, 2>> merged;
    if (cands.size() > 0) {
        merged.push_back(cands[0]);

        for (size_t i = 1; i < cands.size(); ++i) {
            if (cands[i][0] - merged.back()[1] <= 5000) {
                assert(merged.back()[1] <= cands[i][1]);
                merged.back()[1] = cands[i][1];
            } else {
                merged.push_back(cands[i]);
            }
        }
    }
    LOG(INFO)("Merged: %zd", merged.size());
    for (auto w : merged) {
        LOG(INFO)("[%zd, %zd]", w[0], w[1]);
    }
}

void ContigErrorDetector::VerifyCandidates(const std::vector<std::array<size_t, 2>> &merged) {
    std::vector<std::array<size_t, 2>> verified;

    for (const auto& win : merged) {
        
    }
}

std::vector<size_t> ContigErrorDetector::GetBigInserts(size_t start, size_t end) {
    const int BIG_INSERT = 1000;
    std::vector<size_t> pos;
    for (size_t i = start; i < end; ++i) {
        double s = ctg_cov_[i].inssize * 1.0 / ctg_cov_[i].bases[5] ;
        if (s > BIG_INSERT) {
            if (pos.size() == 0 || pos.back() + BIG_INSERT < i) {
                pos.push_back(i);
            } 

        }
    }

    for (auto p : pos) {
        LOG(INFO)("BIGINSERT %zd %zd %zd", p, ctg_cov_[p].bases[5], ctg_cov_[p].inssize);
    }

    return pos;
}


void ContigErrorDetector::EvaluateQuality() {
    size_t ref_mis = 0;
    size_t rd_mat = 0;
    size_t rd_total = 0;
    //for (const auto& c : ctg_cov_) {
    for (size_t i = 0; i < ctg_cov_.size(); ++i) {
        const auto& c = ctg_cov_[i];
        auto m = std::max_element(c.bases, c.bases+6);
        //LOG(INFO)("S(%zd): %s | %zd %zd", i, c.ToString().c_str(), m-c.bases, *m);

        if ((m-c.bases) != c.ref) {
            ref_mis ++;
        }

        rd_total += c.bases[0] + c.bases[1] + c.bases[2] + c.bases[3] + c.bases[5];

        switch ((m-c.bases)) {
        case 0:
        case 1:
        case 2:
        case 3:
            rd_mat += *m;
            break;
        case 4:
            // 
            break;
        case 5:
            rd_mat += c.inssize;
        default:
            assert(!"never coming here");
        }
    
    }
    LOG(INFO)("%zd, %zd, %zd, %zd", ref_mis, ctg_cov_.size(), rd_mat, rd_total);
    LOG(INFO)("%.02f, %.02f", ref_mis*1.0/ctg_cov_.size(), rd_mat*1.0 / rd_total);
}

}