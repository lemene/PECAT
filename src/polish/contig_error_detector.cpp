#include "contig_error_detector.hpp"

namespace fsa {

ContigErrorDetector::ContigErrorDetector(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds) {
    size_t tlen = ds.read_store_.GetSeqLength(tid);
    ctg_cov_.assign(tlen, BaseCoverage());
}

void ContigErrorDetector::Detect() {
    ComputeCoverage();

}

void ContigErrorDetector::ComputeCoverage() {
    // short name
    const ReadStore& rd_store = dataset_.read_store_;
    const std::unordered_map<int, std::vector<const Overlap*>>& overlaps = dataset_.groups_.find(tid_)->second;
    // 
    
    //const DnaSeq& target = owner_.dataset_.read_store_.GetSeq(tid);
    for (auto it : overlaps) {
        //const Overlap* ol = it.second;
        for (auto ol : it.second) {
            const DnaSeq& qseq = rd_store.GetSeq(ol->a_.id);
            const DnaSeq& tseq = rd_store.GetSeq(ol->b_.id);
            assert(ol->b_.strand == 0);
        
            std::vector<uint8_t> tal;   
            tal.reserve(ol->AlignedLength()*2);
            std::vector<uint8_t> qal;   
            qal.reserve(ol->AlignedLength()*2);
        
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
                        ctg_cov_[tidx].bases[ct]++;
                        if (cq != ct) {
                            distance++;
                        }
                    }
                    qidx += d.len;
                    tidx += d.len;
                    break;
                case 'I':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        char cq = get_base(ol->a_, qseq, qidx+i);
                        qal.push_back(cq+1);
                        tal.push_back(0);
                    }
                    ctg_cov_[tidx].ins++;
                    ctg_cov_[tidx].inssize += d.len;
                    qidx += d.len;
                    distance += d.len;
                    break;
                case 'D':
                    for (size_t i = 0; i < (size_t)d.len; ++i) {
                        char ct = tseq[tidx+i];
                        qal.push_back(0);
                        tal.push_back(ct+1);
                        ctg_cov_[tidx].del++;
                    }
                    distance += d.len;
                    tidx += d.len;
                    break; 
                default:
                    LOG(ERROR)("never come here");
                }
            } 
        }

    }

    for (size_t i = 0; i < ctg_cov_.size(); ++i) {
        //LOG(INFO)("ctgcov(%d): %s", i, ctg_cov_[i].ToString().c_str());
    }

}


std::vector<size_t> ContigErrorDetector::GetBigInserts(size_t start, size_t end) {
    const int BIG_INSERT = 1000;
    std::vector<size_t> pos;
    for (size_t i = start; i < end; ++i) {
        double s = ctg_cov_[i].inssize * 1.0 / ctg_cov_[i].ins ;
        if (s > BIG_INSERT) {
            pos.push_back(i);

        }
    }
    for (auto p : pos) {
        LOG(INFO)("BIGINSERT %zd %zd %zd", p, ctg_cov_[p].ins, ctg_cov_[p].inssize);
    }
}
}