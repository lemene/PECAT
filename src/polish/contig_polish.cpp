#include "contig_polish.hpp"

#include <edlib.h>
#include <iostream>
#include "./utils/logger.hpp"
#include "utility.hpp"

namespace fsa {


ArgumentParser ContigPolish::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}

void ContigPolish::Running() {
    dataset_.Load();

    LOG(INFO)("Start polishing");
    Correct();
}




void ContigPolish::Correct() {
    std::mutex mutex;
    std::ofstream of_cread(opts_.cread_fname_);

   auto save_contig = [&](const std::string& name, const std::string &seq) {
        std::lock_guard<std::mutex> lock(mutex);
        of_cread << ">" << name << "\n" << seq << "\n";
    };

    std::vector<std::shared_ptr<WindowJob>> windows;
    std::vector<std::shared_ptr<ContigJob>> jobs;
    for (auto i : dataset_.read_ids_) {        
        jobs.push_back(std::shared_ptr<ContigJob>(new ContigJob(i, dataset_, opts_.window_size_, opts_.overlap_size_)));
        for (auto &s : jobs.back()->windows) {
            windows.push_back(s);
        }
    }

    std::atomic<size_t> index {0};

    auto work_func = [&](size_t i) {
        Worker worker(*this);
        auto curr = index.fetch_add(1);
        while (curr < windows.size()) {
            
            //if (curr != 0) { curr = index.fetch_add(1);continue; }
            auto& wjob = *windows[curr];
            worker.Correct(wjob);
            worker.Clear();
            //if (true || wjob.owner->IsDone()) {
            if ( wjob.owner->Savable()) {
                LOG(INFO)("Write contig: %s", dataset_.read_store_.QueryNameById(wjob.owner->tid).c_str());
                save_contig(dataset_.read_store_.QueryNameById(wjob.owner->tid), wjob.owner->GetSeq());
                LOG(INFO)("Write contig: %s", dataset_.read_store_.QueryNameById(wjob.owner->tid).c_str());
            }
            if (curr % 100 == 0) {
                LOG(INFO)("Jobs done: %d/%d", curr, windows.size());
            }
            curr = index.fetch_add(1);
        }
    };

 
    LOG(INFO)("thread size %zd, jobsize %d", opts_.thread_size, windows.size());
    if (of_cread.is_open()) {
        MultiThreadRun((size_t)opts_.thread_size, work_func);
    } else {
        LOG(INFO)("Failed to open file: %s", opts_.rread_fname_.c_str());
    }

}


bool ContigPolish::Worker::ExactFilter(const Alignment &r) {
    if (r.AlignSize() < (size_t)owner_.opts_.filter1_.min_aligned_length && 
        r.AlignSize() < r.QuerySize() * owner_.opts_.filter1_.min_aligned_length) return true;
    
    if (r.Identity() < owner_.opts_.min_identity_) return true;

    if (r.AlignSize() >= owner_.opts_.filter1_.min_accept_aligned_length) return false;

    const double oh_rate = owner_.opts_.filter1_.max_overhang_rate;


    size_t t_overhang = std::max(size_t(r.TargetSize()*oh_rate), (size_t)owner_.opts_.filter1_.max_overhang);
    size_t q_overhang = std::max(size_t(r.QuerySize()*oh_rate), (size_t)owner_.opts_.filter1_.max_overhang);

    if (r.target_start > t_overhang && r.query_start > q_overhang) return true;
    if (r.target_end + t_overhang < r.TargetSize() && r.query_end + q_overhang < r.QuerySize()) return true;

    return false;
}

bool ContigPolish::Worker::GetAlignment(Seq::Id tid, const Overlap& ol, Alignment& al, int ctgstart) {
    const auto& tread = ol.GetRead(tid);
    const auto& qread = ol.GetOtherRead(tid);

  
    if (ol.detail_.size() != 0) {
        GetAlignmentFromCigar(tid, ol, al);
        al.target_start -= ctgstart;
        al.target_end -= ctgstart;
        DEBUG_printf("T:%s\nQ:%s\n", al.aligned_target.c_str(), al.aligned_query.c_str());
        return true;
    } else {

        std::array<int, 4> range = {qread.start, qread.end, tread.start-ctgstart, tread.end-ctgstart};
        return aligner_.Align(owner_.dataset_.read_store_.GetSeq(qread.id), !ol.SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
    }


}



void ContigPolish::Worker::GetAlignmentFromCigar(Seq::Id tid, const Overlap& ol, Alignment &al) {
    assert(ol.detail_.size() > 0);
    const DnaSeq& qseq = owner_.dataset_.read_store_.GetSeq(ol.a_.id);
    const DnaSeq& tseq = owner_.dataset_.read_store_.GetSeq(ol.b_.id);
    assert(ol.b_.strand == 0);

    std::vector<uint8_t> tal;   
    tal.reserve(ol.AlignedLength()*2);
    std::vector<uint8_t> qal;   
    qal.reserve(ol.AlignedLength()*2);

    auto get_base = [](const Overlap::Read &r, const DnaSeq& seq, size_t idx) {
        return r.strand == 0 ? seq[r.start+idx] : (3 - seq[r.end - idx - 1]);
    };

    size_t qidx = 0;        // not from ol.b_.start;
    size_t tidx = ol.b_.start;
    size_t distance = 0;
    for (const auto &d : ol.detail_) {
        switch (d.type){
        case 'M':
        case '=':
            for (size_t i = 0; i < (size_t)d.len; ++i) {
                uint8_t cq = get_base(ol.a_, qseq, qidx+i);
                uint8_t ct = tseq[tidx+i];
                qal.push_back(cq+1);
                tal.push_back(ct+1);
                if (cq != ct) {
                    distance++;
                }
            }
            qidx += d.len;
            tidx += d.len;
            break;
        case 'I':
            for (size_t i = 0; i < (size_t)d.len; ++i) {
                char cq = get_base(ol.a_, qseq, qidx+i);
                qal.push_back(cq+1);
                tal.push_back(0);
            }
            qidx += d.len;
            distance += d.len;
            break;
        case 'D':
            for (size_t i = 0; i < (size_t)d.len; ++i) {
                char ct = tseq[tidx+i];
                qal.push_back(0);
                tal.push_back(ct+1);
            }
            distance += d.len;
            tidx += d.len;
            break; 
        default:
            LOG(ERROR)("never come here");
        }
    }
    
    const auto& tread = ol.GetRead(tid);
    const auto& qread = ol.GetOtherRead(tid);

    al.query = &owner_.dataset_.read_store_.GetSeq(qread.id);
    al.target = &owner_.dataset_.read_store_.GetSeq(tread.id);
    assert(al.target!= nullptr);
    al.target_start = tread.start;
    al.target_end = tread.end;
    al.query_start = qread.start;
    al.query_end = qread.end;
    al.distance = distance;

    const char* ACGT = "-ACGT-";
    if (tread.id == ol.b_.id) {
        assert(tread.strand == 0);
        for (size_t i = 0; i < tal.size(); ++i) {
            al.aligned_target.push_back(ACGT[tal[i]]);
            al.aligned_query.push_back(ACGT[qal[i]]);
        }
    } else {
        if (tread.strand == 0) {
            for (size_t i = 0; i < tal.size(); ++i) {
                al.aligned_target.push_back(ACGT[qal[i]]);
                al.aligned_query.push_back(ACGT[tal[i]]);
            }
        } else {
            for (size_t i = 0; i < tal.size(); ++i) {
                al.aligned_target.push_back(ACGT[5 - qal[tal.size()-i-1]]);
                al.aligned_query.push_back(ACGT[5 - tal[tal.size()-i-1]]);
            }
        }
    }
    // for (size_t i = 0, it = 0; i < al.aligned_target.size(); ++i) {
    //     if (al.aligned_target[i] != '-') {
    //         assert(al.aligned_target[i] == "ACGT"[(*al.target)[it+al.target_start]]);
    //         it ++;
    //     } 
    // }

}
      
 std::array<size_t,2> MostEffectiveCoverage(size_t tsize, const std::vector<Alignment> &aligns, size_t stub, int min_coverage) {
    if (aligns.size() == 0) return {0, 0};

    std::vector<int> coverage(tsize+1, 0);
    for (const auto& al : aligns) {
        //assert(al.target_end - al.target_start > 2*stub);
        if (al.target_end - al.target_start > 2*stub) {
            coverage[al.target_start+stub] += 1;
            coverage[al.target_end - stub] -= 1;
        }
     }

     for (size_t i=1; i< coverage.size(); ++i) {
         coverage[i] += coverage[i-1];
     }

    std::vector<std::array<size_t, 2>> ranges;
    int start = -1;
    for (size_t i=0; i<coverage.size(); i++) {
        if (start >= 0) {
            if (coverage[i] < min_coverage) {
                ranges.push_back({(size_t)start, i});
                start = -1;
            }
        } else {
            if (coverage[i] >= min_coverage) {
                start = i;
            }
        }
    }
    if (start > 0) {
        ranges.push_back({(size_t)start, coverage.size()});
    }
    if (ranges.size() > 0) {
        std::sort(ranges.begin(), ranges.end(), [](const std::array<size_t,2> &a, const std::array<size_t,2> &b) {
            return a[1] - a[0] > b[1] - b[0];
        });
        return ranges[0];
    } else {
        
        return {0, 0}; 
    }

}


bool ContigPolish::Worker::Correct(WindowJob &job) {
 
    auto id = job.GetTId();
    std::vector<const Overlap*> cands = job.GetOverlaps();
    if (cands.size() == 0) {    // if the area is not coveraged by any reads.
        job.seqs.push_back(*DnaSeq(owner_.dataset_.read_store_.GetSeq(id), job.start, job.end-job.start).ToString());
        job.ranges.push_back({job.start, job.end});
        return true;
    }
    if (job.start != 449500) return false;
    // 寻找
    int ctgstart = job.start;
    int ctgend = job.end;
    job.owner->ctg_err_dt.GetBigInserts(ctgstart, ctgend);

    LOG(INFO)("start correct: %d - %d", ctgstart, ctgend);
    for (auto o : cands) {
        auto& r = o->GetRead(id);
        if (r.start < ctgstart)  ctgstart = r.start;
        if (r.end > ctgend) ctgend = r.end;
    }
    const DnaSeq target(owner_.dataset_.read_store_.GetSeq(id), ctgstart, ctgend - ctgstart);
    
    CalculateWeight(id, target, cands, ctgstart, {job.start, job.end});

    std::make_heap(cands.begin(), cands.end(), [](const Overlap* a, const Overlap* b) {
       return a->attached < b->attached;    // CAUTION
    });
    
    size_t heap_size = cands.size();
    aligner_.SetTarget(target); 
    std::vector<int> coverage(target.Size(), 0);
    Alignment al;
    while (heap_size > 0) {
        auto o = cands[0];
        al.Reset();
        const auto& tread = o->GetRead(id);
        const auto& qread = o->GetOtherRead(id);

        std::array<int, 4> range = {qread.start, qread.end, tread.start-ctgstart, tread.end-ctgstart};
        assert(range[2] >= 0);
        auto r = GetAlignment(id, *o, al, ctgstart);
        DEBUG_printf("alignment(%s<->%s): (%d, %d, %d) - (%d, %d, %d) %f\n", owner_.dataset_.QueryStringById(o->a_.id).c_str(), owner_.dataset_.QueryStringById(o->b_.id).c_str(),
            al.query_start, al.query_end, al.QuerySize(), al.target_start, al.target_end, al.TargetSize(), al.Identity());
        if (r && !ExactFilter(al) && al.Identity() >= owner_.opts_.min_identity_) {
            //al.Rearrange();
            aligned_.push_back(al);

            std::for_each(coverage.begin()+al.target_start, coverage.begin()+al.target_end, [](int& c) {c++;} );
            if (IsCoverageEnough(coverage) ) {
                break;
            }
        }
        
        std::pop_heap(cands.begin(), cands.begin()+heap_size, [](const Overlap* a, const Overlap* b) {
            return a->attached < b->attached;    // CAUTION
        });
        heap_size--;
    }


    //LOG(INFO)("al size %zd/%zd\n", aligned_.size(), cands.size());

    //auto range = MostEffectiveCoverage(target.Size(), aligned_, 500, owner_.opts_.min_coverage_);
    std::array<size_t,2> range = {job.start - ctgstart, job.end - ctgstart};
   
    graph_.Build(target, range, aligned_);
    graph_.Consensus();
    job.seqs = graph_.GetSequence();
    job.ranges = graph_.GetSequenceRange();
    for (size_t i = 0; i < job.seqs.size(); ++i) {
        auto seq = job.seqs[i];
        auto s = job.ranges[i];
        LOG(INFO)("rrr(%zd/%zd): %zd %zd, %zd, %zd -> %zd(%d-%d) %zd-%zd %d-%d",i, job.seqs.size(),aligned_.size(), s[0], s[1], target.Size(), seq.size(), job.start, job.end,range[0], range[1], ctgstart, ctgend);
    }

    job.done = true;

    return true;
}

void ContigPolish::Worker::CalculateWeight(Seq::Id id,  const DnaSeq& target, const std::vector<const Overlap*> & cands, int offset, const std::array<int,2> &range) {
    std::vector<double> cand_cov_wts (target.Size()+1);

    for (auto o : cands) {
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        const double w = 0.0;
        auto mr = o->MappingTo<2>(t, {0, q.end});
        if (mr[0] < mr[1]) {
            assert(mr[0] <= t.start && mr[1] >= t.end);
            cand_cov_wts[std::max(range[0], mr[0]) - offset] += w;
            cand_cov_wts[std::min(range[1], mr[1]) - offset] -= w;
        } else {
            assert(mr[1] <= t.start && mr[0] >= t.end);
            cand_cov_wts[std::max(range[0], mr[1]) - offset] += w;
            cand_cov_wts[std::min(range[1], mr[0]) - offset] -= w;
        }
        cand_cov_wts[std::max(range[0],t.start) - offset] += 1-w;
        cand_cov_wts[std::min(range[1],t.end) - offset] -= 1-w;
        //cand_cov_wts[t.start] ++;
        //cand_cov_wts[t.end] --;
    }
    for (size_t i=1; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] += cand_cov_wts[i-1];
    }
    
    //for (auto c : cand_cov_wts) { printf("c: %d\n", int(c)); }
    //printf("\n");

    for (size_t i=0; i<cand_cov_wts.size(); ++i) {
        //cand_cov_wts[i] = (cand_cov_wts[i] >= owner_.coverage_ ? 0 : owner_.coverage_ - cand_cov_wts[i]) + 1;
        cand_cov_wts[i] = cand_cov_wts[i] > 0 ? 1/cand_cov_wts[i] : 0;
    }

    std::for_each(cands.begin(), cands.end(), [id, &cand_cov_wts, &range, offset](const Overlap* o) {
        auto &t = o->GetRead(id);
        o->attached = 1000*o->identity_*std::accumulate(cand_cov_wts.begin()+std::max(t.start, range[0])-offset, cand_cov_wts.begin()+std::min(t.end, range[1])-offset, 0.0);
    });
}

bool ContigPolish::Worker::IsCoverageEnough(const std::vector<int> &cov) {

    auto CoverageStatus = [](const std::vector<int> &cov, int th) {
        return std::accumulate(cov.begin(), cov.end(), 0, [th](int a, int c) { return a + (c<th ? (th-c) : 0); }) * 1.0 / th / cov.size();
    };
            
    return CoverageStatus(cov, owner_.opts_.coverage_) < 0.05;
}

ContigPolish::ContigJob::ContigJob(Seq::Id id, const PolDataset& ds, size_t wsize, size_t osize) 
 : tid(id), tlen(ds.read_store_.GetSeqLength(id)), dataset(ds), overlaps(ds.groups_.find(id)->second)
 , win_size(wsize), ovl_size(osize), ctg_err_dt(tid, ds) {
    assert(win_size > ovl_size);

    ctg_err_dt.Detect();

    size_t curr = 0;
    while (curr < tlen) {
        int start = curr > 0 ? curr - ovl_size : 0;
        int end   = curr + 1.5* win_size < tlen ? curr + win_size : tlen;

        windows.push_back(std::shared_ptr<WindowJob>(new WindowJob(this, start, end)));
        curr = end;
    }
} 

std::string ContigPolish::ContigJob::GetSeq() const {
    assert(windows.size() > 0);

    std::string seq(windows.front()->GetSeq());

    for (size_t i=1; i < windows.size(); ++i) {
        // 取后一节窗口的overlap的中间二分一的数据，在前一个窗口的overlap中寻找。

        const std::string next = windows[i]->GetSeq();

        EdlibAlignResult r = edlibAlign(next.c_str(), ovl_size, seq.c_str()+seq.size()-ovl_size, ovl_size, 
            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

        if (r.status == EDLIB_STATUS_OK) {
            assert(r.numLocations >= 1);

            int offq = 0; int offt = 0; int count = 0;
            for (int j = 0; j < r.alignmentLength; ++j) {
                if (r.alignment[j] == EDLIB_EDOP_MATCH) {
                    offq ++; offt ++;
                    count ++;

                } else if (r.alignment[j] == EDLIB_EDOP_INSERT) {
                    offq ++;
                    count = 0;
                } else if (r.alignment[j] == EDLIB_EDOP_DELETE) {
                    offt ++;
                    count = 0;
                } else {
                    assert(r.alignment[j] == EDLIB_EDOP_MISMATCH);
                    offq ++; offt++;
                    count = 0;
                }

                if (count >= 4) break;
            }

            if (count < 4) {
                LOG(WARNING)("Don't find common part in the overlap of windows %d and %d", i-1, i);
            }
            seq.erase(seq.begin() + seq.size() - ovl_size + r.startLocations[0] + offt, seq.end());
            seq.insert(seq.end(), next.begin()+offq, next.end());
            printf("ctg seq: %zd\n", seq.size());

        } else {
            LOG(WARNING)("Don't find common part in the overlap of windows %d and %d", i-1, i);
            seq.insert(seq.end(), next.begin()+ovl_size, next.end());
        }
        
        edlibFreeAlignResult(r);
    }

    return seq;
}

} // namespace fsa {
