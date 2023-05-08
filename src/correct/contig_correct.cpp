#include "contig_correct.hpp"

#include <edlib.h>
#include <iostream>
#include "./utils/logger.hpp"
#include "utility.hpp"

namespace fsa {


ArgumentParser ContigCorrect::GetArgumentParser() {
    ArgumentParser ap;
    ap.AddPositionOption(overlap_fname_, "ol_fname", "overlap file name");
    ap.AddPositionOption(ctg_fname_, "ctg_fname", "contig file name");
    ap.AddPositionOption(rread_fname_, "rr_fname", "raw read file name");
    ap.AddPositionOption(cread_fname_, "cr_fname", "corrected read file name");

    ap.AddNamedOption(filter1_opts_, "filter1", "overlap filtering options using at loading step", "");
    ap.AddNamedOption(thread_size_, "thread_size", "thread size");
    ap.AddNamedOption(output_directory_, "output_directory", "The directory for temporary files");
    ap.AddNamedOption(read_name_, "read_name", "Set read name for correcting");
    ap.AddNamedOption(read_name_fname_, "read_name_fname", "Set read name for correcting");
    ap.AddNamedOption(aligner_, "aligner", "method for local alignment, diff|edlib.");
    ap.AddNamedOption(score_, "score", "");
    ap.AddNamedOption(coverage_, "coverage", "");
    ap.AddNamedOption(min_length_, "min_length", "");
    ap.AddNamedOption(min_aligned_length_, "min_aligned_length", "");
    ap.AddNamedOption(min_identity_, "min_identity", "");
    ap.AddNamedOption(min_local_identity_, "min_local_identity", "");
    ap.AddNamedOption(min_coverage_, "min_coverage", "");
    ap.AddNamedOption(max_overhang_, "max_overhang", "");
    ap.AddNamedOption(max_overhang_rate_, "max_overhang_rate", "");

    return ap;
}

void ContigCorrect::CheckArguments() {
    filter0_.From(filter0_opts_);
    filter1_.From(filter1_opts_);
    
    filter0_opts_ = filter0_.ToString();
    filter1_opts_ = filter1_.ToString();

}

void ContigCorrect::Running() {
    read_store_.Load(ctg_fname_, "", true);
    LoadReadIds();
    read_store_.Load(rread_fname_, "", true);

    LoadOverlaps(overlap_fname_);

    ol_store_.GroupTarget(groups_, thread_size_);

    LOG(INFO)("Start Correcting");
    Correct();
}

void ContigCorrect::LoadOverlaps(const std::string &fname) {
    std::unordered_set<Seq::Id> ids(read_ids_.begin(), read_ids_.end());

    auto filter = [this, &ids](Overlap &o) {
        bool rel = ids.find(o.b_.id) != ids.end();
        return rel;// && filter0_.ValidQuery(o);

    };
    
    ol_store_.LoadFast(fname, "", (size_t)thread_size_, filter);

    LOG(INFO)("Load %zd overlaps from file %s", ol_store_.Size(), fname.c_str());

}



void ContigCorrect::LoadReadIds() {
    if (!read_name_.empty()) {
        read_ids_.push_back(read_store_.QueryIdByName(read_name_));
    } else if (!read_name_fname_.empty()) {
        std::ifstream file(read_name_fname_);
        std::string line;
        while (std::getline(file, line)) {
            read_ids_.push_back(read_store_.QueryIdByName(line));
        }

    } else {
        std::array<int, 2> range = read_store_.GetIdRange();
        for (int i=range[0]; i<range[1]; ++i) {
            read_ids_.push_back(i);
        }
    }
}


void ContigCorrect::Correct() {
    std::mutex mutex;
    std::ofstream of_cread(cread_fname_);

   auto save_contig = [&](const std::string& name, const std::string &seq) {
        std::lock_guard<std::mutex> lock(mutex);
        of_cread << ">" << name << "\n" << seq << "\n";
    };

    std::vector<std::shared_ptr<WindowJob>> windows;
    std::vector<std::shared_ptr<ContigJob>> jobs;
    for (auto i : read_ids_) {        
        jobs.push_back(std::shared_ptr<ContigJob>(new ContigJob(i, read_store_.GetSeqLength(i), groups_[i], window_size_, overlap_size_)));
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
                LOG(INFO)("Write contig: %s", read_store_.QueryNameById(wjob.owner->tid).c_str());
                save_contig(read_store_.QueryNameById(wjob.owner->tid), wjob.owner->GetSeq());
                //save_contig(read_store_.QueryNameById(wjob.owner->tid), wjob.seq);
                LOG(INFO)("Write contig: %s", read_store_.QueryNameById(wjob.owner->tid).c_str());
                //break;
            }
            if (curr % 100 == 0) {
                LOG(INFO)("Jobs done: %d/%d", curr, windows.size());
            }
            curr = index.fetch_add(1);
        }
    };

 
    LOG(INFO)("thread size %zd, jobsize %d", thread_size_, windows.size());
    if (of_cread.is_open()) {
        MultiThreadRun((size_t)thread_size_, work_func);
    } else {
        LOG(INFO)("Failed to open file: %s", rread_fname_.c_str());
    }

}


bool ContigCorrect::Worker::ExactFilter(const Alignment &r) {
    if (r.AlignSize() < (size_t)owner_.filter1_.min_aligned_length && 
        r.AlignSize() < r.QuerySize() * owner_.filter1_.min_aligned_length) return true;
    
    if (r.Identity() < owner_.min_identity_) return true;

    if (r.AlignSize() >= owner_.filter1_.min_accept_aligned_length) return false;

    const double oh_rate = owner_.filter1_.max_overhang_rate;


    size_t t_overhang = std::max(size_t(r.TargetSize()*oh_rate), (size_t)owner_.filter1_.max_overhang);
    size_t q_overhang = std::max(size_t(r.QuerySize()*oh_rate), (size_t)owner_.filter1_.max_overhang);

    if (r.target_start > t_overhang && r.query_start > q_overhang) return true;
    if (r.target_end + t_overhang < r.TargetSize() && r.query_end + q_overhang < r.QuerySize()) return true;

    return false;
}

bool ContigCorrect::Worker::GetAlignment(Seq::Id id, const Overlap* o, Alignment& al) {
    const auto& tread = o->GetRead(id);
    const auto& qread = o->GetOtherRead(id);

    std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
    return aligner_.Align(owner_.read_store_.GetSeq(qread.id), !o->SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。

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


bool ContigCorrect::Worker::Correct(WindowJob &job) {
 
    auto id = job.GetTId();

    std::vector<const Overlap*> cands = job.GetOverlaps();
    if (cands.size() == 0) {    // if the area is not coveraged by any reads.
        job.seq = *DnaSeq(owner_.read_store_.GetSeq(id), job.start, job.end-job.start).ToString();
        return true;
    }

    // 寻找
    int ctgstart = job.start;
    int ctgend = job.end;
    for (auto o : cands) {
        auto& r = o->GetRead(id);
        if (r.start < ctgstart)  ctgstart = r.start;
        if (r.end > ctgend) ctgend = r.end;
    }
    const DnaSeq target(owner_.read_store_.GetSeq(id), ctgstart, ctgend - ctgstart);
    
    CalculateWeight(id, target, cands, ctgstart, {job.start, job.end});

    // 根据权重排序
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
        auto r = aligner_.Align(owner_.read_store_.GetSeq(qread.id), !o->SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
        if (r && !ExactFilter(al)) {
            al.Rearrange();
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
    printf("al size %zd/%zd\n", aligned_.size(), cands.size());

    // auto range = MostEffectiveCoverage(target.Size(), aligned_, 500, owner_.min_coverage_);

    std::array<size_t, 2> range = {(size_t)(job.start - ctgstart), (size_t)(job.end - ctgstart)};
   
    graph_.Build(target, range, aligned_);
    graph_.Consensus();
    job.seq = graph_.GetSequence();

    job.done = true;

    return true;
}

void ContigCorrect::Worker::CalculateWeight(Seq::Id id,  const DnaSeq& target, const std::vector<const Overlap*> & cands, int offset, const std::array<int,2> &range) {
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

bool ContigCorrect::Worker::IsCoverageEnough(const std::vector<int> &cov) {

    auto CoverageStatus = [](const std::vector<int> &cov, int th) {
        return std::accumulate(cov.begin(), cov.end(), 0, [th](int a, int c) { return a + (c<th ? (th-c) : 0); }) * 1.0 / th / cov.size();
    };
            
    return CoverageStatus(cov, owner_.coverage_) < 0.05;
}

ContigCorrect::ContigJob::ContigJob(Seq::Id id, size_t len, const std::unordered_map<int, std::vector<const Overlap*>>& ols, size_t wsize, size_t osize) 
 : tid(id), tlen(len), overlaps(ols), win_size(wsize), ovl_size(osize) {
    
    assert(win_size > ovl_size);

    size_t curr = 0;
    while (curr < tlen) {
        int start = curr > 0 ? curr - ovl_size : 0;
        int end   = curr + 1.5* win_size < tlen ? curr + win_size : tlen;

        windows.push_back(std::shared_ptr<WindowJob>(new WindowJob(this, start, end)));
        curr = end;
    }
} 

std::string ContigCorrect::ContigJob::GetSeq() const {
    assert(windows.size() > 0);

    std::string seq(windows.front()->seq);

    for (size_t i=1; i < windows.size(); ++i) {
        // 取后一节窗口的overlap的中间二分一的数据，在前一个窗口的overlap中寻找。

        const std::string& next = windows[i]->seq; // alias
        printf("next %zd %d %d\n", next.size(), windows[i]->start, windows[i]->end);

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

        } else {
            LOG(WARNING)("Don't find common part in the overlap of windows %d and %d", i-1, i);
            seq.insert(seq.end(), next.begin()+ovl_size, next.end());
        }
        
        edlibFreeAlignResult(r);

    }

    return seq;
}

} // namespace fsa {
