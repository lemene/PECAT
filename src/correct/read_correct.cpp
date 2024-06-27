#include "read_correct.hpp"

#include <iostream>
#include <atomic>
#include <ctime>
#include "./utils/logger.hpp"
#include "../utility.hpp"

namespace fsa {

TimeCounter tc_align("get_al");
TimeCounter tc_graph("graph");

ReadCorrect::ReadCorrect() {
}

ArgumentParser ReadCorrect::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}

void ReadCorrect::Running() {
    dataset_.Load();
    Correct();
}

void ReadCorrect::Correct() {
    LOG(INFO)("Start Correcting");

    std::mutex mutex;
    
    std::ofstream of_cread(opts_.cread_fname_);
    std::ofstream of_infos(opts_.infos_fname_);
    const size_t flush_block = 20*1024*1024;
    StatInfo stat_info;

    Progress progress(5000, dataset_.read_ids_.size());

    auto combine_func = [&](std::ostringstream &oss_cread, std::ostringstream &oss_scores, StatInfo &si) {
        std::lock_guard<std::mutex> lock(mutex);
        of_cread << oss_cread.str();              
        oss_cread.str("");

        if (of_infos.is_open()) {
            of_infos << oss_scores.str();
        }
        oss_scores.str("");
        stat_info.Merge(si);
        si.Clear();
    };

    auto dispatcher = dataset_.GetDispatcher();

    auto work_func = [&](size_t i) {
        Worker worker(*this);

        std::ostringstream oss_cread;
        std::ostringstream oss_scores;

        for (auto ids = dispatcher->Get(); ids.size() > 0; ids = dispatcher->Get()){
            worker.ResetCache(ids, ids.size());
            for (auto tid : ids) {
                if (worker.Correct(tid)) {
                    if ( worker.GetCorrected().size() > 0) {
                        SaveCRead(oss_cread, tid, worker.GetCorrected(), worker.GetTrueRange());
                        if (of_infos.is_open()) worker.SaveReadInfos(oss_scores, tid, dataset_.read_store_);
                    } else {
                        LOG(WARNING)("Failed to correct read(%s)", dataset_.read_store_.QueryNameById(tid).c_str());
                    }
                }
                worker.Clear();
            }
            
            if (oss_cread.tellp() > (int)flush_block) {
                combine_func(oss_cread, oss_scores, worker.stat_info);
            }
            progress.Forward(ids.size());
        }

        if (oss_cread.tellp() > 0) {
            combine_func(oss_cread, oss_scores, worker.stat_info);
        }
    };

 
    LOG(INFO)("thread size %zd, totalsize %d", opts_.thread_size, dataset_.read_ids_.size());
    if (of_cread.is_open()) {
        MultiThreadRun((size_t)opts_.thread_size, work_func);
    } else {
        LOG(INFO)("Failed to open file: %s", opts_.rread_fname_.c_str());
    }

    stat_info.Report();
}


void ReadCorrect::SaveCRead(std::ostream &os, int tid, const std::string &cread, const std::array<size_t,2> &range) {
    os << ">" << dataset_.read_store_.QueryNameById(tid) << " range=" << range[0] << "-" << range[1] << "\n" 
       <<  cread << "\n";
}


bool ReadCorrect::ExactFilter(const Alignment &r) {
    if (r.AlignSize() < (size_t)opts_.filter1_.min_aligned_length && 
        r.AlignSize() < r.TargetSize() * opts_.filter1_.min_aligned_length) return true;

    if (r.Identity() < opts_.min_identity_) return true;

    if (r.AlignSize() >= (size_t)opts_.filter1_.min_accept_aligned_length) return false;

    const double oh_rate = opts_.filter1_.max_overhang_rate;

    size_t t_overhang = std::max(size_t(r.TargetSize()*oh_rate), (size_t)opts_.filter1_.max_overhang);
    size_t q_overhang = std::max(size_t(r.QuerySize()*oh_rate), (size_t)opts_.filter1_.max_overhang);

    if (r.target_start > t_overhang && r.query_start > q_overhang) return true;
    if (r.target_end + t_overhang < r.TargetSize() && r.query_end + q_overhang < r.QuerySize()) return true;

    return false;
}

bool ReadCorrect::ExactFilter(const Alignment &r, const std::array<size_t,2> &trange) {
    size_t start = std::max(trange[0], r.target_start);
    size_t end = std::min(trange[1], r.target_end);
    auto align_size = start < end ? end - start : 0;

    if (align_size < (size_t)opts_.filter1_.min_aligned_length && 
        align_size < (trange[1]-trange[0]) * opts_.filter1_.min_aligned_rate) return true;
    
    if (r.Identity() < opts_.min_identity_) return true;

    if (align_size >= (size_t)opts_.filter1_.min_accept_aligned_length) return false;

    const double oh_rate = opts_.filter1_.max_overhang_rate;

    size_t t_overhang = std::max(size_t((trange[1]-trange[0])*oh_rate), (size_t)opts_.filter1_.max_overhang);
    size_t q_overhang = std::max(size_t(r.QuerySize()*oh_rate), (size_t)opts_.filter1_.max_overhang);

    if ( r.target_start > trange[0] + t_overhang && r.query_start > q_overhang) return true;
    if (r.target_end + t_overhang < trange[1] && r.query_end + q_overhang < r.QuerySize()) return true;

    return false;
}

bool ReadCorrect::Worker::GetAlignment(Seq::Id id, const Overlap* o, Alignment& al) {
    TimeCounter::Mark m(tc_align);
    const auto& tread = o->GetRead(id);
    const auto& qread = o->GetOtherRead(id);

    stat_info.total++;
    DEBUG_printf("start align %s %s\n", owner_.dataset_.read_store_.QueryNameById(qread.id).c_str(), owner_.dataset_.read_store_.QueryNameById(tread.id).c_str());
    if (owner_.opts_.use_cache) {
        if (!cache_.GetAlignment(qread.id, tread.id, o->SameDirect(), al)) {
            std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
            auto r = aligner_.Align(owner_.dataset_.read_store_.GetSeq(qread.id), !o->SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
            cache_.SetAlignment(qread.id, tread.id, o->SameDirect(), al);
            return r;

        } else {
            stat_info.cache++;
            return al.Valid();
        }
    } else {
        if (o->detail_.size() == 0) {
            std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
            return aligner_.Align(owner_.dataset_.read_store_.GetSeq(qread.id), !o->SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
        } else {
            return Cigar2Alignment(id, o, al);
        }
    }
}
      
std::array<size_t,2> MostEffectiveCoverage(size_t tsize, const std::vector<Alignment> &aligns, size_t stub, int min_coverage) {
    if (aligns.size() == 0) return {0, 0};

    assert(tsize >= stub*2);
    if (min_coverage == 0) return {stub, tsize - stub} ;
    
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
    if (start >= 0) {
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

std::vector<int> CalculateLocalDistanceThreshold(const std::vector<Alignment>& als, size_t cov, double identity) {

    std::vector<size_t> positions;
    for (auto &al : als ) {
        if (al.MaxLocalIdentity_100(1000) <= 98) {
            positions.push_back(al.MaxLocalDistancePosition());
        }
    }
    
    for (auto al : als ) {
        printf("pppp: %zd %0.02f\n", al.MaxLocalDistancePosition(), al.MaxLocalIdentity_100(1000));
    }

    std::sort(positions.begin(), positions.end());

    for (auto p : positions) {
        printf("pppp %zd, \n", p);
    }

    assert(als.size() > 0);
    size_t n = als[0].local_distances.size();
    std::vector<int> thresholds(n, -1);

    for (size_t i = 0; i < n; ++i) {
        std::vector<int> dis;
        for (const auto &al : als) {
            if (al.local_distances[i] >= 0) {
                dis.push_back(al.local_distances[i]);
                DEBUG_printf("ckck add %zd\n", al.local_distances[i]);
            }
        }

        if (dis.size() == 0) continue;


        if (dis.size() <= 10) {
            auto m = ComputeMeanAbsoluteDeviation(dis);
            thresholds[i] = m[0] + 3*1.253*m[1];
            DEBUG_printf("ckck mean th(%zd) = %d, %d, %d, %zd\n", i , thresholds[i], m[0], m[1], dis.size());
        } else {
            std::sort(dis.begin(), dis.end(), [](int a, int b) { return a < b; });
            std::vector<int> oks(dis.begin(), dis.begin() + std::min(dis.size(), cov));
            auto m = ComputeMedianAbsoluteDeviation(oks);
            thresholds[i] = m[0] + 3*1.4826*m[1];
            DEBUG_printf("ckck median th(%zd) = %d, %d, %d, %zd\n", i , thresholds[i], m[0], m[1], dis.size());
        }
    }
    return thresholds;
}

bool CheckLocalDistance(const Alignment &al, const std::vector<int> thresholds) {
    assert(al.local_distances.size() == thresholds.size());

    for (size_t i = 0; i < thresholds.size(); ++i) {
        DEBUG_printf("ckck CMP(%zd) %d < %d\n", i, al.local_distances[i] , thresholds[i]);
        if (thresholds[i] >= 0 && al.local_distances[i] >= 0) {
            if (al.local_distances[i] > thresholds[i]) {
                return false;
            }
        }
    }
    return true;
}

bool ReadCorrect::Worker::Correct(int id) {
    auto group = owner_.dataset_.GetOverlaps(id);
    if (group.Empty()) return false;

    const DnaSeq& target = owner_.dataset_.read_store_.GetSeq(id);
    assert(target.Size() >= (size_t)owner_.opts_.filter0_.min_length); 
    group.Sort(owner_.opts_.cands_opts_.overhang_weight);

    aligner_.SetTarget(target);
    std::vector<int> coverage(target.Size(), 0);

    std::vector<Alignment> first_als;
    for (size_t i = 0; i < group.Size(); ++i) {
        DEBUG_printf("done = %zd, group_size = %zd\n", i, group.Size());
        auto ol = group.Get(i, 0); 
        const auto& tread = ol->GetRead(id);
        const auto& qread = ol->GetOtherRead(id);
        Alignment al(tread.id, qread.id);
        bool r = false;
        double best_identity = 0.0;

        for (size_t j = 0; j < group.Size(i); ++j) {
            auto ol = group.Get(i, j);
            Alignment al_local(tread.id, qread.id);

            auto r_local = GetAlignment(id, ol, al_local);
            DEBUG_printf("alignment(%s-%s): r = %d, q = (%zd %zd %zd),  d=%d, t = (%zd %zd %zd), d=%zd,%f,  %zd\n", 
                owner_.dataset_.QueryStringById(qread.id).c_str(), owner_.dataset_.QueryStringById(tread.id).c_str(),
                r_local,
               al_local.query_start, al_local.query_end, al_local.QuerySize(), ol->SameDirect(),
               al_local.target_start, al_local.target_end, al_local.TargetSize(), al_local.distance, al_local.Identity(), al_local.local_distances.size());

            DEBUG_printf("al_global_ident %.02f <= %.02f\n", al_local.Identity(), owner_.opts_.min_identity_);
            if (r_local && !owner_.ExactFilter(al_local) && al_local.Identity() >= owner_.opts_.min_identity_) {
                al_local.ComputeDistance(owner_.opts_.local_window_size_);
                DEBUG_printf("al_local_ident %.02f <= %.02f\n", al_local.MaxLocalIdentity_100(owner_.opts_.local_window_size_), owner_.opts_.min_local_identity_);
                if (al_local.MaxLocalIdentity_100(owner_.opts_.local_window_size_) >= owner_.opts_.min_local_identity_) {
                    if (best_identity < al_local.Identity()) {
                        best_identity = al_local.Identity();
                        r = r_local;
                        al = al_local;
                        DEBUG_printf("ext d = %d\n", al.distance);
                    }
                    if (j >= 3) break;
                }

            }
        }

        if (r && !owner_.ExactFilter(al)) { 
            stat_info.succ++;
            first_als.push_back(al);
            std::for_each(coverage.begin()+al.target_start, coverage.begin()+al.target_end, [](int& c) {c++;} );
        } 

        if (owner_.opts_.cands_opts_.IsEndCondition(coverage)) break;
    }

    if (first_als.size() > 0) {
    ////////

    if (owner_.opts_.check_local_identity_) {
        std::vector<Alignment> first_als1 = CheckLocalDistance0(first_als);
        
        DEBUG_printf("ckck first_als size %zd\n", first_als.size());
        std::swap(first_als1, first_als);
        DEBUG_printf("ckck first_als size %zd\n", first_als.size());
    }


    size_t stub = 500;
    auto range = MostEffectiveCoverage(target.Size(), first_als, stub, owner_.opts_.min_coverage); 
    DEBUG_printf("Range: %zd - %zd\n", range[0], range[1]);
    if (range[0] < range[1] && range[1] - range[0] + 2*stub >= (size_t)owner_.opts_.filter0_.min_length) {
        assert(range[0] >= stub && range[1] + stub <= target.Size());
        range[0] -= stub;
        range[1] += stub;

        for (const auto &al : first_als) {
            if (!owner_.ExactFilter(al, range)) {
               aligned_.push_back(al);
            }
        }
        DEBUG_printf("aligned_.size: %zd\n", aligned_.size());
        for (auto &al : aligned_) { al.Rearrange(); }
        {
            TimeCounter::Mark m(tc_graph);
        graph_.Build(target, range, aligned_);
        graph_.Consensus();
        }
        return true;
    }
    }
    return false;
}

std::vector<Alignment>  ReadCorrect::Worker::CheckLocalDistance0(const std::vector<Alignment>& als) {

    std::vector<size_t> positions;
    for (auto &al : als ) {
        if (al.MaxLocalIdentity_100(1000) <= owner_.opts_.min_identity_) {
            positions.push_back(al.MaxLocalDistancePosition());
        }
    }
    
    std::sort(positions.begin(), positions.end());
    auto groups = GroupPositions(positions);

    std::unordered_set<size_t> removed;
    for (const auto gp : groups) {
        DEBUG_printf("al_local_group: %zd - %zd\n", gp[0], gp[1]);

        std::vector<std::pair<bool, uint16_t>> distances;
        std::vector<uint16_t> vdist;
        for (auto &al : als) {
            auto r = al.MaxLocalDistance(gp[0], gp[1]);
            distances.push_back(r);
            if (r.first) {
                vdist.push_back(r.second);
            }
            DEBUG_printf("al_local_group_i: %d - %zd\n", r.first, r.second);
        }
        
        size_t threshold = 1000;
        size_t cov = 60;
        if (vdist.size() == 0) continue;
        if (vdist.size() <= 10) {
            auto m = ComputeMeanAbsoluteDeviation(vdist);
            threshold = m[0] + 6*1.253*m[1];
            DEBUG_printf("al_local_th mean = %d, %d, %d, %zd\n" , threshold, m[0], m[1], vdist.size());
        } else {
            std::sort(vdist.begin(), vdist.end(), [](int a, int b) { return a < b; });
            std::vector<uint16_t> oks(vdist.begin(), vdist.begin() + std::min(vdist.size(), cov));
            auto m = ComputeMedianAbsoluteDeviation(oks);
            threshold = m[0] + 6*1.4826*m[1];
            DEBUG_printf("al_local_th median = %d, %d, %d, %zd\n", threshold, m[0], m[1], vdist.size());
        }

        for (size_t i = 0; i < distances.size(); ++i) {
            if (distances[i].second > std::max<size_t>(threshold, 10)) {
                removed.insert(i);
            }
        }
        DEBUG_printf("al_local_removed: %zd\n", removed.size());
    }

    DEBUG_printf("al_local_removed: %zd\n", removed.size());
    std::vector<Alignment> new_als;
    for (size_t i = 0; i < als.size(); ++i) {
        if (removed.find(i) == removed.end()) {
            new_als.push_back(als[i]);
        }
    }

    return new_als;
}

std::vector<std::array<size_t, 2>> ReadCorrect::Worker::GroupPositions(const std::vector<size_t> &sorted_positions) {
    std::vector<std::array<size_t, 2>> groups;
    //assert(sorted_positions.size() >= 1);

    std::array<int, 2> curr = {-1, -1 };
    for (auto p : sorted_positions) {
        if (curr[0] == -1) {
            curr = {p, p};
        } else {
            if (p - curr[1] < owner_.opts_.local_window_size_/2 && p - curr[0] <= owner_.opts_.local_window_size_) {
                curr[1] = p;
            } else {
                groups.push_back({curr[0], curr[1]});
                curr = {p , p};
            }

        }
    }
    if (curr[0] != -1) {
        groups.push_back({curr[0], curr[1]});
        curr = {-1, -1};
    }
    return groups;
}



bool ReadCorrect::Worker::Cigar2Alignment(Seq::Id tid, const Overlap* ol, Alignment &al) {
    assert(ol->detail_.size() > 0);

    const DnaSeq& qseq = owner_.dataset_.read_store_.GetSeq(ol->a_.id);
    const DnaSeq& tseq = owner_.dataset_.read_store_.GetSeq(ol->b_.id);
    assert(ol->b_.strand == 0);

    std::vector<uint8_t> tal;   tal.reserve(ol->AlignedLength()*2);
    std::vector<uint8_t> qal;   qal.reserve(ol->AlignedLength()*2);


    auto get_base = [](const Overlap::Read &r, const DnaSeq& seq, size_t idx) {
        return r.strand == 0 ? seq[r.start+idx] : (3 - seq[r.end - idx - 1]);
    };

    size_t qidx = 0;        // not from ol->b_.start;
    size_t tidx = ol->b_.start;
    for (const auto &d : ol->detail_) {
        switch (d.type){
        case 'M':
        case '=':
            for (size_t i = 0; i < d.len; ++i) {
                uint8_t cq = get_base(ol->a_, qseq, qidx+i);
                uint8_t ct = tseq[tidx+i];
                qal.push_back(cq+1);
                tal.push_back(ct+1);
            }
            qidx += d.len;
            tidx += d.len;
            break;
        case 'I':
            for (size_t i = 0; i < d.len; ++i) {
                char cq = get_base(ol->a_, qseq, qidx+i);
                qal.push_back(cq+1);
                tal.push_back(0);
            }
            qidx += d.len;
            break;
        case 'D':
            for (size_t i = 0; i < d.len; ++i) {
                char ct = tseq[tidx+i];
                qal.push_back(0);
                tal.push_back(ct+1);
            }
            tidx += d.len;
            break; 
        default:
            LOG(ERROR)("never come here");
        }
    }
    
    const auto& tread = ol->GetRead(tid);
    const auto& qread = ol->GetOtherRead(tid);
    al.query = &owner_.dataset_.read_store_.GetSeq(qread.id);
    al.target = &owner_.dataset_.read_store_.GetSeq(tread.id);
    al.target_start = tread.start;
    al.target_end = tread.end;
    al.query_start = qread.start;
    al.query_end = tread.end;

    const char* ACGT = "-ACGT-";
    if (tread.id == ol->b_.id) {
        for (size_t i = 0; i < tal.size(); ++i) {
            al.aligned_target.push_back(ACGT[tal[i]]);
            al.aligned_query.push_back(ACGT[qal[i]]);
        }
    printf("q:%s\nt:%s\n", al.aligned_target.c_str(), al.aligned_query.c_str());
    } else {
        for (size_t i = 0; i < tal.size(); ++i) {
            al.aligned_target.push_back(ACGT[5 - tal[tal.size()-i-1]]);
            al.aligned_query.push_back(ACGT[5 - qal[tal.size()-i-1]]);
        }
    }


    return true;
}
   
} // namespace fsa {
    