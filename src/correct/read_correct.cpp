#include "read_correct.hpp"

#include <iostream>
#include <atomic>
#include <ctime>
#include "./utils/logger.hpp"
#include "../utility.hpp"
#include "edlib.h"

namespace fsa {

TimeCounter tc_align("get_al");
TimeCounter tc_graph("graph");;
TimeCounter tc_al_ava("al_ava");
TimeCounter tc_al_cigar("al_cigar");
TimeCounter tc_al_map("al_map");
TimeCounter tc_al_head_tail("al_head_tail");

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
    std::ofstream of_flt_ols(opts_.filtered_overlaps_fname_);
    const size_t flush_block = 20*1024*1024;
    StatInfo stat_info;

    Progress progress(5000, dataset_.read_ids_.size());

    auto combine_func = [&](std::ostringstream &oss_cread, std::ostringstream &oss_scores, std::ostringstream& oss_flt_ols, StatInfo &si) {
        std::lock_guard<std::mutex> lock(mutex);
        of_cread << oss_cread.str();              
        oss_cread.str("");

        if (of_infos.is_open()) {
            of_infos << oss_scores.str();
        }
        oss_scores.str("");

        if (of_flt_ols.is_open()) {
            of_flt_ols << oss_flt_ols.str();
        }
        oss_flt_ols.str("");
        
        stat_info.Merge(si);
        si.Clear();
    };

    auto dispatcher = dataset_.GetDispatcher();

    auto work_func = [&](size_t i) {
        LOG(INFO)("Start thread %zd", i);
        Worker worker(*this);

        std::ostringstream oss_cread;
        std::ostringstream oss_scores;
        std::ostringstream oss_flt_ols;

        for (auto ids = dispatcher->Get(); ids.size() > 0; ids = dispatcher->Get()){
            worker.ResetCache(ids, ids.size());
            for (auto tid : ids) {
                
                //LOG(INFO)("START: %s", dataset_.QueryStringById(tid).c_str());
                if (Correct(tid, worker)) {
                    if ( worker.GetCorrected().size() > 0) {
                        SaveCRead(oss_cread, tid, worker.GetCorrected(), worker.GetTrueRange());
                        if (of_infos.is_open()) worker.SaveReadInfos(oss_scores, tid, dataset_.read_store_);
                        if (of_flt_ols.is_open()) worker.DumpFilteredOverlaps(oss_flt_ols);
                    } else {
                        LOG(WARNING)("Failed to correct read(%s)", dataset_.read_store_.QueryNameById(tid).c_str());
                    }
                }
                worker.Clear();
                //LOG(INFO)("End: %s", dataset_.QueryStringById(tid).c_str());
            }
            
            if (oss_cread.tellp() > (int)flush_block) {
                combine_func(oss_cread, oss_scores, oss_flt_ols, worker.stat_info);
            }
            progress.Forward(ids.size());
        }

        if (oss_cread.tellp() > 0) {
            combine_func(oss_cread, oss_scores, oss_flt_ols, worker.stat_info);
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

    DEBUG_printf("cov(0): %d\n", coverage[0]);
    for (size_t i=1; i< coverage.size(); ++i) {
        coverage[i] += coverage[i-1];
        DEBUG_printf("cov(%zd): %d\n", i, coverage[i]);
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

bool ReadCorrect::Correct(Seq::Id id, Worker& wrk) {
    auto group = dataset_.GetOverlaps(id);
    if (group.Empty()) return false;
    
    const DnaSeq& target = dataset_.read_store_.GetSeq(id);
    assert(target.Size() >= (size_t)opts_.filter0_.min_length); 
    group.Sort(opts_.cands_opts_.overhang_weight);

    wrk.aligner_.SetTarget(target);
    std::vector<int> coverage(target.Size(), 0);

    std::vector<Alignment> first_als;
    
    DEBUG_printf("start_correcting tid = %s(%d), groupsize = %zd\n", dataset_.QueryStringById(id).c_str(), id, group.Size());
    //for (size_t i = 0; i < std::min<size_t>(1000, group.Size()); ++i) {
    for (size_t i = 0; i < group.Size(); ++i) {
        
        auto al = GetAlignmentWithCache(id, group, i, wrk);

        if (al.Valid()) { 
            first_als.push_back(al);
            std::for_each(coverage.begin()+al.target_start, coverage.begin()+al.target_end, [](int& c) {c++;} );
        } 

        if (opts_.cands_opts_.IsEndCondition(coverage)) break;
    }
    
    DEBUG_printf("alignsize = %zd, %zd\n", first_als.size(), wrk.aligned_.size());

    if (first_als.size() > 0) {
    ////////

    if (opts_.check_local_identity_) {
        std::vector<Alignment> first_als1 = wrk.CheckLocalDistance0(first_als);
        
        DEBUG_printf("ckck first_als size %zd\n", first_als.size());
        std::swap(first_als1, first_als);
        DEBUG_printf("ckck first_als size %zd\n", first_als.size());
    }


    size_t stub = 500;
    auto range = MostEffectiveCoverage(target.Size(), first_als, stub, opts_.min_coverage); 
    DEBUG_printf("Range: %zd - %zd\n", range[0], range[1]);
    if (range[0] < range[1] && range[1] - range[0] + 2*stub >= (size_t)opts_.filter0_.min_length) {
        assert(range[0] >= stub && range[1] + stub <= target.Size());
        range[0] -= stub;
        range[1] += stub;

        for (const auto &al : first_als) {
            if (!ExactFilter(al, range)) {
               wrk.aligned_.push_back(al);
            }
        }
        DEBUG_printf("aligned_.size: %zd\n", wrk.aligned_.size());
        for (auto &al : wrk.aligned_) { al.Rearrange(); }
        {
            
            TimeCounter::Mark m(tc_graph);
        wrk.graph_.Build(target, range, wrk.aligned_);
        wrk.graph_.Consensus();
        }
        
        return true;
    }
    }
    return false;
}


Alignment ReadCorrect::GetAlignmentWithCache(Seq::Id tid, const CrrDataset::OlGroup& group, size_t ig, Worker& wrk) {
    
    TimeCounter::Mark m(tc_align);
    auto ol = group.Get(ig, 0); 
    const auto& tread = ol->GetRead(tid);
    const auto& qread = ol->GetOtherRead(tid);

    DEBUG_printf("doing=%zd/%zd, qid=%s\n", ig, group.Size(), dataset_.QueryStringById(qread.id).c_str());
    Alignment al(tread.id, qread.id);
    al.strand = ol->SameDirect() ? 0 : 1;

    wrk.stat_info.total++;
    // 从cache查询
    if (opts_.use_cache) {
        
        if (!wrk.cache_.GetAlignment(qread.id, tread.id, al)) {
            al = GetAlignmentOnes(tid, group, ig, wrk);
            wrk.cache_.SetAlignment(qread.id, tread.id, al);
        } else {
            wrk.stat_info.cache++;
        }
    } else {{
        al = GetAlignmentOnes(tid, group, ig, wrk);
    }}
    return al;
}

Alignment ReadCorrect::GetAlignmentOnes(Seq::Id tid, const CrrDataset::OlGroup& group, size_t ig, Worker& wrk) {

    auto ol = group.Get(ig, 0); 
    const auto& tread = ol->GetRead(tid);
    const auto& qread = ol->GetOtherRead(tid);

    Alignment al_best(tread.id, qread.id);

    for (size_t j = 0; j < group.Size(ig); ++j) {
        const Overlap& ol = *group.Get(ig, j);
        DEBUG_printf("subgroup(%zd) = %zd: %s <-> %s: %lld, %s\n", j, group.Size(ig), dataset_.QueryStringById(ol.a_.id).c_str(),
             dataset_.QueryStringById(ol.b_.id).c_str(), ol, ol.ToM4Line().c_str());
        Alignment al = GetAlignmentOne(tid, ol, wrk);
        DEBUG_printf("alignment(%s<->%s): (%d, %d, %d) - (%d, %d, %d) %f\n", dataset_.QueryStringById(ol.a_.id).c_str(), dataset_.QueryStringById(ol.b_.id).c_str(),
            al.query_start, al.query_end, al.QuerySize(), al.target_start, al.target_end, al.TargetSize(), al.Identity());
        
        DEBUG_printf("q:%s\nt:%s\n", al.aligned_query.c_str(), al.aligned_target.c_str());
        DEBUG_printf("check_global: %d, %f > %f\n", ExactFilter(al), al.Identity(), opts_.min_identity_);
        if (!ExactFilter(al) && al.Identity() >= opts_.min_identity_) {     
            al.ComputeDistance(opts_.local_window_size_);
            DEBUG_printf("check_local: %f > %f\n", al.MaxLocalIdentity_100(opts_.local_window_size_), opts_.min_local_identity_);
            if (al.MaxLocalIdentity_100(opts_.local_window_size_) >= opts_.min_local_identity_) {
                // pass check
                if (al_best.Identity() < al.Identity()) {
                DEBUG_printf("exchange(%zd) = %zd: %s <-> %s\n", j, group.Size(ig), dataset_.QueryStringById(ol.a_.id).c_str(),
                    dataset_.QueryStringById(ol.b_.id).c_str());
                //if (al_best.AlignSize() < al.AlignSize()) {
                    al_best = al;
                }
            }
        }
        if (j >= 15) {
            break;
        }
    }
    return al_best;
}

Alignment ReadCorrect::GetAlignmentOne(Seq::Id tid, const Overlap &ol, Worker& wrk) {
    const auto& tread = ol.GetRead(tid);
    const auto& qread = ol.GetOtherRead(tid);
    Alignment al(tread.id, qread.id);
    al.strand = ol.SameDirect() ? 0 : 1;
    al.query = &dataset_.read_store_.GetSeq(ol.a_.id);
    al.target = &dataset_.read_store_.GetSeq(ol.b_.id);
    
    if (ol.detail_.size() != 0) {
        TimeCounter::Mark m(tc_al_cigar);
        GetAlignmentFromCigar(tid, ol, al);
    } else if (false && CrrDataset::OlGroup::GetType(ol) == CrrDataset::OlGroup::MAP) {
        TimeCounter::Mark m(tc_al_map);
        GetAlignmentFromMapping1(tid, ol, wrk, al);
    } else {
        TimeCounter::Mark m(tc_al_ava);
        std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
        // TODO target 由调用者设置，可能存在不一致，需要优化。
        wrk.aligner_.Align(dataset_.read_store_.GetSeq(qread.id), !ol.SameDirect(), range, al);  
    }
    return al;
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
            DEBUG_printf("al_local_group_i(%d-%d): %d - %zd\n", al.qid, al.tid, r.first, r.second);
        }
        
        size_t threshold = 1000;
        size_t cov = 20;
        if (vdist.size() == 0) continue;
        if (vdist.size() <= 10) {
            auto m = ComputeMeanAbsoluteDeviation(vdist);
            threshold = m[0] + 3*1.253*m[1];
            DEBUG_printf("al_local_th mean = %d, %d, %d, %zd\n" , threshold, m[0], m[1], vdist.size());
        } else {
            std::sort(vdist.begin(), vdist.end(), [](int a, int b) { return a < b; });
            std::vector<uint16_t> oks(vdist.begin(), vdist.begin() + std::min(vdist.size(), cov));
            for (auto o : oks) {
                DEBUG_printf("al_local_th ok %zd\n", o);
            }
            auto m = ComputeMedianAbsoluteDeviation(oks);
            threshold = m[0] + 3*1.4826*m[1];
            DEBUG_printf("al_local_th median1 = %d, %d, %d, %zd\n", threshold, m[0], m[1], vdist.size());
            for (auto t : vdist) {
                if (t > threshold) {
                    if (t <= threshold + std::max(m[1], m[0])) {
                        threshold = t;
                    } else {
                        break;
                    }
                }
            }
            DEBUG_printf("al_local_th median2 = %d, %d, %d, %zd\n", threshold, m[0], m[1], vdist.size());
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
    const int half_win = owner_.opts_.local_window_size_ / 2;
    std::vector<std::array<size_t, 2>> groups;

    if (sorted_positions.size() > 0) {
    //assert(sorted_positions.size() >= 1);

        std::array<int, 2> curr = {sorted_positions[0], sorted_positions[0] };
        for (auto p : sorted_positions) {

            //if (p - curr[1] < half_win && p - curr[0] <= half_win*2) {
            if (p - curr[1] < half_win) {
                curr[1] = p;
            } else {
                groups.push_back({curr[0], curr[1]});
                curr = {p , p};
            }
        }
        groups.push_back({curr[0] - half_win, curr[1] + half_win});
    }
    return groups;
}

void ReadCorrect::Worker::DumpFilteredOverlaps(std::ostream &os) {
    StringPool::UnsafeNameId ni(owner_.dataset_.GetStringPool());
    for (const auto &al : aligned_) {
        Overlap ol;
        ol.a_.id = al.qid;
        ol.a_.len = al.query->Size();
        if (al.strand == 0) {
            ol.a_.start = al.query_start;
            ol.a_.end = al.query_end;
        } else {
            ol.a_.start = al.query->Size() - al.query_end;
            ol.a_.end = al.query->Size() - al.query_start;
        }
        ol.a_.strand = al.strand;

        ol.b_.id = al.tid;
        ol.b_.len = al.target->Size();
        ol.b_.start = al.target_start;
        ol.b_.end = al.target_end;
        ol.b_.strand = 0;

        ol.identity_ = al.Identity();
        
        os << OverlapStore::ToPafLine(ol, ni) << "\n";
    }
}

void ReadCorrect::GetAlignmentFromCigar(Seq::Id tid, const Overlap& ol, Alignment &al) {
    assert(ol.detail_.size() > 0);

    const DnaSeq& qseq = dataset_.read_store_.GetSeq(ol.a_.id);
    const DnaSeq& tseq = dataset_.read_store_.GetSeq(ol.b_.id);
    assert(ol.b_.strand == 0);

    std::vector<uint8_t> tal;   tal.reserve(ol.AlignedLength()*2);
    std::vector<uint8_t> qal;   qal.reserve(ol.AlignedLength()*2);


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
    al.query = &dataset_.read_store_.GetSeq(qread.id);
    al.target = &dataset_.read_store_.GetSeq(tread.id);
    assert(al.target!= nullptr);
    al.target_start = tread.start;
    al.target_end = tread.end;
    al.query_start = qread.start;
    al.query_end = tread.end;
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
TimeCounter tc_map("mapping");
void ReadCorrect::GetAlignmentFromMapping0(Seq::Id tid, const Overlap& ol, Worker& wrk, Alignment &al) {
    TimeCounter::Mark m(tc_map);
    const auto& tread = ol.GetRead(tid);
    const auto& qread = ol.GetOtherRead(tid);

    std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
    wrk.aligner_.Align(dataset_.read_store_.GetSeq(qread.id), !ol.SameDirect(), range, al); 
}

void ReadCorrect::GetAlignmentFromMapping1(Seq::Id tid, const Overlap& ol, Worker& wrk, Alignment &al) {
    const auto& tread = ol.GetRead(tid);
    const auto& qread = ol.GetOtherRead(tid);
    const auto& qseq = dataset_.read_store_.GetSeq(qread.id);
    const auto& tseq = dataset_.read_store_.GetSeq(tread.id);
    
    assert(CrrDataset::OlGroup::GetType(ol) == CrrDataset::OlGroup::MAP);
    auto pair = static_cast<const Mapping::Pair&>(ol);

    al.target_start = tread.start;
    al.target_end = tread.end;
    // the query position is at reverse-complement sequence
    al.query_start = ol.SameDirect() ? qread.start : qread.len - qread.end;
    al.query_end = ol.SameDirect() ? qread.end : qread.len - qread.start;

    const char* ACGT = "ACGT";
    size_t it = al.target_start;
    size_t iq = 0;
    auto get_query_base = [&qseq, &qread, &pair](size_t p) {
        return pair.SameDirect() ? qseq[qread.start + p] : 3 - qseq[qread.end - 1 - p];
    };

    auto aligned = pair.AlignBases(tid, qseq, tseq);
    al.aligned_query.reserve(aligned.size());
    al.aligned_target.reserve(aligned.size());
    for (auto i : aligned) {
        if (i == EDLIB_EDOP_MATCH) {
            al.aligned_query .push_back(ACGT[get_query_base(iq++)]);
            al.aligned_target.push_back(ACGT[tseq[it++]]);
        } else if (i == EDLIB_EDOP_INSERT) {
            al.aligned_query .push_back(ACGT[get_query_base(iq++)]);
            al.aligned_target.push_back('-');
            al.distance++;

        } else if (i == EDLIB_EDOP_DELETE) {
            al.aligned_query .push_back('-');
            al.aligned_target.push_back(ACGT[tseq[it++]]);
            al.distance++;
            
        } else {
            assert(i == EDLIB_EDOP_MISMATCH);
            al.aligned_query .push_back(ACGT[get_query_base(iq++)]);
            al.aligned_target.push_back(ACGT[tseq[it++]]);
            al.distance++;
        }
    }
    assert(it == al.target_end);
    assert(iq == al.query_end - al.query_start);

    if (ol.SameDirect()) {
        TimeCounter::Mark m(tc_al_head_tail);
        // head
        auto thead = tseq.ToUInt8(0, al.target_start);
        auto qhead = qseq.ToUInt8(0, al.query_start);
        Alignment al_h;
        auto r_h = wrk.aligner_.GetWorker()->Align(
            (const char*)&qhead[0], qhead.size(), 
            (const char*)&thead[0], thead.size(), 
            {al.query_start, al.query_start},
            {al.target_start, al.target_start}, 
            al_h);
        if (r_h) {
            al.distance += al_h.distance;
            al.query_start = al_h.query_start;
            al.target_start = al_h.target_start;
            al.aligned_query = al_h.aligned_query + al.aligned_query;
            al.aligned_target = al_h.aligned_target + al.aligned_target;
        }

        // tail
        auto tseq_t = tseq.ToUInt8(al.target_end);
        auto qseq_t = qseq.ToUInt8(al.query_end);
        Alignment al_t;
        auto r_t = wrk.aligner_.GetWorker()->Align(
            (const char*)&qseq_t[0], qseq_t.size(), 
            (const char*)&tseq_t[0], tseq_t.size(), 
            {0, 0},
            {0, 0}, 
            al_t);
        if (r_t) {
            al.distance += al_t.distance;
            al.query_end += al_t.query_end;
            al.target_end += al_t.target_end;
            al.aligned_query += al_t.aligned_query;
            al.aligned_target += al_t.aligned_target;
        }

    } else {
        TimeCounter::Mark m(tc_al_head_tail);
        // head
        auto tseq_h = tseq.ToUInt8(0, tread.start);
        auto qseq_h = qseq.ToUInt8(qread.end, -1, true);
        Alignment al_h;
        auto r_h = wrk.aligner_.GetWorker()->Align(
            (const char*)&qseq_h[0], qseq_h.size(), 
            (const char*)&tseq_h[0], tseq_h.size(), 
            {qread.len - qread.end, qread.len - qread.end},
            {tread.start, tread.start}, 
            al_h);
        if (r_h) {
            al.distance += al_h.distance;
            al.query_start = al_h.query_start;
            al.target_start = al_h.target_start;
            al.aligned_query = al_h.aligned_query + al.aligned_query;
            al.aligned_target = al_h.aligned_target + al.aligned_target;
        }

        // tail
        auto tseq_t = tseq.ToUInt8(tread.end);
        auto qseq_t = qseq.ToUInt8(0, qread.start, true);
        Alignment al_t;
        auto r_t = wrk.aligner_.GetWorker()->Align(
            (const char*)&qseq_t[0], qseq_t.size(), 
            (const char*)&tseq_t[0], tseq_t.size(), 
            {0, 0},
            {0, 0}, 
            al_t);
        if (r_t) {
            al.distance += al_t.distance;
            al.query_end += al_t.query_end;
            al.target_end += al_t.target_end;
            al.aligned_query += al_t.aligned_query;
            al.aligned_target += al_t.aligned_target;
        }

    }
    
}


} // namespace fsa {
    