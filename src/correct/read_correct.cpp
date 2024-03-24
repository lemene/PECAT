#include "read_correct.hpp"

#include <iostream>
#include <atomic>
#include <ctime>
#include <random>
#include "./utils/logger.hpp"
#include "../utility.hpp"

namespace fsa {


ReadCorrect::ReadCorrect() {
}

ArgumentParser ReadCorrect::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}

void ReadCorrect::CheckArguments() {
    opts_.filter0_.From(opts_.filter0_opts_);
    opts_.filter1_.From(opts_.filter1_opts_);
    
    opts_.filter0_opts_ = opts_.filter0_.ToString();
    opts_.filter1_opts_ = opts_.filter1_.ToString();

    opts_.cands_opts_.From(opts_.cands_opts_str_);          // 合并用户设置
    opts_.cands_opts_str_ = opts_.cands_opts_.ToString();   // 输出所有参数
    if (opts_.debug) SetDebug();
}

void ReadCorrect::Running() {
    LoadReadIds();

    LoadOverlaps(opts_.overlap_fname_);
    LoadReads();

    dataset_.grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>(dataset_.read_ids_.begin(), dataset_.read_ids_.end()));

    if (opts_.use_cache) GroupReadIds();

    EstimateParameters();
    Correct();
}

void ReadCorrect::LoadReads() {

    std::unordered_set<Seq::Id> ids;
    for (size_t i = 0; i < dataset_.ol_store_.Size(); ++i) {
        const Overlap& o = dataset_.ol_store_.Get(i);
        ids.insert(o.a_.id);
        ids.insert(o.b_.id);
    }
    dataset_.read_store_.Load(opts_.rread_fname_, "", false, ids);
    
    if (dataset_.read_ids_.empty()) {
        dataset_.read_ids_.assign(ids.begin(), ids.end());
    }
}

void ReadCorrect::LoadOverlaps(const std::string &fname) {
    std::unordered_set<Seq::Id> ids(dataset_.read_ids_.begin(), dataset_.read_ids_.end());
    
    dataset_.ol_store_.Load(fname, "", (size_t)opts_.thread_size, [this, &ids](Overlap &o) {
        bool rel = ids.empty() || ids.find(o.a_.id) != ids.end() || ids.find(o.b_.id) != ids.end();
        return rel && opts_.filter0_.Valid(o);
    });

    LOG(INFO)("Load %zd overlaps from file %s", dataset_.ol_store_.Size(), fname.c_str());

    dataset_.Load(dataset_.ol_store_.GetStringPool());
}



void ReadCorrect::LoadReadIds() {
    // read_store_.Load(rread_fname_, "", false);
    std::unordered_set<Seq::Id> ids;    // Remove duplicate names
    if (!opts_.read_name_.empty()) {
        ids.insert(dataset_.string_pool_.GetIdByStringUnsafe(opts_.read_name_));
    } else if (!opts_.read_name_fname_.empty()) {
        std::ifstream file(opts_.read_name_fname_);
        std::string line;
        while (std::getline(file, line)) {
            ids.insert(dataset_.string_pool_.GetIdByStringUnsafe(line));
        }
    } else {
        // correct all reads in read file
    }
    dataset_.read_ids_.assign(ids.begin(), ids.end());
}

void ReadCorrect::Correct() {
    LOG(INFO)("Start Correcting");

    std::mutex mutex;
    
    std::ofstream of_cread(opts_.cread_fname_);
    std::ofstream of_infos(opts_.infos_fname_);
    const bool save_infos = of_infos.is_open();
    const size_t flush_block = 20*1024*1024;

    Progress progress(*this);

    auto save_oss = [&](std::ostringstream &oss_cread, std::ostringstream &oss_scores) {
        std::lock_guard<std::mutex> lock(mutex);
        of_cread << oss_cread.str();              
        oss_cread.str("");

        if (save_infos) {
            of_infos << oss_scores.str();
            oss_scores.str("");
        }
    };

    struct Dispatcher { virtual Seq::Id Get(bool &uc) = 0; };
    struct SimpleDispatcher : public Dispatcher {
        SimpleDispatcher(Progress &p) : progress(p) {}
        Seq::Id Get(bool &uc) { uc = false; return progress.Get(); }
        Progress &progress;
    };

    struct GroupDispatcher : public Dispatcher {
        GroupDispatcher(Progress &p, ReadCorrect::Worker &w) : progress(p), worker(w), ids(10000){}

        Seq::Id Get(bool &uc) {
            uc = true; 
            if (index < size) {
                return ids[index++];
            } else {
                index = 0;
                size = progress.Get(ids);
                worker.ResetCache(ids, size);
                return index < size ? ids[index++] : -1;
            }
        }

        Progress &progress;
        ReadCorrect::Worker &worker;

        std::vector<Seq::Id> ids;
        size_t index { 0 };
        size_t size { 0 };

    };

    auto work_func = [&](size_t i) {
        Worker worker(*this);
        std::unique_ptr<Dispatcher> dispatcher(opts_.use_cache ? 
            (Dispatcher*)new GroupDispatcher(progress, worker) : 
            (Dispatcher*)new SimpleDispatcher(progress));

        std::ostringstream oss_cread;
        std::ostringstream oss_scores;
        bool uc = false;
        for (auto tid=dispatcher->Get(uc); tid != Seq::NID; tid = dispatcher->Get(uc)) {
            if (worker.Correct(tid, uc)) {
                if ( worker.GetCorrected().size() > 0) {
                    SaveCRead(oss_cread, tid, worker.GetCorrected(), worker.GetTrueRange());
                    if (save_infos) worker.SaveReadInfos(oss_scores, tid, dataset_.read_store_);
                } else {
                    LOG(WARNING)("Corrected Read(%s) is emtpy", dataset_.read_store_.QueryNameById(tid).c_str());
                }
            }
            worker.Clear();
            
            if (oss_cread.tellp() > (int)flush_block) {
                save_oss(oss_cread, oss_scores);
            }
        }

        if (oss_cread.tellp() > 0) {
            save_oss(oss_cread, oss_scores);
        }
        CollectWorkerInfo(worker, mutex);
    };

 
    LOG(INFO)("thread size %zd, totalsize %d", opts_.thread_size, dataset_.read_ids_.size());
    if (of_cread.is_open()) {
        MultiThreadRun((size_t)opts_.thread_size, work_func);
    } else {
        LOG(INFO)("Failed to open file: %s", opts_.rread_fname_.c_str());
    }

    Report();
}


void ReadCorrect::SaveCRead(std::ostream &os, int tid, const std::string &cread, const std::array<size_t,2> &range) {
    
    os << ">" << dataset_.read_store_.QueryNameById(tid) << " range=" << range[0] << "-" << range[1] << "\n" 
       <<  cread << "\n";
    
}

void ReadCorrect::GroupReadIds() {
    std::unordered_map<Seq::Id, bool> done;
    std::unordered_map<Seq::Id, int> lens;

    for (auto i : dataset_.read_ids_) {
        done[i] = false;
        lens[i] = dataset_.read_store_.GetSeqLength(i);
    }

    std::sort(dataset_.read_ids_.begin(), dataset_.read_ids_.end(), [&lens](int a, int b) {return lens[a] > lens[b]; });

    dataset_.group_ticks.push_back(0);

    for (auto i : dataset_.read_ids_) {
        if (!done[i]) {
            size_t s = dataset_.grouped_ids_.size();
            dataset_.grouped_ids_.push_back(i);
            done[i] = true;

            while (s < dataset_.grouped_ids_.size()) {
                auto gp = dataset_.grouper_.Get(dataset_.grouped_ids_[s]);
                if (!gp.Empty()) {
                    for (size_t igp = 0; igp < gp.Size(); ++igp) {
                        auto o = gp.Get(igp, 0);
                        auto d = done.find(o->GetOtherRead(gp.id).id); 
                        if (d != done.end() && !d->second) {
                            dataset_.grouped_ids_.push_back(o->GetOtherRead(gp.id).id);
                            d->second = true;
                        }
                    }
                }
                s++;
                if (dataset_.grouped_ids_.size() >= dataset_.group_ticks.back() + dataset_.group_size) {
                    break;
                }
            }
            if (dataset_.grouped_ids_.size() >= dataset_.group_ticks.back() + dataset_.group_size / 2) {
                dataset_.group_ticks.push_back(dataset_.grouped_ids_.size());
            }
        }
    }

    if (dataset_.grouped_ids_.size() > dataset_.group_ticks.back()) {
        dataset_.group_ticks.push_back(dataset_.grouped_ids_.size());
    }

    for (auto t : dataset_.group_ticks) {
        LOG(INFO)("TICK: %zd", t);
    }

}

void ReadCorrect::EstimateParameters() {

    size_t count = std::min<size_t>(10, dataset_.read_ids_.size());

    std::unordered_set<int> tests;

    std::default_random_engine e;
    std::uniform_int_distribution<int> u(0, dataset_.read_ids_.size()-1);
    e.seed(time(0));
    
    while (tests.size() < count) {
        tests.insert(dataset_.read_ids_[u(e)]);
    }

    auto aligner = ToolAligner::Create("edlib");

    double max_idt = 0.0;
    
    for (auto id : tests) {
        auto group = dataset_.grouper_.Get(id);
        if (group.Empty()) continue;

        for (size_t i = 0; i < group.Size() && i < 10; ++i) {
            auto o = group.Get(i, 0);
            const auto& tread = o->GetRead(id);
            const auto& qread = o->GetOtherRead(id);
            std::vector<uint8_t> tseq = dataset_.read_store_.GetSeq(tread.id).ToUInt8(tread.start, tread.end, false);
            std::vector<uint8_t> qseq = dataset_.read_store_.GetSeq(qread.id).ToUInt8(qread.start, qread.end, !o->SameDirect());
            Alignment al;
            auto r = aligner->Align((const char*)&qseq[0], qseq.size(), (const char*)&tseq[0], tseq.size(), {0, qseq.size()}, {0, tseq.size()}, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
            if (r && al.AlignSize() > al.TargetSize() / 2) {
                if (al.Identity() > max_idt) {
                    max_idt = al.Identity();
                }

            }
        }
    }

    double min_idt = 0.0, min_lc_idt = 0.0;
    if (max_idt >= 0.99) {
        min_idt = 95;
        min_lc_idt = 90;
    } else if (max_idt >= 0.95) {
        min_idt = 90;
        min_lc_idt = 80;
    } else if (max_idt >= 0.85) {
        min_idt = 75;
        min_lc_idt = 65;
    } else {
        min_idt = 60;
        min_lc_idt = 50;
    }
    opts_.min_identity_ = opts_.min_identity_ < 0 ? min_idt : opts_.min_identity_;
    opts_.min_local_identity_ = opts_.min_local_identity_ < 0 ? min_lc_idt : opts_.min_local_identity_;
    LOG(INFO)("Estimate parameters(%0.02f): min_identity = %f min_local_identity = %f", max_idt, opts_.min_identity_, opts_.min_local_identity_);
}

bool ReadCorrect::Worker::ExactFilter(const Alignment &r, const std::array<size_t,2> &trange) {
    size_t start = std::max(trange[0], r.target_start);
    size_t end = std::min(trange[1], r.target_end);
    auto align_size = start < end ? end - start : 0;

    if (align_size < (size_t)owner_.opts_.filter1_.min_aligned_length && 
        align_size < (trange[1]-trange[0]) * owner_.opts_.filter1_.min_aligned_rate) return true;
    
    if (r.Identity() < owner_.opts_.min_identity_) return true;

    if (align_size >= (size_t)owner_.opts_.filter1_.min_accept_aligned_length) return false;

    const double oh_rate = owner_.opts_.filter1_.max_overhang_rate;

    size_t t_overhang = std::max(size_t((trange[1]-trange[0])*oh_rate), (size_t)owner_.opts_.filter1_.max_overhang);
    size_t q_overhang = std::max(size_t(r.QuerySize()*oh_rate), (size_t)owner_.opts_.filter1_.max_overhang);

    if ( r.target_start > trange[0] + t_overhang && r.query_start > q_overhang) return true;
    if (r.target_end + t_overhang < trange[1] && r.query_end + q_overhang < r.QuerySize()) return true;

    return false;
}

bool ReadCorrect::Worker::ExactFilter(const Alignment &r) {
    if (r.AlignSize() < (size_t)owner_.opts_.filter1_.min_aligned_length && 
        r.AlignSize() < r.TargetSize() * owner_.opts_.filter1_.min_aligned_length) return true;

    if (r.Identity() < owner_.opts_.min_identity_) return true;

    if (r.AlignSize() >= (size_t)owner_.opts_.filter1_.min_accept_aligned_length) return false;

    const double oh_rate = owner_.opts_.filter1_.max_overhang_rate;

    size_t t_overhang = std::max(size_t(r.TargetSize()*oh_rate), (size_t)owner_.opts_.filter1_.max_overhang);
    size_t q_overhang = std::max(size_t(r.QuerySize()*oh_rate), (size_t)owner_.opts_.filter1_.max_overhang);

    if (r.target_start > t_overhang && r.query_start > q_overhang) return true;
    if (r.target_end + t_overhang < r.TargetSize() && r.query_end + q_overhang < r.QuerySize()) return true;

    return false;
}

bool ReadCorrect::Worker::GetAlignment(Seq::Id id, const Overlap* o, bool uc, Alignment& al) {
    const auto& tread = o->GetRead(id);
    const auto& qread = o->GetOtherRead(id);

    DEBUG_printf("start align %s %s\n", owner_.dataset_.read_store_.QueryNameById(qread.id).c_str(), owner_.dataset_.read_store_.QueryNameById(tread.id).c_str());
    if (owner_.opts_.use_cache && uc) {
        if (!cache_.GetAlignment(qread.id, tread.id, o->SameDirect(), al)) {
            std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
            auto r = aligner_.Align(owner_.dataset_.read_store_.GetSeq(qread.id), !o->SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
            cache_.SetAlignment(qread.id, tread.id, o->SameDirect(), al);
            return r;

        } else {
            return al.Valid();
        }
    } else {
        std::array<int, 4> range = {qread.start, qread.end, tread.start, tread.end};
        return aligner_.Align(owner_.dataset_.read_store_.GetSeq(qread.id), !o->SameDirect(), range, al);  // TODO target 由调用者设置，可能存在不一致，需要优化。
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
        if (thresholds[i] >= 0 && al.local_distances[i] >= 0) {
            DEBUG_printf("ckck CMP %d %d\n", al.local_distances[i] , thresholds[i]);
            if (al.local_distances[i] > thresholds[i]) {
                return false;
            }
        }
    }
    return true;
}


bool ReadCorrect::Worker::Correct(int id, bool uc) {
    auto group = owner_.dataset_.grouper_.Get(id);
    if (group.Empty()) return false;

    const DnaSeq& target = owner_.dataset_.read_store_.GetSeq(id);
    assert(target.Size() >= (size_t)owner_.opts_.filter0_.min_length); 

    std::vector<std::tuple<const Overlap*, double, size_t>> cands(group.Size());
    for (size_t i = 0; i < group.Size(); ++i) {
        std::get<0>(cands[i]) = group.Get(i, 0);
        std::get<1>(cands[i]) = 0.0;
        std::get<2>(cands[i]) = i;
    }
    CalculateWeight(id, target, cands, owner_.opts_.cands_opts_.overhang_weight);

    std::make_heap(cands.begin(), cands.end(), [](const std::tuple<const Overlap*, double, size_t>& a, const std::tuple<const Overlap*, double, size_t>& b) {
       return std::get<1>(a) < std::get<1>(b);    // CAUTION, calulated by CalculateWeight
    });
    size_t heap_size = cands.size();
    aligner_.SetTarget(target);
    std::vector<int> coverage(target.Size(), 0);

    std::vector<Alignment> first_als;
    std::vector<Alignment> flt_als;    // 
    int num_consecu_fail = 0;
    int accu_base = 0;
    while (heap_size > 0) {
        auto i = std::get<2>(cands[0]);
        auto ol = group.Get(i, 0);
        const auto& tread = ol->GetRead(id);
        const auto& qread = ol->GetOtherRead(id);
        Alignment al(tread.id, qread.id);
        bool r = false;
        double best_identity = 0.0;

        for (size_t j = 0; j < group.Size(i); ++j) {
            auto ol = group.Get(i, j);
            Alignment al_local(tread.id, qread.id);
        
            auto r_local = GetAlignment(id, ol, uc, al_local);
            DEBUG_printf("alignment: r = %d, q = (%zd %zd %zd),  d=%d, t = (%zd %zd %zd), d=%zd,%f\n", r_local,
               al_local.query_start, al_local.query_end, al_local.QuerySize(), ol->SameDirect(),
               al_local.target_start, al_local.target_end, al_local.TargetSize(), al_local.distance, al_local.Identity());

            if (r_local && !ExactFilter(al_local)) {
                if (best_identity < al_local.Identity()) {
                    best_identity = al_local.Identity();
                    r = r_local;
                    al = al_local;
                }
                if (j >= 3) break;
            }
        }

        if (r) accu_base += al.AlignSize();
        if (r && !ExactFilter(al)) { 
            stat_info.aligns[0]++;
            first_als.push_back(al);
            std::for_each(coverage.begin()+al.target_start, coverage.begin()+al.target_end, [](int& c) {c++;} );

            num_consecu_fail = 0;
        } else {
            stat_info.aligns[1]++;
            if (r) { flt_als.push_back(al); }
            num_consecu_fail++;
        }

        if (owner_.opts_.cands_opts_.IsEndCondition(coverage, first_als.size(), num_consecu_fail)) break;
        
        std::pop_heap(cands.begin(), cands.begin()+heap_size, [](std::tuple<const Overlap*, double, size_t>& a, std::tuple<const Overlap*, double, size_t>& b) {
            return std::get<1>(a) < std::get<1>(b);    // CAUTION, calulated by CalculateWeight
        });
        heap_size--;
    }

    if (first_als.size() > 0) {
    ////////

    if (owner_.opts_.check_local_identity_) {
        auto local_thresholds = CalculateLocalDistanceThreshold(first_als, 40, owner_.opts_.min_identity_);
        
        std::vector<Alignment> first_als1;
        for (size_t i = 0; i < first_als.size(); ++i) {
            if (CheckLocalDistance(first_als[i], local_thresholds)) {
                first_als1.push_back(first_als[i]);
                DEBUG_printf("ckck ADD\n");
            } else {
                
                DEBUG_printf("ckck 00\n");
            }
        }
        
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
            if (!ExactFilter(al, range)) {
               aligned_.push_back(al);
            }
        }

        for (const auto &al : flt_als) {
            if (!ExactFilter(al, range)) {
                aligned_.push_back(al);
            }
            if (owner_.opts_.cands_opts_.IsEnough(aligned_.size())) break;
        }

        for (auto &al : aligned_) { al.Rearrange(); }
        graph_.Build(target, range, aligned_);
        graph_.Consensus();
        return true;
    }
    }
    return false;
}

void ReadCorrect::Worker::CalculateWeight(Seq::Id id,  const DnaSeq& target, std::vector<std::tuple<const Overlap*, double, size_t>>& cands, double opt_ohwt) {
    std::vector<double> cand_cov_wts (target.Size()+1);

    double wtsum = 0.0;
    for (auto &it : cands) {
        auto o = std::get<0>(it); // it.first;
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        double ohwt = opt_ohwt * o->identity_ / 100;
        double olwt = o->identity_ / 100;

        auto mr = o->MappingTo<2>(t, {0, q.len});
        auto start = std::max(0, mr[0] < mr[1] ? mr[0] : mr[1]);
        auto end =   std::min(t.len, mr[0] >= mr[1] ? mr[0] : mr[1]);
        // start -- t.start -- t.end -- end
        assert(start <= t.start && t.end <= end);

        cand_cov_wts[start]   += ohwt;
        cand_cov_wts[t.start] += (olwt - ohwt);
        cand_cov_wts[t.end]   -= (olwt - ohwt);
        cand_cov_wts[end]     -= ohwt;

        wtsum += olwt;
    }

    for (size_t i=1; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] += cand_cov_wts[i-1];
    }
    assert(std::abs(cand_cov_wts.back()) < 0.0000001);  // cand_cov_wts.back() == 0

    for (size_t i=0; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] = wtsum - cand_cov_wts[i];
    }

    std::for_each(cands.begin(), cands.end(), [opt_ohwt, id, &cand_cov_wts, wtsum](std::tuple<const Overlap*, double, size_t>& it) {
        auto o = std::get<0>(it);//it.first;
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        double ohwt = opt_ohwt * o->identity_ / 100;
        double olwt = o->identity_ / 100;

        auto mr = o->MappingTo<2>(t, {0, q.len});
        auto start = std::max(0, mr[0] < mr[1] ? mr[0] : mr[1]);
        auto end =   std::min(t.len, mr[0] >= mr[1] ? mr[0] : mr[1]);
        // start -- t.start -- t.end -- end
        assert(start <= t.start && t.end <= end);

        double wt = std::accumulate(cand_cov_wts.begin()+start, cand_cov_wts.begin()+t.start, 0.0) * ohwt +
                    std::accumulate(cand_cov_wts.begin()+t.start, cand_cov_wts.begin()+t.end, 0.0) * olwt + 
                    std::accumulate(cand_cov_wts.begin()+t.end, cand_cov_wts.begin()+end, 0.0) * ohwt;
        
        //it.second = wt;
        std::get<1>(it) = wt;
    });
}

} // namespace fsa {
    