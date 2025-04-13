#include "asm_dataset.hpp"

#include <algorithm>
#include <memory>
#include <random>

#include "../utility.hpp"
#include "../kmer/kmer.hpp"

#include "edlib.h"

namespace fsa {
    
void AsmDataset::Load() {
    assert(!opts_.read_file.empty());
    rd_store_.Load(opts_.read_file, "", true);
    
    LoadOverlaps(opts_.ifname);
}


void AsmDataset::Purge() {
    
    GroupOverlaps();
    
    //FilterOverhang();

    if (is_ol_accurate) FilterLowQuality();

    EstimateCoverage();

    ExtendOverlapToEnd();

    FilterCoverage();

    FilterDuplicate();

    EstimateGenomeSize();

    FilterContained();
    
    Dump();
}

void AsmDataset::LoadOverlaps(const std::string &fname) {
    LOG(INFO)("Load overlap file");

    std::atomic<size_t> count {0};
    auto check = [&](Overlap& o) {
        count++;
        return opts_.filter0.Valid(o);
    };

    assert(rd_store_.Size() > 0);
    ol_store_.LoadFast(fname, "", (size_t)opts_.thread_size, check);

    LOG(INFO)("Overlap size: %zd/%zd", ol_store_.Size(), count.load());

    for (size_t i=0; i<rd_store_.Size(); ++i) {
        read_infos_[i] = ReadStatInfo();
        read_infos_[i].len = rd_store_.GetSeqLength(i);
        read_infos_[i].id = i; // TODO
    }

    TestOverlapIdentity();

}


double AsmDataset::CalcLocalOverhangThreshold(std::vector<std::array<double,2>> &overhangs) {
        
    std::sort(overhangs.begin(), overhangs.end(), [](const std::array<double,2>& a, const std::array<double,2> &b){
        return  a[0]*1.0/ a[1] < b[0] *1.0/ b[1];
    });
    double overhang_threshold = opts_.filter0.max_overhang;

    
    double median = 0;
    double mad = 0;
    size_t size = std::min((int)overhangs.size(), 40);  // TODO
    std::vector<std::array<double,2>> ohs(size);
    std::transform(overhangs.begin(), overhangs.begin()+size, ohs.begin(),  [](const std::array<double,2> &a) {
        return std::array<double,2>({a[0], a[1]});
    });

    if (ohs.size() > 0) {
        ComputeMedianAbsoluteDeviation(ohs, median, mad);
        if (mad >= 1.0e-10) {
            overhang_threshold = median + 6*1.4826*mad;
            overhang_threshold = std::min<double>(overhang_threshold, opts_.filter0.max_overhang);

        } else {
            ComputeMeanAbsoluteDeviation(ohs, median, mad);
            overhang_threshold = median + 6*1.253*mad;
            overhang_threshold = std::min<double>(overhang_threshold, opts_.filter0.max_overhang);
        }
    }
    return overhang_threshold;
}

void AsmDataset::FilterLowQuality() {
    LOG(INFO)("Filter low-quality overlaps");

    std::mutex mutex;
    std::unordered_set<const Overlap*> exist;
    auto combine_func = [this, &mutex, &exist](std::vector<const Overlap*> &flt) {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto o : flt) {
            if (exist.find(o) != exist.end()) {
                SetOlReason(*o, OlReason::Simple());
            } else {
                exist.insert(o);
            }
        }
    };


    std::atomic<size_t> index {0};
    auto work_func = [this, &index, combine_func](size_t tid) {
        std::vector<const Overlap*> flt;
        for (size_t i = index.fetch_add(1); i < rd_store_.Size(); i = index.fetch_add(1)) {
            FilterLowQuality(i, grouper_.Get(i), flt);
        }
        combine_func(flt);
    };

    MultiThreadRun(opts_.thread_size, work_func);
}

void AsmDataset::FilterLowQuality(int id, const OverlapGrouper::Group &group, std::vector<const Overlap*>& ignored) {

    const int WIN_SIZE = 4000;      // TODO param: 
    const int MIN_COV = 20;         // TODO param;

    auto& rinfo = read_infos_[id];
    assert(rinfo.len > 0);
    const size_t win_count = std::max(1, (rinfo.len + WIN_SIZE / 2) / WIN_SIZE);
    const size_t win_size = (rinfo.len + win_count - 1 ) /  win_count;

    std::vector<std::vector<std::array<double,2>>> winidents(win_count);

    // collect target read information: 
    for (size_t i = 0; i < group.Size(); ++i) {
        for (size_t j = 0; j < group.Size(i); ++j) {
            auto &ol = *group.Get(i, j);

            if (IsReserved(ol)) {
                auto &tr = ol.GetRead(id);

                size_t s = (tr.start + win_size / 2 ) / win_size;
                size_t e = (tr.end + win_size / 2 ) / win_size;
                bool end0 = s * win_size >= tr.start;
                bool end1 = e * win_size <= tr.end;
                assert(e >= s && s >= 0 && e <= winidents.size());

                for (size_t iw = s; iw < e; ++iw) {
                    assert(iw >= 0 && iw < winidents.size());
                    if (iw == s && end0) {
                        winidents[iw].push_back({ ol.Identity(), 0.5});
                    } else if (iw + 1 == e && end1) {
                        winidents[iw].push_back({ ol.Identity(), 0.5});
                    } else {
                        winidents[iw].push_back({ ol.Identity(), 1 + 1.0*(e-s) / win_count});
                    }
                } 

                // std::for_each(winidents.begin()+s, winidents.begin()+e, [&ol, s, e](std::vector<std::array<double,2>>& v) {
                //     v.push_back({ ol.Identity(), double(e-s)});
                // });  
                
                if (rd_store_.QueryNameById(id) == opts_.debug_name) {
                    LOG(INFO)("add idt: %s: (%zd, %zd) %.04f", rd_store_.QueryNameById(ol.GetOtherRead(id).id).c_str(), s, e, ol.Identity());
                }
            }
        }
    }

    std::vector<double> identity_threshold(winidents.size());
    std::transform(winidents.begin(), winidents.end(), identity_threshold.begin(), [this, id](std::vector<std::array<double,2>>& ident) {
        std::sort(ident.begin(), ident.end(), [](std::array<double,2>& a, std::array<double,2>& b) { return a[0] > b[0]; });
        if (ident.size() >= MIN_COV) {

            double median, mad;
            ComputeMedianAbsoluteDeviation(std::vector<std::array<double,2>>(ident.begin(), ident.begin()+MIN_COV),  median, mad);
            if (rd_store_.QueryNameById(id) == opts_.debug_name) {
                LOG(INFO)("filter_low_quality: size=%zd %0.04f, %0.04f, %0.04f",ident.size(), median, mad, std::max(opts_.filter0.min_identity, median-6*1.4826*mad));
                for (size_t i = 0; i <ident.size(); ++i) {
                    LOG(INFO)("detial: %zd = %.04f", i, ident[i][0]);
                }
            }
            return std::max(opts_.filter0.min_identity /100, median-6*1.4826*mad);

        } else {
            return ident.size() > 0 ? std::max(opts_.filter0.min_identity /100, ident.back()[0]) : opts_.filter0.min_identity /100;
        }
    });

    for (auto& th : identity_threshold) {
        th = std::floor(th * 100) / 100;
    }

    rinfo.identity_threshold = identity_threshold;

    auto check_ol_identity = [win_size, this, id](size_t start, size_t end, double ident, const std::vector<double>& thresholds) {

        size_t s = (start + win_size / 2) / win_size;
        size_t e = (end + win_size / 2) / win_size;
        assert (s >= 0 && s <= e && e <= thresholds.size());

        return e > s && ident >= std::accumulate(thresholds.begin()+s, thresholds.begin()+e, 0.0) / (e - s);
    };

    for (size_t i = 0; i < group.Size(); ++i) {
        for (size_t j = 0; j < group.Size(i); ++j) {
            auto &ol = *group.Get(i, j);

            if (IsReserved(ol)) {
                auto &tr = ol.GetRead(id);
                
                if (!check_ol_identity(tr.start, tr.end, ol.Identity(), identity_threshold)) {
                    ignored.push_back(&ol);
                }

            }
        }
    }
    
}


void AsmDataset::FilterOverhang() {
    LOG(INFO)("Filter overlaps with long overhang");

    std::mutex mutex;
    auto combine_func = [this, &mutex](std::unordered_set<const Overlap*> &flt) {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto o : flt) {
            SetOlReason(*o, OlReason::Simple());
        }
    };

    std::atomic<size_t> index {0};
    auto work_func = [this, &index, combine_func](size_t tid) {
        std::unordered_set<const Overlap*> flt;
        for (size_t i = index.fetch_add(1); i < rd_store_.Size(); i = index.fetch_add(1)) {
            FilterOverhang(i, grouper_.Get(i), flt);
        }
        combine_func(flt);
    };

    MultiThreadRun(opts_.thread_size, work_func);
}

void AsmDataset::FilterOverhang(Seq::Id id, const OverlapGrouper::Group& group, std::unordered_set<const Overlap*>& ignored) {
    const int MIN_COV = 30;         // param;

    auto& rinfo = read_infos_[id];

    std::vector<double> lohs;
    std::vector<double> rohs;

    // collect overhang information: 
    for (size_t i = 0; i < group.Size(); ++i) {
        for (size_t j = 0; j < group.Size(i); ++j) {
            auto &ol = *group.Get(i, j);

            if (IsReserved(ol)) {
                auto &tr = ol.GetRead(id);

                auto oh = ol.Overhang2();
                if (tr.id == ol.a_.id) {    
                    if (oh[0] & 0x1) {
                        lohs.push_back(tr.start - 0);
                    }
                    if (oh[0] & 0x2) {
                        rohs.push_back(tr.len - tr.end);
                    }
                } else {
                    if (oh[1] & 0x1) {
                        lohs.push_back(tr.start - 0);
                    }
                    if (oh[1] & 0x2) {
                        rohs.push_back(tr.len - tr.end);
                    }
                }
            }
        }
    }


    auto calc_oh_threshold = [this](std::vector<double>& ohs) -> double {
        if (ohs.size() > MIN_COV) {
            double median, mad;
            std::sort(ohs.begin(), ohs.end());
            ComputeMedianAbsoluteDeviation(ohs, median, mad);
            return std::min<double>(median + 6*1.4826*mad, opts_.filter0.max_overhang);
        } else {
            return opts_.filter0.max_overhang;
        }
    };
    rinfo.overhang_l_threshold = calc_oh_threshold(lohs);
    rinfo.overhang_r_threshold = calc_oh_threshold(rohs);
}

void AsmDataset::DumpOverlaps(const std::string &fname) const {
    
    auto filter = [&](const Overlap& o)->bool {
        return GetOlReason(o).type == OlReason::RS_OK;
    };

    ol_store_.Save(fname, "", opts_.thread_size, filter);
}


void AsmDataset::ExtendOverlapToEnd(const Overlap &oldone, int maxoh) {
    if (oldone.Location(0) != Overlap::Loc::Abnormal) return;
    if (oldone.Location(maxoh) == Overlap::Loc::Abnormal) return;

    Overlap& o = const_cast<Overlap&>(oldone);
    size_t alsize = o.AlignedSize();
    if (o.a_.strand == o.b_.strand) {
        if (o.a_.start <= maxoh && o.b_.start <= maxoh) {
            o.a_.start = 0;
            o.b_.start = 0;
        }
        else if (o.a_.start <= maxoh) {
            o.b_.start -= o.a_.start;
            o.a_.start = 0;
        } 
        else if (o.b_.start <= maxoh) {
            o.a_.start -= o.b_.start;
            o.b_.start = 0;
        }

        if (o.a_.end >= o.a_.len - maxoh && o.b_.end >= o.b_.len - maxoh) {
            o.a_.end = o.a_.len;
            o.b_.end = o.b_.len;
        }
        else if (o.a_.end >= o.a_.len - maxoh) {
            o.b_.end += o.a_.len - o.a_.end;
            o.a_.end = o.a_.len;
        }
        else if (o.b_.end >= o.b_.len - maxoh) {
            o.a_.end += o.b_.len - o.b_.end;
            o.b_.end = o.b_.len;
        }

    }
    else {
        if (o.a_.start <= maxoh && o.b_.end >= o.b_.len - maxoh) {
            o.a_.start = 0;
            o.b_.end = o.b_.len;
        }
        else if (o.a_.start <= maxoh) {
            o.b_.end += o.a_.start;
            o.a_.start = 0;
        }
        else if (o.b_.end >= o.b_.len - maxoh) {
            o.a_.start -= o.b_.len - o.b_.end;
            o.b_.end = o.b_.len;
        }

        if (o.b_.start <= maxoh && o.a_.end >= o.a_.len - maxoh) {
            o.b_.start = 0;
            o.a_.end = o.a_.len;
        }
        else if (o.b_.start <= maxoh) {
            o.a_.end += o.b_.start;
            o.b_.start = 0;

        }
        else if (o.a_.end >= o.a_.len - maxoh) {
            o.b_.start -= o.a_.len - o.a_.end;
            o.a_.end = o.a_.len;
        }
    }
    o.UpdateIdentity(alsize);

    assert(o.Location(0) != Overlap::Loc::Abnormal);

}

void AsmDataset::FilterContained() {
    LOG(INFO)("Remove contained reads");

    std::atomic<size_t> index { 0 };

    auto set_contained = [this](int contained, int containing) {
        auto iter = read_infos_.find(contained);
        assert(iter != read_infos_.end());
        iter->second.filtered = RdReason::Contained(containing);
    };
    
    auto work_func = [&](size_t tid) {
        std::unordered_set<const Overlap*> ignored;
        
        for (size_t i = index.fetch_add(1); i < ol_store_.Size(); i = index.fetch_add(1)) {
            const Overlap& o = ol_store_.Get(i);

            if (IsReserved(o)) {
                auto loc = o.Location(0);
                Seq::Id contained = Seq::NID;
                Seq::Id containing = Seq::NID;
                if (loc == Overlap::Loc::Equal) {
                    contained = std::max(o.a_.id, o.b_.id);
                    containing = std::min(o.a_.id, o.b_.id);
                } else if (loc == Overlap::Loc::Contained) {
                    contained = o.a_.id;
                    containing = o.b_.id;
                } else if (loc == Overlap::Loc::Containing) {
                    contained = o.b_.id;
                    containing = o.a_.id;
                }
                if (containing != Seq::NID && contained != Seq::NID) {
                    // auto contained_others = GetOverlapReads(contained);
                    // auto containing_others = GetOverlapReads(containing);

                    bool contain_all = true;
                    // for (auto i : contained_others) {
                    //     if (containing_others.find(i) == containing_others.end()) {
                    //         contain_all = false;
                    //         break;
                    //     }
                    // }

                    if (contain_all) {
                        set_contained(contained, containing);
                    } else {
                        ignored.insert(&o);
                    }

                }
                

            }
            for (auto o : ignored) {
                SetOlReason(*o, OlReason::Simple());
            }
        }
    };

    MultiThreadRun((size_t)opts_.thread_size, work_func);

    LOG(INFO)("Start filtering contained reads and relative overlaps");
    UpdateFilteredRead();
}


void AsmDataset::FilterCoverage() {
    LOG(INFO)("Check Coverage");

    std::atomic<size_t> index { 0 };
    auto work_func = [this, &index](size_t tid) {
        for (size_t i = index.fetch_add(1); i < rd_store_.Size(); i = index.fetch_add(1)) {
            auto rinfo = read_infos_.find(i);
            auto group = grouper_.Get(i);
            if (rinfo != read_infos_.end() && !group.Empty()) {
                AnalyzeCoverage(i, group, rinfo->second);
            }
        }
    };

    MultiThreadRun(opts_.thread_size, work_func);  

    auto avg_cov = EstimateCoverage(read_infos_);
    LOG(INFO)("The avarage coverage is %.02f", avg_cov);
    
    // std::unordered_set<Seq::Id> done;
    // for (const auto &i : read_infos_) {
    //     if (done.find(i.first) != done.end()) continue;
    //     if (i.second.minmax_coverage[1] >= avg_cov * 4 ) {
    //         std::unordered_set<Seq::Id> group {i.first};
            
    //         std::vector<Seq::Id> check_list { i.first };
    //         size_t index = 0;
    //         while (index < check_list.size()) {
    //             auto id = check_list[index];
    //             index++;
    //             auto ols = grouper_.Get(id);
    //             if (!ols.Empty()) {
    //                 for (size_t i = 0; i < ols.Size(); ++i) {
    //                     auto ol = ols.Get(i, 0);
    //                     auto oid = ol->GetOtherRead(id).id;
    //                     if (group.find(oid) == group.end()) {
    //                         if (read_infos_[oid].minmax_coverage[1] >= 90) {
    //                             check_list.push_back(oid);
    //                             group.insert(oid);
    //                         }
    //                     }
    //                 }
    //             }
    //         }

    //         //LOG(INFO)("Group: %zd", group.size());
    //         for (auto i : group) {
    //             DUMPER["data"]("g %zd %s %d\n", group.size(), rd_store_.QueryNameById(i).c_str(), read_infos_[i].minmax_coverage[1]);
    //         }
    //         //ClusterBundle(group);
            
    //         done.insert(group.begin(), group.end());

    //     }
    //}

    auto threshold = CalcCoverageThreshold();
    int mincov = opts_.min_coverage < 0 ? threshold[0] : opts_.min_coverage;
    int maxcov = threshold[1];
    int maxdiff = threshold[2];
    LOG(INFO)("min_coverage = %d(%d), max_coverage = %d, max_diff_coverage = %d", mincov, threshold[0], maxcov, maxdiff);
    for (auto &c : read_infos_) {
        const auto &minmax = c.second.minmax_coverage;
        if (minmax[0] < mincov || minmax[1] > maxcov || minmax[1]-minmax[1] > maxdiff) {
            c.second.filtered = RdReason::Coverage(minmax);
        }

        if (c.second.covtype != 0) {
            c.second.filtered = RdReason::CoverageType(c.second.covtype);
        }
    }

    UpdateFilteredRead();

}


double AsmDataset::GetOverlapQuality0(const Overlap &ol) {
    const auto &rd_store = rd_store_;

    const auto & query = rd_store.GetSeq(ol.a_.id);
    const auto & target = rd_store.GetSeq(ol.b_.id);

    auto tseq = target.ToUInt8(ol.b_.start, ol.b_.end);
    auto qseq = query.ToUInt8(ol.a_.start, ol.a_.end, !ol.SameDirect());

    auto r = edlibAlign((const char*)&qseq[0], qseq.size(), (const char*)&tseq[0], tseq.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    if (r.status == EDLIB_STATUS_OK) {
        return 1.0 - r.editDistance * 1.0 / ol.AlignedLength();
    } else {
        return 0.0;
    }
}

void AsmDataset::GroupOverlaps() {
    LOG(INFO)("Group overlaps");
    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<Seq::Id>(), [](const Overlap &a, const Overlap &b) {
        if (a.AlignedSize() > b.AlignedSize()) return 1;
        else if (a.AlignedSize() == b.AlignedSize()) return 0;
        else return -1;
    }); 
}

void AsmDataset::FilterDuplicate() {
    LOG(INFO)("Group overlaps and remove duplicated");

    auto better = [](const Overlap* a, const Overlap *b) {
        assert(a != nullptr && b != nullptr);
        return a->Identity() > b->Identity() || 
               a->Identity() ==  b->Identity() && a->AlignedSize() > b->AlignedSize();
    };

    for (size_t rid = 0; rid < rd_store_.Size(); ++rid) {
        auto group = grouper_.Get(rid);
        for (size_t i = 0; i < group.Size(); ++i) {
            const Overlap* best = nullptr;
            for (size_t j = 0; j < group.Size(i); ++j) {
                const Overlap &ol = *group.Get(i,j);
                if (IsReserved(ol)) {
                    if (best != nullptr) {
                        if (better(&ol, best)) {
                            SetOlReason(*best, OlReason::Duplicate());
                            best = &ol;
                        } else {
                            SetOlReason(ol, OlReason::Duplicate());
                        }

                    } else {
                        best = &ol;
                    }
                }
            }
        }
    }

    LOG(INFO)("Overlap size: %zd/%zd", ReservedSize(), ol_store_.Size());
}
        

void AsmDataset::ExtendOverlapToEnd() {
    LOG(INFO)("Extend Overlaps to ends");

    std::atomic<size_t> index { 0 };
    auto work_func = [this, &index](size_t tid) {
        for (size_t i = index.fetch_add(1); i <  ol_store_.Size(); i = index.fetch_add(1)) {
            const Overlap& o = ol_store_.Get(i);
            if (IsReserved(o) || GetOlReason(o).type == OlReason::RS_DUPLICATE) {
                auto oh = o.Overhang();
                int th = std::max(oh[0], oh[1]);
                if (th > 0) {
                    ExtendOverlapToEnd(o, th);
                }
            }
        }
    };

    MultiThreadRun(opts_.thread_size, work_func);
}


void AsmDataset::AnalyzeCoverage(int id, const OverlapGrouper::Group& group, ReadStatInfo &rinfo) {

    if (group.Size() > 0) {
        std::vector<int> cov(group.Get(0, 0)->GetRead(id).len + 1, 0);
        const int redundance = - std::min<int>(500, cov.size()/10);

        for (size_t i = 0; i < group.Size(); ++i) {
            const Overlap& o = *group.Get(i, 0);

            if (IsReserved(o)) {

                for (size_t j = 0; j < group.Size(i); ++j) {
                    const Overlap& ol = *group.Get(i, j);
                    auto& r = ol.GetRead(id);
                    if (std::max<size_t>(0, r.start-redundance) < std::min<size_t>(cov.size()-1, r.end+redundance)) {

                        cov[std::max<size_t>(0, r.start-redundance)] ++;
                        cov[std::min<size_t>(cov.size()-1, r.end+redundance)] --;
                    }

                }
            }
        }

        // 展开coverage
        for (size_t i = 1; i < cov.size(); ++i) {
            cov[i] += cov[i - 1];
        }
        assert(cov.back() == 0);

        if ( opts_.debug_name == rd_store_.QueryNameById(id)) {
           for (size_t i = 0; i < cov.size(); ++i) {
               printf("cov %d\n", cov[i]);
           }
        }

        // TODO 根据Coverage分布检查序列是否异常

        // 识别Coverage剧烈变换的位置
        if (cov.size() > 3000) {
            rinfo.cliff = CoverageConfidencePoints1(std::vector<int>(cov.begin()-redundance, cov.end()+redundance), rd_store_.QueryNameById(id) == opts_.debug_name);
            if (rinfo.cliff[0] > 0) rinfo.cliff[0] -= redundance;
            if (rinfo.cliff[1] > 0) rinfo.cliff[1] -= redundance;

            if (rinfo.cliff[0] > 0) {
                if (rinfo.cliff[0] >= std::min<int>(opts_.max_unreliable_length, opts_.max_unreliable_rate*cov.size())) {
                    rinfo.cliff[0] = -1;
                }
            }
            
            if (rinfo.cliff[1] > 0) {
                if ((int)cov.size() - rinfo.cliff[1] >= std::min<int>(opts_.max_unreliable_length, opts_.max_unreliable_rate*cov.size())) {
                    rinfo.cliff[1] = -1;
                }
            }
        }

        // 计算序列的左中右三个部分的平均覆盖度
        int inv = cov.size() / 3;
        auto a0 = std::accumulate(cov.begin(), cov.begin()+inv, 0) / inv;
        auto a1 = std::accumulate(cov.begin()+inv, cov.begin()+2*inv, 0) / inv;
        auto a2 = std::accumulate(cov.begin()+2*inv, cov.end(), 0) / inv;
        rinfo.coverage = {a0, a1, a2};

        // 
        // 读数的最大覆盖度和最小覆盖度，两端不计算
        int oh = std::max(0, -redundance) + opts_.filter0.max_overhang;
        auto c_minmax = std::minmax_element(oh + cov.begin(), cov.end() -1 - oh);
        rinfo.minmax_coverage = {*c_minmax.first, *c_minmax.second };
        rinfo.covtype = 0;

    } 
}

int AsmDataset::AnalyzeCoverageType(const std::vector<int>& cov, bool log) {
    const int INV = 100;
    assert(cov.size() > INV);

    std::vector<double> smooth(cov.size()-INV+1);
    smooth[0] = std::accumulate(cov.begin(), cov.begin()+INV, 0.0);
    for (size_t i = 1; i < smooth.size(); ++i) {
        smooth[i] = smooth[i-1] + cov[i+INV-1] - cov[i-1];
    }

    for (size_t i = 0; i < smooth.size(); ++i) {
        smooth[i] =  smooth[i] / INV;
    }

    auto minmax = std::minmax_element(smooth.begin(), smooth.end());
    auto diff = *minmax.second - *minmax.first;
    if (diff >= 10 && diff > *minmax.second / 3) {
        std::array<double, 4> values;
        std::array<decltype(minmax.first), 4> positions;
        if (minmax.first < minmax.second) {
            values[1] = *minmax.first;
            positions[1] = minmax.first;
            values[2] = *minmax.second;
            positions[2] = minmax.second;
        } else {
            values[1] = *minmax.second;
            positions[1] = minmax.second;
            values[2] = *minmax.first;
            positions[2] = minmax.first;
        }

        if (values[1] < values[2]) {
            auto ml = std::max_element(smooth.begin(), positions[1]);
            values[0] = *ml;
            positions[0] = ml;

            auto mr = std::min_element(positions[2],smooth.end());
            values[3] = *mr;
            positions[3] = mr;
        } else {
            auto m = std::min_element(smooth.begin(), positions[1]);
            values[0] = *m;
            positions[0] = m;

            auto mr = std::max_element(positions[2],smooth.end());
            values[3] = *mr;
            positions[3] = mr;
        }
        if (values[0] < values[1]) {
            auto diff = values[1] - std::max(values[0], values[2]);
            if (diff > 10 && diff > values[1] * 2/3) {
                return 0;
            }
        } else {

            auto diff = std::min(values[0], values[2]) - values[0];
            if (diff > 10 && diff > std::min(values[0], values[2]) * 2/3) {
                return 2;
            }
        }

        
        if (values[3] < values[2]) {
            auto diff = values[2] - std::max(values[1], values[3]);
            if (diff > 10 && diff > values[2] * 2/3) {
                return 0;
            }
        } else {

            auto diff = std::min(values[1], values[3]) - values[2];
            if (diff > 10 && diff > std::min(values[1], values[3]) * 2/3) {
                return 2;
            }
        }
    }
    return 0;
}

void AsmDataset::CoverageConfidencePoints(const std::vector<int>& cov, bool log) {
    const int INV = 500;
    assert(cov.size() > INV);

    auto is_cliff = [](double a, double b) {
        auto diff = std::abs(b - a);
        return diff >= 10 && diff > std::max(a,b) / 3;
    };

    std::vector<double> smooth(cov.size()-INV+1);
    smooth[0] = std::accumulate(cov.begin(), cov.begin()+INV, 0.0);
    for (size_t i = 1; i < smooth.size(); ++i) {
        smooth[i] = smooth[i-1] + cov[i+INV-1] - cov[i-1];
    }

    for (size_t i = 0; i < smooth.size(); ++i) {
        smooth[i] =  smooth[i] / INV;
    }

    auto analyse_slope = [](double a, double b, int step) {
        double rate = (b - a) * 2 / (a + b);
        if (rate < -0.2) {
            return -1;
        } else if (rate > 0.2) {
            return 1;
        } else {
            return 0;
        }
    };
    
    const int stepsize = 100;
    int state = 0; int start = 0;
    std::vector<std::array<int,2>> states;
    for (size_t i=0; i+stepsize < smooth.size(); i+=stepsize) {
        int s = analyse_slope(smooth[i], smooth[i+stepsize], stepsize);

        if (s != state) {
            if (state != 0) {
                if (is_cliff(smooth[start], smooth[i])) {
                    states.push_back({start, (int)i});
                    if (log) {
                        printf("add point %d(%f) , %zd(%f) \n", start, smooth[start], i, smooth[i]);

                    }
                }
                start = i;
            }
            if (log) {
                printf("xxx add point %d(%f) , %zd(%f) \n", start, smooth[start], i, smooth[i]);

            }
            state = s;
        } else {
            if (s == 0) {
                start = i;
            }
        }
    }

}

template<typename O, typename I>
std::vector<O> SmoothLine(const std::vector<I> &line, size_t winsize, bool log) {
    assert(line.size() >= winsize);
    std::vector<O> smooth(line.size() - winsize + 1);
    smooth[0] = std::accumulate(line.begin(), line.begin()+winsize, 0);
    for (size_t i = 1; i < smooth.size(); ++i) {
        smooth[i] = smooth[i-1] + line[i+winsize-1] - line[i-1];
    }
    for (size_t i = 0; i < smooth.size(); ++i) {
        if (log)
            printf("smooth %zd: %f, %d\n", i, smooth[i], line[i]);
        smooth[i] = smooth[i] / winsize;
    }
    return smooth;
}

template<typename T>
std::vector<std::array<int,2>> FindCorners(const std::vector<T>& line, bool log) {
    const size_t stepsize = 500;
    assert(line.size() > stepsize);

    const double slope_threshold = 0.001;

    auto is_cliff = [](double a, double b, size_t dis, double max_slope) {
        auto diff = std::abs(b - a);
        auto slope = diff / dis;
        
        return diff >= 15 && diff >= std::min(a,b) &&  (diff >= 2*std::min(a,b) ||  std::abs(slope) > 0.01);
    };

    auto slope_state = [slope_threshold](double s) {
        const double S = slope_threshold;
        return s < -S ? -1 : (s > S ? 1 : 0);
    };
    std::vector<std::array<int,2>> cliffs;
    auto max_slope = (line[stepsize] - line[0]) / stepsize;
    auto last = slope_state(max_slope);
    auto s0_count = 0;
    int max_s0_count = 3;
    size_t index = 0;
    for (size_t i = stepsize; i + stepsize < line.size(); i += stepsize) {
        auto slope =  (line[i+stepsize] - line[i]) / stepsize;
        auto curr = slope_state(slope);

        if (log) {
            printf("slope %zd - %d, %f %f\n", i, curr, slope, line[i]);
        }

        if (std::abs(max_slope) < std::abs(slope)) {
            max_slope = slope;
        }

        s0_count = curr == 0 ? s0_count + 1 : 0;

        if (last != curr && (curr != 0 || s0_count >= max_s0_count)) {
            size_t end = i;//i - s0_count*stepsize;
            // if (log) {
            //     printf("slope test %d, %d, %d  %f %zd\n", last != 0, (curr == 0 && s0_count >= max_s0_count || curr != 0) , is_cliff(line[i], line[index], i-index, max_slope), line[index], i-index);
            // }
            if (last != 0 && ((curr == 0 && s0_count >= max_s0_count) || curr != 0) && is_cliff(line[end], line[index], end-index, max_slope)) {
                cliffs.push_back({(int)index, (int)end});
                // if (log) {
                //     printf("cliff %zd - %zd, %f %f\n", end, line[end],  index, line[index]);
                // }
            }
            last = curr;
            index = i;
            max_slope = slope;
        }
    }
    size_t end = line.size() - 1 ;//- s0_count * stepsize;
    if (last != 0 && is_cliff(line[end], line[index], end-index, max_slope)) {
        cliffs.push_back({(int)index, (int)end});
        // if (log) {
        //     printf("cliff %zd - %zd, %f %f\n", end, line[end],  index, line[index]);
        // }
    }


    return cliffs;
}


std::array<int,2> AsmDataset::CoverageConfidencePoints1(const std::vector<int>& cov, bool log) {
    const size_t INV = 500;
    assert(cov.size() > INV);

    std::vector<double> smooth = SmoothLine<double>(cov, INV, log);
    auto corners = FindCorners(smooth, log);

    std::array<int,2> cliff { -1, -1};

    if (log) {
        for (auto c : corners) {
            printf("corners: %d:%f, %d:%f\n", c[0], smooth[c[0]], c[1], smooth[c[1]]);
        }
    
    }

    for (auto c : corners) {
        if (smooth[c[0]] > smooth[c[1]]) {
            cliff[0] = c[1] + 500;
            // 不需要 break 需要找到最后一个 
        }
    }
    
    for (auto c : corners) {
        if (smooth[c[0]] < smooth[c[1]]) {
            cliff[1] = c[0] + 500;
            break; // 找到第一个
        }
    }


    return cliff;

}


std::array<int,3> AsmDataset::CalcCoverageThreshold() const {
    assert(read_infos_.size() > 0);

    std::vector<int> cov_min;
    cov_min.reserve(read_infos_.size());
    std::vector<int> cov_max;
    cov_max.reserve(read_infos_.size());
    std::vector<int> cov_diff;
    cov_diff.reserve(read_infos_.size());
    for (const auto & i : read_infos_) {
        const auto& minmax = i.second.minmax_coverage;
        if (minmax[0] >= 0) {
            cov_min.push_back(minmax[0]);
            cov_max.push_back(minmax[1]);
            cov_diff.push_back(minmax[1]-minmax[0]);
        }
    }
    return {FirstTrough(cov_min, 100, 9), Percentile(cov_max, 100-0.01), Percentile(cov_diff, 100-0.01) };
}


bool AsmDataset::IsContained(const Overlap& o, std::array<int, 2> &rel) {
    auto loc = o.Location(0);
    if (loc == Overlap::Loc::Contained || loc == Overlap::Loc::Containing || loc == Overlap::Loc::Equal) {
        //  contained = rel[0], contain = rel[1]

        if (loc == Overlap::Loc::Equal) {
            rel[0] = std::max(o.a_.id, o.b_.id);
            rel[1] = std::min(o.a_.id, o.b_.id);
        }
        else if (loc == Overlap::Loc::Contained) {
            rel[0] = o.a_.id;
            rel[1] = o.b_.id;
        }
        else {
            assert(loc == Overlap::Loc::Containing);
            rel[0] = o.b_.id;
            rel[1] = o.a_.id;
        }
        return true;

    } else {
        return false;
    }

}

void AsmDataset::UpdateFilteredRead() {

    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);
        if (IsReserved(o)) {
            auto it = read_infos_.find(o.a_.id);
            assert(it != read_infos_.end());
            if (!it->second.filtered.IsOk()) {
                SetOlReason(o, OlReason::FilteredRead(it->first));
            } else {
                it = read_infos_.find(o.b_.id);                
                assert(it != read_infos_.end());
                if (!it->second.filtered.IsOk()) {
                    SetOlReason(o, OlReason::FilteredRead(it->first));
                }

            }
        }
    }
}


std::unordered_set<Seq::Id> AsmDataset::GetNearbyReads(Seq::Id tid) {
    std::unordered_set<Seq::Id> nearby;

    auto group = grouper_.Get(tid);
    if (!group.Empty()) {
        for (size_t i = 0; i < group.Size(); ++i) {
            
            auto &o = *group.Get(i, 0);
            auto &qread = o.GetOtherRead(tid);
            auto oltype = GetOlReason(o).type;
            if (oltype == OlReason::RS_OK) {
                nearby.insert(qread.id);
            } else if (oltype == OlReason::RS_FILTERED_READ) {
                auto ri = read_infos_.find(qread.id);
                if (ri->second.filtered.type == RdReason::RS_CONTAINED) {
                    nearby.insert(qread.id);
                }
            }
        }
    }
    return nearby;
}

std::unordered_set<Seq::Id> AsmDataset::GetOverlapReads(Seq::Id tid) const {
    std::unordered_set<Seq::Id> nearby;

    auto group = grouper_.Get(tid);
    if (!group.Empty()) {
        for (size_t i = 0; i < group.Size(); ++i) {
            auto &o = *group.Get(i, 0);
            if (IsReserved(o)) {
                auto &qread = o.GetOtherRead(tid);
                nearby.insert(qread.id);
            }
        }
    }
    return nearby;
}

std::unordered_set<const Overlap*> AsmDataset::GetExtendOverlaps(Seq::Id tid, int end) const {
    std::unordered_set<const Overlap*> extend;

    auto group = grouper_.Get(tid);
    if (!group.Empty()) {
        for (size_t i = 0; i < group.Size(); ++i) {
            auto &o = *group.Get(i, 0);
            auto &qread = o.GetOtherRead(tid);

            if (!( (end == 0 && o.Location(qread.id, 0) == Overlap::Loc::Left) ||
                   (end == 1 && o.Location(qread.id, 0) == Overlap::Loc::Right) )) continue;

            auto oltype = GetOlReason(o).type;
            if (oltype == OlReason::RS_OK) {
                extend.insert(&o);
            } else if (oltype == OlReason::RS_FILTERED_READ) {
                                
                auto qri = read_infos_.find(qread.id);
                auto tri = read_infos_.find(tid); 
                if (qri->second.filtered.type == RdReason::RS_CONTAINED || tri->second.filtered.type == RdReason::RS_CONTAINED) {
                    extend.insert(&o);
                }
            }
        }
    }

    return extend;
}


std::unordered_set<const Overlap*> AsmDataset::GetExtendOverlapsEx(Seq::Id tid, int end) const {
    std::unordered_set<const Overlap*> extend;

    DUMPER["test"]("extend: %s %d\n", string_pool_.QueryStringById(tid).c_str(), end);
    auto group = grouper_.Get(tid);
    if (!group.Empty()) {
        for (size_t i = 0; i < group.Size(); ++i) {
            auto &o = *group.Get(i, 0);
            auto &qread = o.GetOtherRead(tid);
            
            DUMPER["test"]("extend_checkt: %s %s\n", string_pool_.QueryStringById(tid).c_str(), string_pool_.QueryStringById(qread.id).c_str());

            if (!( (end == 0 && o.Location(qread.id, 0) == Overlap::Loc::Left) ||
                   (end == 1 && o.Location(qread.id, 0) == Overlap::Loc::Right) )) {

                DUMPER["test"]("extend_checkt: 0\n");
                auto oltype = GetOlReason(o).type;
                if (oltype != OlReason::RS_OK) {
                    auto qri = read_infos_.find(qread.id);
                    auto tri = read_infos_.find(tid); 
                    if (qri->second.filtered.type == RdReason::RS_CONTAINED || tri->second.filtered.type == RdReason::RS_CONTAINED) {
                        if (group.Size(i) > 1) {
                            for (size_t j = 0; j < group.Size(i); ++j) {
                                auto idup = group.Get(i, j);
                                DUMPER["test"]("extend_insert dupcheck: %d, %s\n", idup->Location(qread.id, 0), o.ToM4Line().c_str());
                                if ((end == 0 && idup->Location(qread.id, 0) == Overlap::Loc::Left) ||
                                    (end == 1 && idup->Location(qread.id, 0) == Overlap::Loc::Right)) {
                                        
                                        DUMPER["test"]("extend_insert dup: %s\n", o.ToM4Line().c_str());
                                        extend.insert(idup);
                                        break;
                                }
                                
                            }
                        }
                    }

                }
                    
            } else {
                DUMPER["test"]("extend_checkt: 1\n");
                auto oltype = GetOlReason(o).type;
                if (oltype == OlReason::RS_OK) {
                    
                    DUMPER["test"]("extend_insert 0: %s\n", o.ToM4Line().c_str());
                    extend.insert(&o);
                } else if (oltype == OlReason::RS_FILTERED_READ) {
                                    
                    DUMPER["test"]("extend_insert 1\n");
                    auto qri = read_infos_.find(qread.id);
                    auto tri = read_infos_.find(tid); 
                    if (qri->second.filtered.type == RdReason::RS_CONTAINED || tri->second.filtered.type == RdReason::RS_CONTAINED) {
                        extend.insert(&o);
                        DUMPER["test"]("extend_insert 1: %s\n", o.ToM4Line().c_str());
                    }
                }
            }

        }
    }

    return extend;
}

std::unordered_set<const Overlap*> AsmDataset::GetBackOverlaps(Seq::Id tid, int end) const {
    std::unordered_set<const Overlap*> extend;

    auto group = grouper_.Get(tid);
    if (!group.Empty()) {
        for (size_t i = 0; i < group.Size(); ++i) {
            auto &o = *group.Get(i, 0);
            auto &qread = o.GetOtherRead(tid);

           //if ((end == 0 && o.Location(qread.id, 0) == Overlap::Loc::Left) ||
           //     (end == 1 && o.Location(qread.id, 0) == Overlap::Loc::Right)) continue;

           if ((end == 0 && o.Location(qread.id, 0) != Overlap::Loc::Right) ||
               (end == 1 && o.Location(qread.id, 0) != Overlap::Loc::Left)) continue;

            auto oltype = GetOlReason(o).type;
            if (oltype == OlReason::RS_OK) {
                extend.insert(&o);
            } else if (oltype == OlReason::RS_FILTERED_READ) {
                                
                auto qri = read_infos_.find(qread.id);
                auto tri = read_infos_.find(tid); 
                if (qri->second.filtered.type == RdReason::RS_CONTAINED || tri->second.filtered.type == RdReason::RS_CONTAINED) {
                    extend.insert(&o);
                }
            }
        }
    }

    return extend;
}

void AsmDataset::ReplaceOverlapInGroup(const Overlap* new_ol, const Overlap* old_ol) {
    assert(new_ol->a_.id == old_ol->a_.id && new_ol->b_.id == old_ol->b_.id);

    LOG(WARNING)("todo "); 

    //groups_[new_ol->a_.id][new_ol->b_.id] = new_ol;
    //groups_[new_ol->b_.id][new_ol->a_.id] = new_ol;
}

std::unordered_set<Seq::Id> AsmDataset::ReservedReads() {
    std::unordered_set<Seq::Id> reserved;

    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);
        if (GetOlReason(o).type == OlReason::RS_OK) {
            reserved.insert(o.a_.id);
            reserved.insert(o.b_.id);
        }
    }

    return reserved;
}

void AsmDataset::Dump() const {
    LOG(INFO)("Dump infos");
    
    DumpReadInfos(OutputPath("readinfos"), read_infos_); 
    DumpOverlaps(OutputPath("filter.paf"));
    DumpFilteredOverlaps(OutputPath("filtered_overlaps.txt"));

}


void AsmDataset::DumpFilteredOverlaps(const std::string &fname) const {
    GzFileWriter writer(fname);

    std::mutex mutex;
    std::atomic<size_t> index { 0 };

    auto combine_func = [&mutex, &writer](std::ostringstream &oss) {
        std::lock_guard<std::mutex> lock(mutex);
        writer.Flush(oss);
    };

    auto work_func = [&](int tid) {
        std::ostringstream oss;


        size_t curr = index.fetch_add(1);
        while (curr < ol_store_.Size()) {
            
            const auto &o = ol_store_.Get(curr);
            OlReason rs = GetOlReason(o);
            switch(rs.type) {    
            case OlReason::RS_FILTERED_READ:
                oss << rd_store_.QueryNameById(o.a_.id) << " " << rd_store_.QueryNameById(o.b_.id) << " "
                    << rs.ToString() << " " << rd_store_.QueryNameById(rs.sub[0]) << " " << rs.sub[1] << "\n";
                break;

            case OlReason::RS_SIMPLE:
            case OlReason::RS_DUPLICATE:
            case OlReason::RS_LOCAL:
            case OlReason::RS_CONSISTENCY:
            case OlReason::RS_CONSISTENCY1:
            case OlReason::RS_CONTIG:
            case OlReason::RS_UNKNOWN:
                oss << rd_store_.QueryNameById(o.a_.id) << " " << rd_store_.QueryNameById(o.b_.id) << " "
                    << rs.ToString() << " " << rs.sub[0] << " " << rs.sub[1] << "\n";
                break;
            case OlReason::RS_OK:
            default:
                break;
            }

            if (oss.tellp() > 10000000) {
                combine_func(oss);
            }
            curr = index.fetch_add(1);
        }
        combine_func(oss);

    };

    if (writer.Valid()) {
        MultiThreadRun(opts_.thread_size, work_func);
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}

void AsmDataset::DumpReadInfos(const std::string &fname, const std::unordered_map<int, ReadStatInfo> &readInfos) const {
    GzFileWriter writer(fname);

    if (writer.Valid()) {
        for (const auto &i : readInfos) {
            const auto &ri = i.second;// " "  << ri.overhang << << ri.identity << " "
            writer << rd_store_.QueryNameById(i.first) << " " << ri.len <<  " " << ri.count
                   << " " << ri.overhang_l_threshold << " " << ri.overhang_r_threshold << " " 
                   <<  ri.minmax_coverage[0] << " " << ri.minmax_coverage[1] << " " <<  ri.coverage[0] << "," << ri.coverage[1] << "," << ri.coverage[2] << " "
                   << ri.filtered.ToString() << " " << (ri.filtered.type == RdReason::RS_CONTAINED ? rd_store_.QueryNameById(ri.filtered.sub[0]) : "0") << " "
                   << ri.cliff[0] << " " << ri.cliff[1];
            for (const auto &i : ri.identity_threshold) {
                writer << " " << i;
            } 
            writer << "\n";
        }
    } else {
        LOG(ERROR)("Fail to open ReadInfos file %s", fname.c_str());
    }
}

void AsmDataset::SetOlReason(const Overlap &o, OlReason rs) {
    o.attached = ((long long)(rs.type) << 32) + rs.sub[0];
    assert(rs.sub[1] == 0);
}

OlReason AsmDataset::GetOlReason(const Overlap &o) {
    OlReason rs; 
    rs.type = (OlReason::Type)(o.attached >> 32);
    assert(rs.type >= OlReason::RS_OK && rs.type <= OlReason::RS_UNKNOWN);

    rs.sub[0] = o.attached & 0xFFFFFFFF;
    return rs;
}


int AsmDataset::Percentile(const std::vector<int> &data, double percent) {
    assert(0 <= percent && percent <= 100);

    auto minmaxv = std::minmax_element(data.begin(), data.end());

    auto minv = (*minmaxv.first);
    auto maxv = (*minmaxv.second);

    std::vector<int> counts(maxv-minv+1, 0);
    for (auto c : data) {
        counts[c-minv]++;
    }
    
    int accu = 0;
    for (size_t i=0; i<counts.size(); ++i) {
        accu += counts[i];
        if (data.size() * percent / 100 <= accu) {
            return i + minv;
        }
    }
    return maxv;
}

int AsmDataset::FirstTrough(const std::vector<int> &data, size_t last, size_t k) {
    assert(k % 2 == 1 && data.size() >= k);

    auto minmaxv = std::minmax_element(data.begin(), data.end());

    auto minv = (*minmaxv.first);
    auto maxv = (*minmaxv.second);

    std::vector<int> counts(maxv-minv+1, 0);
    for (auto c : data) {
        counts[c-minv]++;
    }
    
    // calc the starting poistion
    size_t s = 0;
    for (size_t i=1; i<last/10; ++i) {
        if (counts[i-1] > counts[i]) {
            s = i - 1;
            break;
        }
    }

    int value = std::accumulate(counts.begin()+s, counts.begin()+s+k, 0);
    std::pair<size_t, int> best(s, value);
    for (size_t i=s+1; i<counts.size()+k-1 && i<last; ++i) {
        value += -counts[i-1] + counts[i+k-1];
        if (value < best.second*1.00) {
            best.first = i;
            best.second = value;
        } else {
            break;
        }
    }

    assert(best.first >= s);
    size_t bestbest = best.first;
    for (auto i=bestbest+1; i<best.first+k; ++i) {
        if (counts[i] < counts[bestbest]) bestbest = i;
    }


    double rate0 = 0.15;
    double rate1 = 0.5;
    int select = bestbest;
    for (; select-1 >= 0; select--) {
        if (counts[select-1] < counts[select] * (1-rate0) || 
            counts[select-1] > counts[select] * (1+rate0) || 
            counts[select-1] < counts[bestbest] * (1-rate1) ||
            counts[select-1] > counts[bestbest] * (1+rate1))

            break;
    }

    return std::max(select, 1);     // avoid selecting 0
}


void AsmDataset::EstimateGenomeSize() {
    std::vector<int> covs;
    long long int size = 0;
    for (auto &ri : read_infos_) {
        covs.push_back(ri.second.coverage[1]);
        size += ri.second.len;
    }
    std::sort(covs.begin(), covs.end());
    int ave_cov = covs[covs.size()/2];
    long long int gsize = size / ave_cov;
    LOG(INFO)("Esitmate genome size(%lld): %lld = %lld / %d", opts_.genome_size, gsize, size, ave_cov);

    opts_.UpdateByGenomeSize(gsize);

}

void AsmDataset::TestOverlapIdentity() {

    size_t count = std::min<size_t>(20, ol_store_.Size());

    std::vector<const Overlap*> ols(count, nullptr);

    std::default_random_engine e;
    std::uniform_int_distribution<int> u(0, ol_store_.Size()-1);
    e.seed(time(0));
    
    std::generate(ols.begin(), ols.end(), [this, &u, &e](){ return &ol_store_.Get(u(e)); });

    double diff = 1;
    for (const auto o : ols) {
        diff = std::min(diff, std::abs(GetOverlapQuality0(*o) - o->Identity()));
    }
    is_ol_accurate = diff < 0.005;

    if (is_ol_accurate) {
        LOG(INFO)("The identity of overlaps is accurate %.02f", diff / count);
    } else {
        LOG(INFO)("The identity of overlaps is inaccurate %.02f", diff / count);
    }
}


void AsmDataset::EstimateCoverage() {
    const int N = 40;
    size_t count = std::min<size_t>(N, rd_store_.Size());
    std::vector<Seq::Id> reads(count);

    std::default_random_engine e;
    std::uniform_int_distribution<int> u(0, rd_store_.Size()-1);
    e.seed(time(0));
    
    std::generate(reads.begin(), reads.end(), [this, &u, &e](){ return u(e); });

    std::vector<double> average_coverages;
    for (const auto rd : reads) {
        auto len = rd_store_.GetSeqLength(rd);
        std::vector<int> covs(len+1);
        auto gp = grouper_.Get(rd);
        for (size_t i = 0; i < gp.Size(); ++i) {
            for (size_t j = 0; j < gp.Size(i); ++j) {
                auto ol = gp.Get(i, j);
                auto read = ol->GetRead(rd);
                covs[read.start] += 1;
                covs[read.end] -= 1;
            }   
        }
        for (size_t i = 1; i <covs.size(); ++i) {
            covs[i] += covs[i-1];
        }
        assert(covs.back() == 0);

        auto ave = std::accumulate(covs.begin(), covs.end(), 0.0) * 1.0 / covs.size();
        average_coverages.push_back(ave);
    }

    std::sort(average_coverages.begin(), average_coverages.end());
    size_t start = average_coverages.size() / 2 - average_coverages.size() / 10 ;
    size_t end = average_coverages.size() / 2  + average_coverages.size() / 10;
    for (auto c : average_coverages) {
        LOG(INFO)("coverage is  %.02f", c);
    }
    average_coverage_ = std::accumulate(average_coverages.begin() + start, average_coverages.begin() + end, 0.0) / (end - start);

    LOG(INFO)("The coverage is  %.02f, %zd %zd", average_coverage_, start, end);

}



double AsmDataset::EstimateCoverage(const std::unordered_map<Seq::Id, ReadStatInfo>& read_infos) {

    std::vector<int> min_covs;
    min_covs.reserve(read_infos.size());
    std::vector<int> max_covs;
    max_covs.reserve(read_infos.size());

    for (const auto &rinfo : read_infos) {
        min_covs.push_back(rinfo.second.minmax_coverage[0]);
        max_covs.push_back(rinfo.second.minmax_coverage[1]);
    }

    std::sort(min_covs.begin(), min_covs.end());
    std::sort(max_covs.begin(), max_covs.end());
    size_t start = 0;
    size_t end = min_covs.size();
    if (min_covs.size() >= 10) {
        start =  min_covs.size() / 2 - min_covs.size() / 6;
        end =    min_covs.size() / 2 + min_covs.size() / 6;
    }

    double accu = std::accumulate(min_covs.begin()+start, min_covs.begin()+end, 0)+
                  std::accumulate(max_covs.begin()+start, max_covs.begin()+end, 0);
    return accu * 1.0 / 2 / (end - start);

}

void AsmDataset::ClusterBundle(const std::unordered_set<Seq::Id> &bundle) {
    // 统计所有kmer
    const int klen = 15;
    LOG(INFO)("Start counting kmer count");
    std::unordered_map<KmerId, size_t> kmer_counts;
    KmerCounter kc(klen);
    for (auto rid : bundle) {
        auto rseq = rd_store_.GetSeq(rid);
        auto kmers = kc.CountAll(rseq);
        for (const auto &k : kmers) {
            auto mink = std::min(k[0], k[1]);
            kmer_counts[mink] ++;
        }
    }
    LOG(INFO)("End counting kmer count");


    // 计算有效kmer
    std::unordered_set<KmerId> valid;
    for (const auto& it : kmer_counts) {
        if (it.second > average_coverage_/2 && it.second < average_coverage_*1.5 ) {
            valid.insert(it.first);
        }
    }
    LOG(INFO)("Valid kmer size = %zd", valid.size());

    std::unordered_map<Seq::Id, std::vector<std::tuple<size_t, KmerId>>> rd_kmers;
    for (auto rid : bundle) {
        auto rseq = rd_store_.GetSeq(rid);
        auto kmers = kc.CountAll(rseq);
        rd_kmers[rid] = std::vector<std::tuple<size_t, KmerId>>();
        for (size_t i = 0; i < kmers.size(); ++i) {
            const auto &k = kmers[i];
            auto mink = std::min(k[0], k[1]);
            if (valid.find(mink) != valid.end()) {
                rd_kmers[rid].push_back(std::make_tuple(i, mink));
            }
        }
    }
    LOG(INFO)("End count read kmer, %d", klen);

    auto jaccard = [](const std::unordered_set<KmerId>& a, const std::unordered_set<KmerId>& b) {
        size_t comm = 0;
        for (auto & i : a) {
            if (b.find(i) != b.end()) comm ++;
        }

        if (a.size() == 0 && b.size() == 0) return 1.1;
        return comm * 1.0 / (a.size() + b.size() - comm);
    };

    auto get_kmers = [klen](const std::vector<std::tuple<size_t, KmerId>>& kmers, int start, int end) {
        std::unordered_set<KmerId> kset;

        for (auto& i : kmers) {
            if (std::get<0>(i) > start + klen) {
                if (std::get<0>(i) < end - klen) {
                    kset.insert(std::get<1>(i));
                } else {
                    break;
                }
            } 
        }
        return kset;

    };

    auto overlap_jaccard = [this, get_kmers, jaccard](const std::vector<std::tuple<size_t, KmerId>>& a, const Overlap::Read &ard, 
                                     const std::vector<std::tuple<size_t, KmerId>> &b, const Overlap::Read &brd) {
        std::unordered_set<KmerId> aset = get_kmers(a, ard.start, ard.end);
        std::unordered_set<KmerId> bset = get_kmers(b, brd.start, brd.end);
        
        // DUMPER["data"]("jac %s(%zd) - %s(%zd)\n", rd_store_.QueryNameById(ard.id).c_str(), aset.size(),
        //     rd_store_.QueryNameById(brd.id).c_str(), bset.size());
        return jaccard(aset, bset);
    };
    

    
    std::mutex mutex;
    auto combine_func = [this, &mutex](const std::unordered_set<const Overlap*> &ignored) {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto ol : ignored) {
            SetOlReason(*ol, OlReason::Consistency(0));
        }
    };
    
    std::vector<int> bundle_list(bundle.begin(), bundle.end());
    std::atomic<size_t> index {0};
    auto work_func = [this, &bundle, &index, &bundle_list, &rd_kmers, combine_func, overlap_jaccard](size_t tid) {
        std::unordered_set<const Overlap*> ignored;
        for (size_t i = index.fetch_add(1); i < bundle_list.size(); i = index.fetch_add(1)) {
            auto rid = bundle_list[i];       
            auto gp = grouper_.Get(rid);
            for (size_t i = 0; i < gp.Size(); ++i) {
                for (size_t j = 0; j < gp.Size(i); ++j) {
                    auto ol = gp.Get(i, j);
                    auto& read = ol->GetRead(rid);
                    auto& other = ol->GetOtherRead(rid);
                    if (bundle.find(other.id) == bundle.end()) continue;
                    auto ord_kmer = rd_kmers.find(other.id);
                    auto rd_kmer = rd_kmers.find(rid);
                    if (ord_kmer != rd_kmers.end() && rd_kmer != rd_kmers.end()) {
                        auto jac = overlap_jaccard(ord_kmer->second, other, rd_kmer->second, read); 
                            
                        DUMPER["data"]("jac  %s - %s %.02f\n", rd_store_.QueryNameById(rid).c_str(),
                            rd_store_.QueryNameById(other.id).c_str(), jac);
                        
                        if (jac < 0.5) {
                            ignored.insert(ol);
                        }
                    }
                }
            }
        }
        combine_func(ignored);
    };
    

    MultiThreadRun(opts_.thread_size, work_func);
}

} // namespace fsa {
