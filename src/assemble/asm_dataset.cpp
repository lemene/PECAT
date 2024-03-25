#include "asm_dataset.hpp"

#include <algorithm>
#include <memory>
#include "../utility.hpp"

#include "edlib.h"

namespace fsa {
    
void AsmDataset::Load() {
    if (!opts_.read_file.empty()) {
        rd_store_.Load(opts_.read_file, "", true);
        rd_store_.SaveIdToName(opts_.OutputPath("id2name.txt.gz"));
    }
    
    LoadOverlaps(opts_.ifname);
    DumpReadInfos(OutputPath("readinfos"), read_infos_); 
}

void AsmDataset::LoadPurged() {
    if (!opts_.read_file.empty()) {
        rd_store_.Load(opts_.read_file, "", true);
        rd_store_.SaveIdToName(opts_.OutputPath("id2name.txt.gz"));
    }
    
    
    if (rd_store_.Size() > 0) {    // read names have been loaded
        ol_store_.LoadFast(OutputPath("filter.m4a"), "", (size_t)opts_.thread_size);
    } else {
        ol_store_.Load(OutputPath("filter.m4a"), "", (size_t)opts_.thread_size);
    }

    LOG(INFO)("Overlap size: %zd ", ol_store_.Size());
}

void AsmDataset::Purge() {
    
    GroupAndFilterDuplicate();

    FilterLowQuality();

    //CheckOverlapEnd();

    ExtendOverlapToEnd();

    FilterCoverage();

    EstimateGenomeSize();

    //FilterConsistency();

    FilterContained();
    
    Dump();
}

void AsmDataset::LoadOverlaps(const std::string &fname) {
    
    bool scan_overlaps_ = false;
    if (scan_overlaps_) {
        LOG(INFO)("Scan overlaps to select params");
        ScanOverlapsToSelectParams(opts_.ifname);
        LOG(INFO)("Load overlap file");
        LoadOverlapsWithoutLowQuality(opts_.ifname);
    } else {
        
        LOG(INFO)("Load overlap file");
        LoadOverlaps2(opts_.ifname);
    }

    
    if (opts_.dump >= 2) {
        DumpOverlaps(OutputPath("load.m4a"));
    }
}

    
void AsmDataset::LoadOverlaps2(const std::string &fname) {
    std::atomic<size_t> count {0};
    auto check = [&](Overlap& o) {
        count++;
        return opts_.filter0.Valid(o);
    };

    if (rd_store_.Size() > 0) {    // read names have been loaded
        ol_store_.LoadFast(fname, "", (size_t)opts_.thread_size, check);
    } else {
        ol_store_.Load(fname, "", (size_t)opts_.thread_size, check);
    }

    LOG(INFO)("Overlap size: %zd/%zd", ol_store_.Size(), count.load());

    //
    
    for (size_t i=0; i<rd_store_.Size(); ++i) {
        read_infos_[i] = ReadStatInfo();
        read_infos_[i].len = rd_store_.GetSeqLength(i);
        read_infos_[i].id = i; // TODO
    }
}


void AsmDataset::LoadOverlapsWithoutLowQuality(const std::string &fname) {
    std::atomic<size_t> total { 0 };
    auto filter_simple = [&](Overlap& o) {
        if (opts_.filter0.Valid(o)) {
            total ++;
            // TODO not filtering
            // const auto& rinfo_a = read_infos_[o.a_.id];
            // const auto& rinfo_b = read_infos_[o.b_.id];

            // if (!rinfo_a.CheckIdentity(o) && !rinfo_b.CheckIdentity(o)) {
            //     return false;
            // }
            
            return true;
        } else {
            return false;
        }
    };

    if (rd_store_.Size() > 0) {    // read names have been loaded
        ol_store_.LoadFast(fname, "", (size_t)opts_.thread_size, filter_simple);
    } else {
        ol_store_.Load(fname, "", (size_t)opts_.thread_size, filter_simple);
    }

    LOG(INFO)("Overlap size: %zd / %zd", ol_store_.Size(), total.load());
}

void AsmDataset::ScanOverlapsToSelectParams(const std::string &fname) {

    std::mutex mutex;

    size_t block_size = 50000;
    struct WorkArea{
        StatReadInfo sri;
        void Clear() { sri.Clear(); }
    };
    
    std::list<std::shared_ptr<WorkArea>> works;  // for each thread. Vector may cause memory reallocating, so list is used.

    auto alloc_work = [&]() -> std::shared_ptr<WorkArea> {
        std::lock_guard<std::mutex> lock(mutex);
        works.push_back(std::shared_ptr<WorkArea>(new WorkArea()));

        return works.back();
    };

    auto combine = [&](std::shared_ptr<WorkArea> work) {
        std::lock_guard<std::mutex> lock(mutex);
        sread_info_.Merge(work->sri);
    };

    auto scan_overlap = [&](Overlap& o) {
        auto thread_local work = alloc_work();
        if (opts_.filter0.Valid(o)) {
            work->sri.Add(o);
        }
        if (work->sri.Size() >= block_size) {
            combine(work);
        }
        return false;   // not load to memory
    };

    OverlapStore ol(rd_store_.GetStringPool());
    if (rd_store_.GetIdRange()[1] > 0) {    // read names have been loaded
        ol.LoadFast(fname, "", (size_t)opts_.thread_size, scan_overlap);
    } else {
        ol.Load(fname, "", (size_t)opts_.thread_size, scan_overlap);
    }

    for (auto & w : works) {
        combine(w);
    }
    
    auto work_func = [&](const std::vector<int>& input) {
        for (auto i : input) {
            read_infos_[i].Stat(i, opts_);
        }
    };


    MultiThreadRun((size_t)opts_.thread_size, read_infos_, SplitMapKeys<std::unordered_map<int, fsa::ReadStatInfo>>, work_func);   
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
    auto work_func = [&](const std::vector<int>& input) {
        std::unordered_set<const Overlap*> filtered;
        for (auto i : input) {
            FilterLowQuality(i, groups_[i], filtered);
        }
        return filtered;
    };

    auto combine_func = [](const std::vector<std::unordered_set<const Overlap*>>& data) {
        std::unordered_set<const Overlap*> filtered;
        for (auto d : data) {
            filtered.insert(d.begin(), d.end());

        }
        return filtered;
    };

    auto filtered = MultiThreadRun(opts_.thread_size, groups_, 
        SplitMapKeys<decltype(groups_)>, 
        work_func, 
        combine_func);   


    LOG(INFO)("Filtered low-quality overlaps: %zd", filtered.size());
    for (auto o : filtered) {
        SetOlReason(*o, OlReason::Simple());
    }
}

void AsmDataset::FilterLowQuality(int id, const std::unordered_map<int, const Overlap*> &group, std::unordered_set<const Overlap*> &ignored) {
    const int winsize = 4000;

    auto& rinfo = read_infos_[id];
    std::vector<std::vector<double>> winidents((rinfo.len + winsize / 2) / winsize);
    std::vector<double> lohs;
    std::vector<double> rohs;

    // collect target read information: 
    for (const auto& g : group) {
        auto &ol = *g.second;

        if (IsReserved(ol)) {
            auto &tr = ol.GetRead(id);

            size_t s = (tr.start + winsize / 2) / winsize;
            size_t e = (tr.end + winsize / 2) / winsize;
            std::for_each(winidents.begin()+s, winidents.begin()+e, [&ol](std::vector<double>& v) {
                v.push_back(ol.identity_);
            });
            
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

    std::vector<double> identity_threshold(winidents.size());
    std::transform(winidents.begin(), winidents.end(), identity_threshold.begin(), [this, id](std::vector<double>& ident) {
        if (ident.size() > 0) {
            std::sort(ident.begin(), ident.end(), [](double a, double b) { return a > b; });

            double median, mad;
            ComputeMedianAbsoluteDeviation(std::vector<double>(ident.begin(), ident.begin()+std::min<size_t>(30, ident.size())),  median, mad);
            if (rd_store_.QueryNameById(id) == "740780") {
                LOG(INFO)("filter_low_quality: size=%zd %0.02f, %0.02f, %0.02f",ident.size(), median, mad, std::max(opts_.filter0.min_identity, median-6*1.4826*mad));
            }
            return std::max(opts_.filter0.min_identity, median-6*1.4826*mad);

        } else {
            return opts_.filter0.min_identity;
        }
    });

    double overhang_l_threshold = 0;
    double overhang_r_threshold = 0;
    
    if (lohs.size() > 0) {
        double median, mad;
        ComputeMedianAbsoluteDeviation(lohs, median, mad);
        overhang_l_threshold = std::min<double>(median + 6*1.4826*mad, opts_.filter0.max_overhang);
    } else {
        overhang_l_threshold = opts_.filter0.max_overhang;
    }

    if (rohs.size() > 0) {
        double median, mad;
        ComputeMedianAbsoluteDeviation(rohs, median, mad);
        overhang_l_threshold = std::min<double>(median + 6*1.4826*mad, opts_.filter0.max_overhang);
    } else {
        overhang_l_threshold = opts_.filter0.max_overhang;
    }

    rinfo.overhang_l_threshold = overhang_l_threshold;
    rinfo.overhang_r_threshold = overhang_r_threshold;

    auto area_threshold = [winsize](const std::vector<double>& idents, int start, int end) {
 
        size_t s = (start + winsize / 2) / winsize;
        size_t e = (end + winsize / 2) / winsize;
        assert (s >= 0 && s <= e && e <= idents.size());

        size_t count = 0;
        double sum = 0;
        for (auto i = s; i< e; ++i) {
            if (idents[i] > 0) {
                count ++;
                sum += idents[i];
            }
        }

        return sum / count;
    };

    for (const auto& g : group) {
        auto &ol = *g.second;
        if (IsReserved(ol)) {
            auto &tr = ol.GetRead(id);
            if (ol.identity_ < area_threshold(identity_threshold, tr.start, tr.end)) {
                ignored.insert(&ol);
            }
        }
    }


}

void AsmDataset::DumpOverlaps(const std::string &fname) const {
    
    auto filter = [&](const Overlap& o)->bool {
        return GetOlReason(o).type == OlReason::RS_OK;
    };

    ol_store_.Save(fname, "", opts_.thread_size, filter);
}


void AsmDataset::ModifyEnd(const Overlap &oldone, int maxoh) {
    if (oldone.Location(0) != Overlap::Loc::Abnormal) return;
    if (oldone.Location(maxoh) == Overlap::Loc::Abnormal) return;

    Overlap& o = const_cast<Overlap&>(oldone);
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
        for (size_t i = index.fetch_add(1); i < ol_store_.Size(); i = index.fetch_add(1)) {
            const Overlap& o = ol_store_.Get(i);

            if (IsReserved(o)) {
                //auto loc = o.Location(opts_.filter0.max_overhang);
                auto loc = o.Location(0);
                if (loc == Overlap::Loc::Equal) {
                    set_contained(std::max(o.a_.id, o.b_.id), std::min(o.a_.id, o.b_.id));
                } else if (loc == Overlap::Loc::Contained) {
                    set_contained(o.a_.id, o.b_.id);
                } else if (loc == Overlap::Loc::Containing) {
                    set_contained(o.b_.id, o.a_.id);
                }
            }
        }
    };

    MultiThreadRun((size_t)opts_.thread_size, work_func);

    LOG(INFO)("Start filtering contained reads and relative overlaps");
    UpdateFilteredRead();
}


void AsmDataset::FilterCoverage() {
    LOG(INFO)("Check Coverage");

    auto work_func = [&](const std::vector<int>& input) {

        for (auto i : input) {
            auto minmax = AnalyzeCoverage(i, groups_[i]);
            auto iter = read_infos_.find(i);
            if (iter != read_infos_.end()) {
                iter->second.minmax_coverage = {minmax[1], minmax[2]};
                iter->second.covtype = minmax[0];
            }
        }
    };

    MultiThreadRun(opts_.thread_size, groups_, SplitMapKeys<decltype(groups_)>, work_func);  

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


void AsmDataset::FilterConsistency() {
    LOG(INFO)("Check consistency");

    auto work_func = [this](const std::vector<int> &input) {
        std::array<std::unordered_set<const Overlap*>, 2> output;  // best, ignored
        for (auto id : input) {
            if (IsReserved(id)) {
                auto &g = groups_[id];
                CalcConsistency(id, g, output[0], output[1]); 
            }
        }
        
        return output;
    }; 
    
    auto comb_func = [this](const std::vector<std::array<std::unordered_set<const Overlap*>, 2>> &sub_outputs) {
        std::array<std::unordered_set<const Overlap*>, 2> output;

        for (const auto& i : sub_outputs) {
            output[0].insert(i[0].begin(), i[0].end());
            output[1].insert(i[1].begin(), i[1].end());
        }
        return output;
    };

    auto output = MultiThreadRun(opts_.thread_size, groups_,
        SplitMapKeys<decltype(groups_)>, 
        work_func, 
        comb_func);

    for (const auto &o : output[1]) {
        SetOlReason(*o, OlReason::Consistency(0));
    }
    
    
    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);
        if (IsReserved(o) && output[0].find(&o) == output[0].end()) {
            SetOlReason(o, OlReason::Consistency1(0));
        }
    }
    
}

void Debug_PrintCluster(const std::vector<std::set<int>>& clusters) {

    for (const auto &i : clusters) {
        printf("type : ");
        for (auto ii : i) {
            printf("%d ", ii);
        }
        printf("\n");
    }
}

void AsmDataset::CalcConsistency(int id, const std::unordered_map<int, const Overlap*> &group, std::unordered_set<const Overlap*> &best, std::unordered_set<const Overlap*> &ignored) {
    std::vector<const Overlap*> left_ols;
    std::vector<const Overlap*> right_ols;

    long min_extension = 100;

    
    for (auto o : group) {
        if (IsReserved(*o.second)) {
            auto &a = o.second->GetRead(id);
            auto &b = o.second->GetOtherRead(id);

            if (a.start == 0 && ((a.strand == b.strand && b.start > min_extension) || (a.strand != b.strand && b.len-b.end > min_extension))) {
                left_ols.push_back(o.second);

            } else if (a.end == a.len && ((a.strand != b.strand && b.start > min_extension) || (a.strand == b.strand && b.len-b.end > min_extension))) {
                right_ols.push_back(o.second);
            } else {
                best.insert(o.second);
            }
        }
    }

    auto left_graph = CalcConsistencyGraph(id, left_ols);
    auto left_clusters = left_graph.Cluster();
    auto right_graph = CalcConsistencyGraph(id, right_ols);
    auto right_clusters = right_graph.Cluster();


    auto max_aligned_length = [id](const std::set<int>& s, const std::vector<const Overlap*>& ols) {
        int a = 0;
        for (auto i : s) {
            auto& r = ols[i]->GetRead(id);
            if (r.end -r.start > a) {
                a = r.end - r.start;
            }
        }
        return a;
    };

    std::sort(left_clusters.begin(), left_clusters.end(), [max_aligned_length, &left_ols](const std::set<int>& a, const std::set<int>& b){
        
        //return max_aligned_length(a, left_ols) > max_aligned_length(b, left_ols);
        return a.size() > b.size();
    });
    
    std::sort(right_clusters.begin(), right_clusters.end(), [max_aligned_length, &right_ols](const std::set<int>& a, const std::set<int>& b){
        //return max_aligned_length(a, right_ols) > max_aligned_length(b, right_ols);
        return a.size() > b.size();
    });

    if ( opts_.debug_name == rd_store_.QueryNameById(id)) {
        printf("Left\n");
        Debug_PrintGraph(opts_.debug_name, left_graph, left_clusters, left_ols);
        printf("Right\n");
        Debug_PrintGraph(opts_.debug_name, right_graph, right_clusters, right_ols);
    }

    SelectBestExtends(id, left_graph, left_clusters, left_ols, best, ignored);
    SelectBestExtends(id, right_graph, right_clusters, right_ols, best, ignored);
}

void AsmDataset::SelectBestExtends(Seq::Id id, const MatrixGraph& graph, const std::vector<std::set<int>>& clusters, const std::vector<const Overlap*>& ols,
                                      std::unordered_set<const Overlap*> &best, std::unordered_set<const Overlap*> &ignored) {
    int cov = read_infos_[id].coverage[1];

    int accu = 0;
    bool pass = false;
    for (const auto &clu : clusters) {
        
    

        // 每一个cluster选择最小的集合
        auto ibest= *clu.begin();
        auto o = ols[ibest];
        auto &b = o->a_.id != id ? o->a_  : o->b_;
        int c = read_infos_[b.id].coverage[1];
        if (opts_.debug_name == rd_store_.QueryNameById(id)) {
            printf("l cov: %d, %d, %d\n", cov, c, accu );
        }
        if (!pass && (clu.size() >= 1 || clusters.size() == 1) && (accu == 0 || accu + c< cov*1.5) && clu.size() >= clusters[0].size() / 3) {

            accu += c;
            if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                printf("+ cov: %d, %d, %d\n", cov, c, accu );
            }

            for (auto i : clu) {
                if (graph.Degree(i) > graph.Degree(ibest)) {
                    ibest = i;
                }
                //if (ols[i]->AlignedLength() > o->AlignedLength()) {
                //    o = ols[i];
                //}
            }
            best.insert(ols[ibest]);

            for (auto i : clu) {
                if (graph.Degree(i) >= clu.size() / 3){//graph.Degree(ibest) ) {   
                    best.insert(ols[i]);
                    if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                            printf("+ best: %s %s\n", rd_store_.QueryNameById(ols[i]->a_.id).c_str(), rd_store_.QueryNameById(ols[i]->b_.id).c_str());
                    }
                } else if (graph.Degree(i)  < graph.Degree(ibest) / 3) {
                    
                    ignored.insert(ols[i]);
                    
                    if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                            printf("+ ignore0: %s %s\n", rd_store_.QueryNameById(ols[i]->a_.id).c_str(), rd_store_.QueryNameById(ols[i]->b_.id).c_str());
                    }
                }

            }
            

        } else {
            
            if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                printf("+ ignore: %zd\n", clu.size() );
                
                for (auto i : clu) {
                    printf("+ ignore: %s %s\n", rd_store_.QueryNameById(ols[i]->a_.id).c_str(), rd_store_.QueryNameById(ols[i]->b_.id).c_str());
                }
            }
            for (auto i : clu) {
                ignored.insert(ols[i]);
            }
            pass = true;
        }
    }

}

void AsmDataset::Debug_PrintGraph(const std::string &name, const MatrixGraph& graph, const std::vector<std::set<int>>& clusters, const std::vector<const Overlap*>& ols) {
    
    printf("debug name %s %zd\n", name.c_str(), ols.size());
    graph.Print();
    Debug_PrintCluster(clusters);

    for (auto o : ols) {
        printf("%s\n", OverlapStore::ToPafLine(*o, StringPool::UnsafeNameId(rd_store_.GetStringPool())).c_str());
    }
}

int a = 0;
MatrixGraph AsmDataset::CalcConsistencyGraph(int id, const std::vector<const Overlap*> ols) {
    double rate = opts_.max_offset_rate;

    MatrixGraph graph(ols.size());
    for (size_t i=0; i<ols.size(); ++i) {
        const Overlap* o0_ab = ols[i];
        auto &a0 = o0_ab->a_.id == id ? o0_ab->a_  : o0_ab->b_;
        auto &b0 = o0_ab->a_.id != id ? o0_ab->a_  : o0_ab->b_;

        for (size_t j=i+1; j<ols.size(); ++j) {
            const Overlap* o1_ac = ols[j];
            auto &a1 = o1_ac->a_.id == id ? o1_ac->a_  : o1_ac->b_;
            auto &c1 = o1_ac->a_.id != id ? o1_ac->a_  : o1_ac->b_;

            assert((a0.start == 0 && a1.start == 0) || (a0.end == a0.len && a1.end == a1.len));
            auto it1 = groups_.find(b0.id);

            if (it1 != groups_.end()) {
                auto it2 = it1->second.find(c1.id);
                if (it2 != it1->second.end() && IsReserved(*it2->second)) {
                    const Overlap* o2_bc = it2->second;
                    if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                        a = 1;
                    }
                    if (Overlap::IsConsistent(*o0_ab, *o1_ac, *o2_bc, 2000)) {
                        graph.AddEdge(i, j, 1);
                    }
                    if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                        a = 0;
                    }
                    
                    if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                        printf("ccc %d %d %s, %s\n", Overlap::IsConsistent(*o0_ab, *o1_ac, *o2_bc, (int)b0.len*rate) ,(int)(b0.len*rate), 
                            rd_store_.QueryNameById(o0_ab->GetOtherRead(id).id).c_str(),
                            rd_store_.QueryNameById(o1_ac->GetOtherRead(id).id).c_str());
                    }
                } else {
                    
                    if (opts_.debug_name == rd_store_.QueryNameById(id)) {
                        printf("ccc x %s, %s\n",  
                            rd_store_.QueryNameById(o0_ab->GetOtherRead(id).id).c_str(),
                            rd_store_.QueryNameById(o1_ac->GetOtherRead(id).id).c_str());
                    }
                }
            }
        }
    }
    return graph;
}


// void AsmDataset::FilterBestn() {
//     LOG(INFO)("Select Best n overlaps");

//     auto work_func = [&](const std::vector<int> &input) {
        
//         std::unordered_set<const Overlap*> output;
//         for (auto &i : input) {
//             auto iter = read_infos_.find(i);
//             assert(iter != read_infos_.end());
//             if (iter->second.filtered.IsOk()) {
//                 auto g = groups_.find(i);
//                 assert(g != groups_.end());
//                 auto k = FindBestN(*g);

//                 output.insert(k.begin(), k.end());
//             }
            
//         }
//         return output;
//     };


//     auto keep = MultiThreadRun(opts_.thread_size, groups_, 
//         SplitMapKeys<decltype(groups_)>, 
//         work_func, 
//         MoveCombineMapOrSet<std::unordered_set<const Overlap*>>);

//     for (size_t i=0; i < ol_store_.Size(); ++i) {
//         const auto &o = ol_store_.Get(i);
//         if (IsReserved(o) && keep.find(&o) == keep.end()) {
//             SetOlReason(o, OlReason::BestN());
//         }

//     }

// }

double AsmDataset::GetOverlapQuality(const Overlap &ol) {
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

void AsmDataset::GroupAndFilterDuplicate() {
    LOG(INFO)("Group overlaps and remove duplicated");

    std::vector<const Overlap*> removed;
    std::mutex mutex;

    auto add_overlap = [this](int low, int a, int b, const Overlap& o, 
                              std::vector<std::unordered_map<Seq::Id, const Overlap*>>& group,
                              std::vector<std::unordered_map<Seq::Id, std::vector<const Overlap*>>> &dups,
                              std::vector<const Overlap*> &rmd) {
        
        auto it = group[a-low].find(b);
        if (it == group[a-low].end()) {
            group[a-low][b] = &o;
        } else {
            if (BetterAlignedLength(o, *(it->second))) {
                //SetOlReason(*(it->second), OlReason::Duplicate());
                rmd.push_back(it->second);
                it->second = &o;
            } else {
                rmd.push_back(&o);
                //SetOlReason(o, OlReason::Duplicate());
            }
        }
        dups[a-low][b].push_back(&o);
    };

    auto split_func = [this]() {
        auto r = ol_store_.GetReadIdRange();
        return SplitRange(opts_.thread_size, r[0], r[1]);
    };
    auto comb_func = [this, &removed, &mutex](int low, std::vector<std::unordered_map<Seq::Id, const Overlap*>>&& group,
                                              std::vector<std::unordered_map<Seq::Id, std::vector<const Overlap*>>> &&dups,
                                              std::vector<const Overlap*> &rmd) {
        std::lock_guard<std::mutex> lock(mutex);
        assert(group.size() == dups.size());
        for (size_t i=0; i<group.size(); ++i) {
            if (group[i].size() > 0) {
                groups_[low+(int)i] = std::move(group[i]);
            }

            for (auto&& d : dups[i]) {
                if (d.second.size() > 1) {
                    std::sort(d.second.begin(), d.second.end(), [](const Overlap* a, const Overlap* b) {
                        return a->AlignedLength() > b->AlignedLength();
                    });
                    dup_groups_[low+(int)i][d.first] = std::move(d.second);
                }
            }
        }
        removed.insert(removed.end(), rmd.begin(), rmd.end());
    };

    auto work_func = [this, add_overlap, comb_func](std::array<Seq::Id, 2> r) {
        std::vector<std::unordered_map<Seq::Id, const Overlap*>> group(r[1] - r[0]);
        std::vector<std::unordered_map<Seq::Id, std::vector<const Overlap*>>> dups(r[1] - r[0]);
        std::vector<const Overlap*> rmd;    // removed

        for (size_t i=0; i < ol_store_.Size(); ++i) {
            const auto &o = ol_store_.Get(i);
            if (IsReserved(o)) {
                if (o.a_.id >= r[0] && o.a_.id < r[1]) {
                    add_overlap(r[0], o.a_.id, o.b_.id, o, group, dups, rmd);
                }
                if (o.b_.id >= r[0] && o.b_.id < r[1]) {
                    add_overlap(r[0], o.b_.id, o.a_.id, o, group, dups, rmd);
                }
            }
        }

        comb_func(r[0], std::move(group), std::move(dups), rmd);
    };

    MultiThreadRun((int)opts_.thread_size, split_func, work_func);

    for (auto r : removed) {
        SetOlReason(*r, OlReason::Duplicate());
    }
    LOG(INFO)("Overlap size: %zd/%zd", ReservedSize(), ol_store_.Size());

    Check_Group();
}

void AsmDataset::CheckOverlapEnd() {
    LOG(INFO)("Check overlap overhang");

    using FilterData = std::pair<std::unordered_set<const Overlap*>,  std::unordered_set<int>>;


    auto work_func = [&](const std::vector<int>& input) {
        FilterData filtered;

        for (auto i : input) {
            CheckOverlapEnd(i, groups_[i], filtered.first, filtered.second);
        }
        return filtered;
    };

    auto combine_func = [](const std::vector<FilterData>& data) {
        FilterData filtered;
        for (auto d : data) {
            filtered.first.insert(d.first.begin(), d.first.end());
            filtered.second.insert(d.second.begin(), d.second.end());

        }
        return filtered;
    };

    auto filtered = MultiThreadRun(opts_.thread_size, groups_, 
        SplitMapKeys<decltype(groups_)>, 
        work_func, 
        combine_func);   

    printf("CheckOverlapEnd %zd, %zd\n", filtered.first.size(), filtered.second.size());

    for (auto o : filtered.first) {
        SetOlReason(*o, OlReason::Consistency(1));
    }
}

void AsmDataset::CheckOverlapEnd(int id, const std::unordered_map<int, const Overlap*> &group, std::unordered_set<const Overlap*> &ignored, std::unordered_set<int>& ignReads) {

    struct Item {
        const Overlap* o;
        int endlen;
    };

    std::vector<Item> left, right;
    for (auto g : group) {
        const Overlap& o = *g.second;

        if (IsReserved(o)) {
            auto & tr = o.GetRead(id);
            auto ti =  o.a_.id == id ? 0 : 1;
            
            auto oh = o.Overhang2();
            if ((oh[ti] & 1) != 0 || tr.start == 0) {
                left.push_back({&o, tr.start});
            } else if ((oh[ti] & 2) != 0) {
                right.push_back({&o, tr.len-tr.end});
            }
        }
    }

    auto verifyEnd = [&](int id, std::vector<Item>& endlist, std::unordered_set<const Overlap*> &ignored, std::unordered_set<int>& ignReads, int end) {
        std::sort(endlist.begin(), endlist.end(), [](const Item &a, const Item &b) { return a.endlen < b.endlen; });

        const auto& rinfo = read_infos_[id];
        int threshold = end == 0 ? rinfo.overhang_l_threshold : rinfo.overhang_r_threshold;

        auto start = std::find_if(endlist.begin(), endlist.end(), [threshold](const Item& a) { return a.endlen > threshold; });

        for (auto i = start; i != endlist.end(); ++i) {

            const Overlap& o = *i->o;

            bool verified = false;
            for (auto j=endlist.begin(); j<start && !verified; ++j) {
                const Overlap& oab = o;
                const Overlap& oac = *j->o;

                // Seq::Id ia = id;
                Seq::Id ib = oab.GetOtherRead(id).id;
                Seq::Id ic = oac.GetOtherRead(id).id;

                auto citer = groups_[ib].find(ic);
                if (citer != groups_[ib].end() && IsReserved(*citer->second)) {
                    const Overlap &obc = *citer->second;

                    if (Overlap::IsConsistent(oab, oac, obc, opts_.filter0.max_overhang)) {
                        auto rb = oab.GetOtherRead(id);
                        int verify_len = (oab.SameDirect() ^ (end != 0)) ? rb.start - obc.GetRead(ib).start : rb.len-rb.end - (obc.GetRead(ib).len - obc.GetRead(ib).end) ;
                        if (i->endlen - verify_len <= threshold) {
                            verified = true;
                        }
                    }

                }
            }

            if (!verified) {
                ignored.insert(i->o);
            }
        }
    };


    verifyEnd(id, left, ignored, ignReads, 0);
    verifyEnd(id, right, ignored, ignReads, 1);

}
                    
// Modify the ends and remove overhangs
void AsmDataset::ExtendOverlapToEnd() {
    // LOG(INFO)("Extend Overlaps");
    MultiThreadRun(opts_.thread_size, 
        [this]() {
            return SplitRange(opts_.thread_size, (size_t)0, ol_store_.Size());
        }, 
        [this](const std::array<size_t, 2> &range) {
            for (size_t i=range[0]; i<range[1]; ++i) {
                const Overlap& o = ol_store_.Get(i);
                if (IsReserved(o) || GetOlReason(o).type == OlReason::RS_DUPLICATE) {
                    auto oh = o.Overhang();
                    int th = std::max(oh[0], oh[1]);
                    if (th > 0) {
                        ModifyEnd(o, th);
                    }
                }
            }
        }
    );
}

std::array<int, 3> AsmDataset::AnalyzeCoverage(int id, const std::unordered_map<int, const Overlap*>& group) {
    if (group.size() > 0) {
        std::vector<int> cov(group.begin()->second->GetRead(id).len + 1, 0);
        const int redundance = - std::min<int>(500, cov.size()/10);

        for (const auto &ig : group) {
            const Overlap& o = *ig.second;
            if (IsReserved(o)) {

                bool added = false;
                // duplication
                auto dup0 = dup_groups_.find(id);
                if (dup0 != dup_groups_.end()) {
                    auto dup1 = dup0->second.find(ig.first);
                    if (dup1 != dup0->second.end()) {
                        // int start = cov.size() - 1;
                        // int end = 0;
                        // for (auto dup_o : dup1->second) {
                        //     auto& r = dup_o->GetRead(id);
                        //     start = std::min(start, r.start);
                        //     end = std::max(end, r.end);
                        // }
                        // if (std::max<size_t>(0, start-redundance) < std::min<size_t>(cov.size()-1, end+redundance)) {

                        //     cov[std::max<size_t>(0, start-redundance)] ++;
                        //     cov[std::min<size_t>(cov.size()-1, end+redundance)] --;
                        // }

                        for (auto dup_o : dup1->second) {
                            auto& r = dup_o->GetRead(id);
                            
                            if (std::max<size_t>(0, r.start-redundance) < std::min<size_t>(cov.size()-1, r.end+redundance)) {

                                cov[std::max<size_t>(0, r.start-redundance)] ++;
                                cov[std::min<size_t>(cov.size()-1, r.end+redundance)] --;
                            }
                        }

                        added = true;
                    }
                } 
                
                if (!added) {
                    
                    auto& r = o.GetRead(id);
                    if (std::max(0, r.start-redundance) < std::min(r.len, r.end+redundance)) {
                        cov[std::max(0, r.start-redundance)] ++;
                        cov[std::min(r.len, r.end+redundance)] --;
                    }
                    if ( opts_.debug_name == rd_store_.QueryNameById(id)) {
                        assert(std::max(0, r.start-redundance) < std::min(r.len, r.end+redundance));
                        printf("cov %s\n", OverlapStore::ToPafLine(o, StringPool::UnsafeNameId(rd_store_.GetStringPool())).c_str());
                    }
                }
            }
        }
        for (size_t i = 1; i < cov.size(); ++i) {
            cov[i] += cov[i - 1];
        }
        assert(cov.back() == 0);

        if ( opts_.debug_name == rd_store_.QueryNameById(id)) {
           for (size_t i = 0; i < cov.size(); ++i) {
               printf("cov %d\n", cov[i]);
           }
        }

        auto covtype = 0;//AnalyzeCoverageType(std::vector<int>(cov.begin()-redundance, cov.end()+redundance), rd_store_.QueryNameById(id) == "178924");
        if (covtype) {
            printf("cov_abn: %s\n", rd_store_.QueryNameById(id).c_str());
        }

        if (cov.size() > 3000) {
            read_infos_[id].cliff = CoverageConfidencePoints1(std::vector<int>(cov.begin()-redundance, cov.end()+redundance), rd_store_.QueryNameById(id) == opts_.debug_name);
            if (read_infos_[id].cliff[0] > 0) read_infos_[id].cliff[0] -= redundance;
            if (read_infos_[id].cliff[1] > 0) read_infos_[id].cliff[1] -= redundance;

            if (read_infos_[id].cliff[0] > 0) {
                if (read_infos_[id].cliff[0] >= std::min<int>(opts_.max_unreliable_length, opts_.max_unreliable_rate*cov.size())) {
                    read_infos_[id].cliff[0] = -1;
                }
            }
            
            if (read_infos_[id].cliff[1] > 0) {
                if ((int)cov.size() - read_infos_[id].cliff[1] >= std::min<int>(opts_.max_unreliable_length, opts_.max_unreliable_rate*cov.size())) {
                    read_infos_[id].cliff[1] = -1;
                }
            }
        }

        int inv = cov.size() / 3;
        auto a0 = std::accumulate(cov.begin(), cov.begin()+inv, 0) / inv;
        auto a1 = std::accumulate(cov.begin()+inv, cov.begin()+2*inv, 0) / inv;
        auto a2 = std::accumulate(cov.begin()+2*inv, cov.end(), 0) / inv;
        read_infos_[id].coverage = {a0, a1, a2};

        int oh = std::max(0, -redundance) + opts_.filter0.max_overhang;
        
        if ( opts_.debug_name == rd_store_.QueryNameById(id)) {
        LOG(INFO)("XXXX %d, %d, %d", oh, redundance, opts_.filter0.max_overhang);
        }
        auto c_minmax = std::minmax_element(oh + cov.begin(), cov.end() -1 - oh);
        return {covtype, *c_minmax.first, *c_minmax.second };

    } else {
        return { 0, 0, 0 };
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


// std::unordered_set<const Overlap*> AsmDataset::FindBestN(const std::pair<int, std::unordered_map<int, const Overlap*>> &g) const {
//     std::unordered_set<const Overlap*> keep;
//     std::vector<const Overlap*> left, right;
//     for (auto &i : g.second) {
//         const Overlap& o = *(i.second);

//         if (IsReserved(o)) {
//             auto loc = o.Location(i.first, 0);
//             if (loc == Overlap::Loc::Left) {
//                 left.push_back(&o);
//             } else {
//                 assert(loc == Overlap::Loc::Right);
//                 right.push_back(&o);
//             }
            
//         }
//     }

//     if (left.size() > (size_t)opts_.bestn) {
//         std::sort(left.begin(), left.end(), [](const Overlap* a, const Overlap *b){ return BetterAlignedLength(*a, *b); });

//         keep.insert(left.begin(), left.begin() + opts_.bestn);
//     }
//     else {
//         keep.insert(left.begin(), left.end());
//     }
//     if (right.size() > (size_t)opts_.bestn) {
//         std::sort(right.begin(), right.end(), [](const Overlap* a, const Overlap *b) { return BetterAlignedLength(*a, *b);} );

//         keep.insert(right.begin(), right.begin() + opts_.bestn);

//     }
//     else {
//         keep.insert(right.begin(), right.end());
//     }

//     return keep;
// }

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

    auto g = groups_.find(tid);
    if (g != groups_.end()) {
        for (auto &i : g->second) {
            auto &o = *i.second;
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

    auto g = groups_.find(tid);
    if (g != groups_.end()) {
        for (auto &i : g->second) {
            auto &o = *i.second;
            auto &qread = o.GetOtherRead(tid);
            nearby.insert(qread.id);
        }
    }
    return nearby;
}

std::unordered_set<const Overlap*> AsmDataset::GetExtendOverlaps(Seq::Id tid, int end) const {
    std::unordered_set<const Overlap*> extend;

    auto g = groups_.find(tid);
    if (g != groups_.end()) {
        for (auto &i : g->second) {
            auto &o = *i.second;
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
    auto g = groups_.find(tid);
    if (g != groups_.end()) {
        for (auto &i : g->second) {
            auto &o = *i.second;
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
                        auto dup0 = dup_groups_.find(tid);
                        if (dup0 != dup_groups_.end()) {
                            auto dup00 = dup0->second.find(qread.id);
                            if (dup00 != dup0->second.end()) {
                                DUMPER["test"]("extend_insert dupsize: %zd\n", dup00->second.size());
                                for (auto idup: dup00->second) {
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

    auto g = groups_.find(tid);
    if (g != groups_.end()) {
        for (auto &i : g->second) {
            auto &o = *i.second;
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

    groups_[new_ol->a_.id][new_ol->b_.id] = new_ol;
    groups_[new_ol->b_.id][new_ol->a_.id] = new_ol;
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

void AsmDataset::Check_Group() const {

    auto get_overlap_in_groups = [this](int a, int b) -> const Overlap* {
        auto it0 = groups_.find(a);
        if (it0 != groups_.end()) {
            auto it1 = it0->second.find(b);
            if (it1 != it0->second.end()) {
                return it1->second;
            }
        }
        return nullptr;
    };

    
    auto get_overlaps_in_dups = [this](int a, int b) -> const std::vector<const fsa::Overlap*>* {
        auto it0 = dup_groups_.find(a);
        if (it0 != dup_groups_.end()) {
            auto it1 = it0->second.find(b);
            if (it1 != it0->second.end()) {
                return &it1->second;
            }
        }
        return nullptr;
    };
    
    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);

        assert(get_overlap_in_groups(o.a_.id, o.b_.id) == get_overlap_in_groups(o.b_.id, o.a_.id));
        
        auto dup_a = get_overlaps_in_dups(o.a_.id, o.b_.id);
        auto dup_b = get_overlaps_in_dups(o.b_.id, o.a_.id);
        if (dup_a != nullptr && dup_b != nullptr) {
            assert(dup_a->size() == dup_b->size());
            assert(dup_a->size() > 1);
        } else {
            assert(dup_a == dup_b); // dup_a == nullptr && dup_b == nullptr
        }

    }
}

void AsmDataset::Dump() const {
    LOG(INFO)("Dump infos");
    
    DumpReadInfos(OutputPath("readinfos"), read_infos_); 
    DumpOverlaps(OutputPath("filter.m4a"));
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

    // if (of != nullptr) {
    //     for (size_t i=0; i < ol_store_.Size(); ++i) {
    //         const auto &o = ol_store_.Get(i);
    //         OlReason rs = GetOlReason(o);
    //         switch(rs.type) {    
    //         case OlReason::RS_FILTERED_READ:
    //             gzprintf(of, "%s %s %s %s %d\n", rd_store_.QueryNameById(o.a_.id).c_str(), rd_store_.QueryNameById(o.b_.id).c_str(), 
    //                 rs.ToString(), rd_store_.QueryNameById(rs.sub[0]).c_str(), rs.sub[1]);
    //             break;

    //         case OlReason::RS_SIMPLE:
    //         case OlReason::RS_DUPLICATE:
    //         case OlReason::RS_LOCAL:
    //         case OlReason::RS_CONSISTENCY:
    //         case OlReason::RS_CONSISTENCY1:
    //         case OlReason::RS_CONTIG:
    //         case OlReason::RS_UNKNOWN:
    //             gzprintf(of, "%s %s %s %d %d\n", rd_store_.QueryNameById(o.a_.id).c_str(), rd_store_.QueryNameById(o.b_.id).c_str(), 
    //                 rs.ToString(), rs.sub[0], rs.sub[1]);
    //             break;
    //         case OlReason::RS_OK:
    //         default:
    //             break;
    //         }
    //     }
    //     gzclose(of);
    // } else {
    //     LOG(ERROR)("Fail to open filterd reads file %s", fname.c_str());
    // }
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
                   << ri.cliff[0] << " " << ri.cliff[1] << "\n";
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

} // namespace fsa {
