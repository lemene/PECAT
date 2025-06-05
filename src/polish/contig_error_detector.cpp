#include "contig_error_detector.hpp"

#include "align/match_info.hpp"

namespace fsa {

ContigErrorDetector::ContigErrorDetector(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds) {
    size_t tlen = ds.seq_store_.GetSeqLength(tid);
    ctg_cov_.assign(tlen, BaseCoverage());
}

void ContigErrorDetector::Detect() {

    ComputeCoverage();
    //EvaluateQuality();
 
    win_cov_ = CalculateCoverage();

    auto errors = DetectUncoveredRegions(ctg_cov_);
    for (auto& reg : errors) {
        LOG(INFO)("Uncovered: %s:%zd-%zd", dataset_.QueryStringById(tid_).c_str(), reg.start, reg.end);
    }

    auto unclear = DetectUnclearRegions(ctg_cov_);
    for (auto& reg : MergeRegions(unclear)) {
        LOG(INFO)("Ubi: %zd-%zd %d", reg.start, reg.end);
    }
    errors.insert(errors.end(), unclear.begin(), unclear.end());

    std::sort(errors.begin(), errors.end(), [](const ErrorRegion& a, const ErrorRegion& b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });

    LOG(INFO)("eRROR %zd", errors.size());
    CheckRegions(MergeRegions(errors));
    for (auto reg : MergeRegions(errors)) {
        if (CheckRegion(reg, win_cov_)) {
            errors_.push_back(reg);
        }
    }

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
            const auto &ol = *ol_group.Get(i, j);
            const DnaSeq& query = seq_store.GetSeq(ol.a_.id);
            assert(query.Size() >= 2000);
            
            match_.push_back(MatchInfo(&ol, query, target));
            if ( match_.back().Identity() >= 0.70) {
                MergeCoverage(match_.back());
            }
            break;  // 

        }
    }
}


void ContigErrorDetector::MergeCoverage(const MatchInfo &match) {
    size_t s = match.Start();
    for (size_t i = 0; i < match.Size(); ++i) {
        const auto& info = match.Get(i);
        ctg_cov_[s+i].ref = info.ref;
        ctg_cov_[s+i].abases[info.base] ++;
        auto inssize = match.GetInssize(info.ins);
        if (inssize > 0) {
            ctg_cov_[s+i].abases[5] ++;
            ctg_cov_[s+i].inssize += inssize;
        }
    }
    // if (match.Start() < 1502600 && match.End()-1 < 4804600) {
    //     LOG(INFO)("xxx: %s, %zd %zd %zd %zd\n", dataset_.QueryStringById(match.GetOverlap()->a_.id).c_str(), match.LClip(), match.RClip(), s, match.End()-1);
    // }
        
    if (match.LClip() >= MIN_CLIP) {
        ctg_cov_[s].clips += 1;
    }
    if (match.RClip() >= MIN_CLIP) {
        ctg_cov_[match.End()-1].clips += 1;
    }

    auto regs = match.GetHighQualityRegions(1000, 0.3, MIN_CLIP, 1000);
    for (auto& r : regs) {
        for (size_t i = r[0]; i < r[1]; ++i) {
            const auto& info = match.Get(i);
            ctg_cov_[s+i].ref = info.ref;
            ctg_cov_[s+i].cbases[info.base] ++;
            auto inssize = match.GetInssize(info.ins);
            if (inssize > 0) {
                ctg_cov_[s+i].cbases[5] ++;
            }
        }
        if (r[0]+s < 1502600 && r[1]+s >1502800 ) {
             LOG(INFO)("xxx: %d, %s, %zd %zd %f\n", match.GetOverlap()->a_.id, dataset_.QueryStringById(match.GetOverlap()->a_.id).c_str(), r[0], r[1], match.MaxLocalDistance(1000));

        }

        // if (dataset_.QueryStringById(match.GetOverlap()->a_.id) == "a5a6d038-8f9f-956f-fe26-4a0ded8a29a6") {
        //     LOG(INFO)("xxx: %d %s, %zd %zd %f\n", match.GetOverlap()->a_.id, dataset_.QueryStringById(match.GetOverlap()->a_.id).c_str(), r[0], r[1], match.MaxLocalDistance(1000));
        // }

    }
}



class WinDivider {
public:
    WinDivider(size_t len, size_t win_size, size_t stride) : length_(len), win_size_(win_size), stride_(stride) {
        count_ = (std::max<int>(0, (int)length_ - (int)win_size) + stride_ -1) / stride_ + 1;
    }
    size_t Size() const {  return count_ ; }
    std::array<size_t, 2> Get(size_t i) const { 
        return {stride_*i, std::min<size_t>(length_, stride_*(i) + win_size_)}; 
    }


protected:
    size_t length_;
    size_t win_size_;
    size_t stride_;
    size_t count_;
};

    

auto ContigErrorDetector::DetectUncoveredRegions(const std::vector<BaseCoverage>& ctg_cov, size_t min_cov) -> std::vector<ErrorRegion> {
    std::vector<ErrorRegion> cands;

    WinDivider witr = WinDivider(ctg_cov_.size(), WIN_SIZE, STRIDE);
    for (size_t i = 0; i < witr.Size(); ++i) {
        auto win = witr.Get(i);
        LOG(INFO)("WIN: %zd-%zd", win[0], win[1]);

        std::vector<int> ss(win[1] - win[0], 0);

        for (size_t i = win[0]; i < win[1]; ++i) {
            const auto& c = ctg_cov[i];
            ss[i - win[0]] = std::accumulate(c.cbases, c.cbases+5, 0);
        }

        auto min = std::min_element(ss.begin(), ss.end());
        if (*min <= 0) {
            cands.push_back({win[0], win[1], 0});
        }

    }

    return cands;
}


auto ContigErrorDetector::DetectUnclearRegions(const std::vector<BaseCoverage>& ctg_cov) -> std::vector<ErrorRegion>  {
    std::vector<ErrorRegion> cands;

    WinDivider witr = WinDivider(ctg_cov_.size(), WIN_SIZE, STRIDE);
    for (size_t i = 0; i < witr.Size(); ++i) {
        auto win = witr.Get(i);


        size_t clips = 0;
        std::vector<size_t> cov0(win[1] - win[0], 0);
        std::vector<size_t> cov1(win[1] - win[0], 0);

        for (size_t i = win[0]; i < win[1]; ++i) {
            const auto& c = ctg_cov[i];
            clips += c.clips;
            cov0[i-win[0]] = std::accumulate(c.cbases, c.cbases+5, 0);
            cov1[i-win[0]] = std::accumulate(c.abases, c.abases+5, 0);
        }

        size_t sum_cov0 = std::accumulate(cov0.begin(), cov0.end(), 0);
        size_t sum_cov1 = std::accumulate(cov1.begin(), cov1.end(), 0);
        LOG(INFO)("UBI (%s) %zd-%zd %zd, %zd, %zd, %zd", dataset_.QueryStringById(tid_).c_str(), win[0], win[1], clips, sum_cov0, sum_cov1, sum_cov1- sum_cov0);

        if ((sum_cov1- sum_cov0)*1.0 > sum_cov0*0.2) {
            cands.push_back({win[0], win[1], 1}); 
        }
    }
    return cands;
}


std::vector<std::array<uint32_t, 2>> ContigErrorDetector::CalculateCoverage() {
    WinDivider witr = WinDivider(ctg_cov_.size(), WIN_SIZE, STRIDE);
    std::vector<std::array<uint32_t, 2>> covs(witr.Size(), {0, 0});

    for (size_t i = 0; i < witr.Size(); ++i) {
        auto win = witr.Get(i);

        size_t clips = 0;
        std::vector<size_t> cov0(win[1] - win[0], 0);
        std::vector<size_t> cov1(win[1] - win[0], 0);

        for (size_t i = win[0]; i < win[1]; ++i) {
            const auto& c = ctg_cov_[i];
            clips += c.clips;
            cov0[i-win[0]] = std::accumulate(c.cbases, c.cbases+5, 0);
            cov1[i-win[0]] = std::accumulate(c.abases, c.abases+5, 0);
        }
        covs[i] = {std::accumulate(cov0.begin(), cov0.end(), 0), std::accumulate(cov1.begin(), cov1.end(), 0) };

    }
    return covs;
}

void ContigErrorDetector::AnalyzeCoverage(const std::vector<std::array<uint32_t, 2>> &covs) {
    for (size_t i = 0 ; i < covs.size(); ++i) {
        
        if ((covs[i][1]- covs[i][0])*1.0 > covs[i][0]*0.2) {
            
        }

    }
 
}

void ContigErrorDetector::SaveErrors(std::ofstream& of) {

    auto ctg_name = dataset_.QueryStringById(tid_);
    for (auto w : errors_) {
        of << ctg_name << " " << w.start << " " << w.end << " " << w.type << "\n";    
    }
}


std::vector<ErrorRegion> ContigErrorDetector::MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) {

    // merge
    std::vector<ErrorRegion> merged;
    if (regs.size() > 0) {
        merged.push_back(regs[0]);

        for (size_t i = 1; i < regs.size(); ++i) {
            if (regs[i].start <= merged.back().end + max_gap) {
                assert(merged.back().end <= regs[i].end);
                merged.back().end = regs[i].end;
            } else {
                merged.push_back(regs[i]);
            }
        }
    }
    return merged;
}



// void ContigErrorDetector::EvaluateQuality() {
//     size_t ref_mis = 0;
//     size_t rd_mat = 0;
//     size_t rd_total = 0;
//     //for (const auto& c : ctg_cov_) {
//     for (size_t i = 0; i < ctg_cov_.size(); ++i) {
//         const auto& c = ctg_cov_[i];
//         auto m = std::max_element(c.bases, c.bases+6);
//         //LOG(INFO)("S(%zd): %s | %zd %zd", i, c.ToString().c_str(), m-c.bases, *m);

//         if ((m-c.bases) != c.ref) {
//             ref_mis ++;
//         }

//         rd_total += c.bases[0] + c.bases[1] + c.bases[2] + c.bases[3] + c.bases[5];

//         switch ((m-c.bases)) {
//         case 0:
//         case 1:
//         case 2:
//         case 3:
//             rd_mat += *m;
//             break;
//         case 4:
//             // 
//             break;
//         case 5:
//             rd_mat += c.inssize;
//             break;
//         default:
//             LOG(INFO)("XXX %zd", (m-c.bases));
//             assert(!"never coming here");
//         }
    
//     }
//     LOG(INFO)("%zd, %zd, %zd, %zd", ref_mis, ctg_cov_.size(), rd_mat, rd_total);
//     LOG(INFO)("%.02f, %.02f", ref_mis*1.0/ctg_cov_.size(), rd_mat*1.0 / rd_total);
// }

auto ContigErrorDetector::Split() -> std::vector<ContigRegion> {
    std::vector<ContigRegion> regs;

    return regs;
}

bool ContigErrorDetector::CheckRegion(const ErrorRegion& reg, const std::vector<std::array<uint32_t,2>>& covs) {
    size_t count = 0;
    auto s = reg.start / STRIDE;
    auto e = (reg.end - WIN_SIZE + 1) / STRIDE;
    if (s > 1 && e + 2< covs.size()) {
        LOG(INFO)("check reg start %zd %zd %zd %zd, %zd %zd: ", reg.start, reg.end, covs[s-1][0], covs[s-1][1], covs[s][0], covs[s][1]);
        LOG(INFO)("check reg end   %zd %zd, %zd %zd: ", covs[e][0], covs[e][1], covs[e+1][0], covs[e+1][1]);
        // if (!(covs[s-2][1] < covs[s][1] && covs[s][1] - covs[s-2][1] > covs[s-2][1] *0.2 && 
        //     covs[e][1] > covs[e+2][1] &&  covs[e][1] - covs[e+2][1] > covs[e+2][1] *0.2)) {
        // //if (std::abs((int)covs[e+1][0] - (int)covs[s-1][0]) > std::min<size_t>(covs[e+1][0] , covs[s-1][0])*0.2) {
        //     return true;
        // }
    }
    for (auto& m : match_) {
        if (m.Start() + 100 < reg.start && m.End() > reg.end + 100 && 
            m.MaxLocalDistance(1000) < 0.2 && m.LClip() < 1000 && m.RClip() < 1000) {
            count++;
        }
    }
    LOG(INFO)("check reg: %zd-%zd %zd", reg.start, reg.end, count);
    return count < 2;
}

void ContigErrorDetector::CheckRegions(const std::vector<ErrorRegion>& regs) {
    // regs 必须排过序

    const std::vector<ErrorRegion> checked;

    auto covs = CalculateCoverage();

    for (const auto &r : regs) {
        auto s = r.start / STRIDE;
        auto e = (r.end - WIN_SIZE + 1)  / STRIDE;

        if (s > 1) {
            LOG(INFO)("xxxx , (%zd,%zd)", covs[s-1][0], covs[s-1][1]);
        }
        LOG(INFO)("%zd %zd, (%zd,%zd) (%zd,%zd)", r.start, r.end, covs[s][0], covs[s][1],  covs[e][0], covs[e][1]);

        if (e + 1 < covs.size()) {
            LOG(INFO)("yyyy , (%zd,%zd)", covs[e+1][0], covs[e+1][1]);
        }
    }
}

void ContigErrorDetector::DumpCoverage(std::ofstream& of) {
    const std::string& ctg_name = dataset_.QueryStringById(tid_);
    for (size_t i = 0; i < ctg_cov_.size(); ++i) {
        const auto &c = ctg_cov_[i];
        of << ctg_name << " " << i << " " << c.ToString() << "\n";
    }
}

}