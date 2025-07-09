#include "coverage_info.hpp"
#include "align/match_info.hpp"

namespace fsa {



void CoverageInfo::Merge(const MatchInfo &match, size_t offsize, double local_threashold, size_t max_clip, double wt) {

    size_t s = match.Start();
    for (size_t i = 0; i < match.Size(); ++i) {
        const auto& info = match.Get(i);

        assert(base_cov_[s+i].ref == info.ref);

        base_cov_[s+i].bases0[info.base] += wt;
        base_cov_[s+i].c[0] += 1;
        auto inssize = match.GetInssize(info.ins);
        if (inssize > 0) {
            base_cov_[s+i].bases0[5] += wt;
            base_cov_[s+i].inssize += inssize*wt;
        }
    }
        
    if (match.LClip() >= max_clip) {
        base_cov_[s].clips += wt;
    }
    if (match.RClip() >= max_clip) {
        base_cov_[match.End()-1].clips += wt;
    }

    auto regs = match.GetHighQualityRegions(offsize, local_threashold, max_clip, offsize);
    for (auto& r : regs) {
        LOG(INFO)("XXX reg %d %d", s+r[0], s+r[1]);
        for (size_t i = r[0]; i < r[1]; ++i) {
            const auto& info = match.Get(i);
            base_cov_[s+i].bases1[info.base] += wt;
            base_cov_[s+i].c[1] += 1;
            auto inssize = match.GetInssize(info.ins);
            if (inssize > 0) {
                base_cov_[s+i].bases1[5] += wt;
            }
        }        
    }

    // read 是否整段全部比对到组装结果上
    if (regs.size() == 1) {        
        if (regs[0][0] + s == match.GetOverlap()->b_.start && regs[0][1] + s == match.GetOverlap()->b_.end) {
            for (size_t i = regs[0][0]; i < regs[0][1]; ++i) {       
                const auto& info = match.Get(i);
                base_cov_[s+i].bases2[info.base] += wt;
                base_cov_[s+i].c[2] += 1;
                auto inssize = match.GetInssize(info.ins);
                if (inssize > 0) {
                    base_cov_[s+i].bases2[5] += wt;
                }
            }
        }
    }
}
 
void CoverageInfo::Scan() {
    for (size_t i = 0; i < base_cov_.size(); ++i) {
        auto &bc = base_cov_[i];
        auto mx = std::max_element(bc.bases1, bc.bases1+6);
        auto ss = std::accumulate(bc.bases1, bc.bases1+6, 0);
        bc.top = mx - bc.bases1;

        if (*mx > ss * 0.7 && bc.top != 5 && bc.clips == 0 && bc.inssize == 0) {
            bc.level = 1;
        } else {
            bc.level = 0;
        }
    }
}


void CoverageInfo::Stat() {
    std::vector<uint32_t> covs0;
    std::vector<uint32_t> covs1;
    std::vector<uint32_t> covs2;
    for (size_t i = 0; i < base_cov_.size(); ++i) {
        auto &bc = base_cov_[i];
        auto cc = bc.Coverage();
        covs0.push_back(cc[0]);
        covs1.push_back(cc[1]);
        covs2.push_back(cc[2]);
    }
    std::sort(covs0.begin(), covs0.end());
    auto iqr0 = covs0[covs0.size()*3/4] - covs0[covs0.size()/4];
    std::sort(covs1.begin(), covs1.end());
    auto iqr1 = covs1[covs1.size()*3/4] - covs1[covs1.size()/4];
    std::sort(covs2.begin(), covs2.end());
    auto iqr2 = covs2[covs2.size()*3/4] - covs2[covs0.size()/4];
    LOG(INFO)("Coverage iqr0: %zd %zd %zd", iqr0, iqr1, iqr2);


    double ava_cov0 = std::accumulate(covs0.begin() + covs0.size() / 4, covs0.begin() + covs0.size()*3/4, 0) * 2/ covs0.size();
    double ava_cov1 = std::accumulate(covs1.begin() + covs1.size() / 4, covs1.begin() + covs1.size()*3/4, 0) * 2/ covs1.size();
    double ava_cov2 = std::accumulate(covs2.begin() + covs2.size() / 4, covs2.begin() + covs2.size()*3/4, 0) * 2/ covs2.size();
    LOG(INFO)("Average Coverage: %.02f %.02f %.02f", ava_cov0, ava_cov1, ava_cov2);
    average_coverages_ = {ava_cov0, ava_cov1, ava_cov2};
}

void CoverageInfo::Dump(std::ofstream& of, const std::string& ctg_name) {
    for (size_t i = 0; i < base_cov_.size(); ++i) {
        const auto &c = base_cov_[i];
        of << ctg_name << " " << i << " " << c.ToString() << "\n";
    }
}


}