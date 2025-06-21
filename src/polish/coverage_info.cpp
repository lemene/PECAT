#include "coverage_info.hpp"
#include "align/match_info.hpp"

namespace fsa {



void CoverageInfo::Merge(const MatchInfo &match) {
    int MIN_CLIP = 500;
    size_t s = match.Start();
    for (size_t i = 0; i < match.Size(); ++i) {
        const auto& info = match.Get(i);
        base_cov_[s+i].ref = info.ref;
        base_cov_[s+i].abases[info.base] ++;
        auto inssize = match.GetInssize(info.ins);
        if (inssize > 0) {
            base_cov_[s+i].abases[5] ++;
            base_cov_[s+i].inssize += inssize;
        }
    }
    if (match.Start() < 984800 && match.End()-1 > 984800) {
        LOG(INFO)("XXXall: %d, %zd %zd %zd %zd\n", match.GetOverlap()->a_.id, match.LClip(), match.RClip(), s, match.End()-1);
    }
        
    if (match.LClip() >= MIN_CLIP) {
        base_cov_[s].clips += 1;
    }
    if (match.RClip() >= MIN_CLIP) {
        base_cov_[match.End()-1].clips += 1;
    }

    auto regs = match.GetHighQualityRegions(1000, 0.2, MIN_CLIP, 1000);
    for (auto& r : regs) {
        LOG(INFO)("XXX reg %d %d", s+r[0], s+r[1]);
        for (size_t i = r[0]; i < r[1]; ++i) {
            const auto& info = match.Get(i);
            base_cov_[s+i].cbases[info.base] ++;
            auto inssize = match.GetInssize(info.ins);
            if (inssize > 0) {
                base_cov_[s+i].cbases[5] ++;
            }
        }        
    }
    
    // regs = match.GetHighQualityRegions(1000, 0.2, MIN_CLIP, 1000000);
    // for (auto& r : regs) {
    //     LOG(INFO)("XXX reg %d %d", s+r[0], s+r[1]);
    //     for (size_t i = r[0]; i < r[1]; ++i) {
    //         const auto& info = match.Get(i);
    //         base_cov_[s+i].ref = info.ref;
    //         base_cov_[s+i].cbases[info.base] ++;
    //         auto inssize = match.GetInssize(info.ins);
    //         if (inssize > 0) {
    //             base_cov_[s+i].cbases[5] ++;
    //         }
    //     }

    //     if (s+r[0] < 686444 && s+r[1] > 686444) {
    //         LOG(INFO)("XXXerr %d", match.GetOverlap()->a_.id);
    //     }
        
    // }
}
 
void CoverageInfo::Scan() {
    for (size_t i = 0; i < base_cov_.size(); ++i) {
        auto &bc = base_cov_[i];
        auto mx = std::max_element(bc.cbases, bc.cbases+6);
        auto ss = std::accumulate(bc.cbases, bc.cbases+6, 0);
        bc.top = mx - bc.cbases;

        if (*mx > ss * 0.7 && bc.top != 5 && bc.clips == 0 && bc.inssize == 0) {
            bc.level = 1;
        } else {
            bc.level = 0;
        }
    }


}

void CoverageInfo::Dump(std::ofstream& of, const std::string& ctg_name) {
    for (size_t i = 0; i < base_cov_.size(); ++i) {
        const auto &c = base_cov_[i];
        of << ctg_name << " " << i << " " << c.ToString() << "\n";
    }
}


}