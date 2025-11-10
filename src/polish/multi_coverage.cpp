#include "multi_coverage.hpp"

#include "align/match_info.hpp"

namespace fsa {


void MultiCoverage::Merge(const MatchInfo &match, double wt) {

    size_t s = match.Start();
    for (size_t i = 0; i < match.Size(); ++i) {
        covs_[CovType::ALL][s+i] += 1;
        covs_[CovType::ALL_WT][s+i] += wt;
        const auto& info = match.Get(i);
        if (info.ref == info.base) {
            covs_[CovType::MATCHED][s+i] += 1;
            covs_[CovType::MATCHED_WT][s+i] += 1;
        }
    }
    size_t half_win_size = 500;
    if (match.Size() >= half_win_size*2+1) {
        auto local_idents = match.LocalIdentity(half_win_size*2+1);
        for (size_t i = 0; i < local_idents.size(); ++i) {
            for (size_t t = 0; t < locaL_thresholds_.size(); ++t) {
                if (local_idents[i] >= locaL_thresholds_[t]) {
                    assert(half_win_size + i < match.Size());
                    covs_[CovType::THRESHOLD+t][s + half_win_size + i] += wt;
                }
            }
        }
    }
}

void MultiCoverage::Dump(std::ofstream& of) {
    // of << "#Pos\tAll\tAll_Wt\t\Matched\tMatched_Wt";
    // for (size_t t = 0; t < locaL_thresholds_.size(); ++t) {
    //     of << "\tThreshold_" << locaL_thresholds_[t];
    // }
    // of << "\n";
    for (size_t i = 0; i < covs_[0].size(); ++i) {
        of << i << "\t" << covs_[CovType::ALL][i] << "\t" << covs_[CovType::ALL_WT][i] << "\t"
           << covs_[CovType::MATCHED][i];

        for (size_t t = 0; t < locaL_thresholds_.size(); ++t) {
            of << "\t" << covs_[CovType::THRESHOLD+t][i];
        }
        of << "\n";
    }
}


}