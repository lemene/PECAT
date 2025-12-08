#include "multi_coverage.hpp"

#include "align/match_info.hpp"

namespace fsa {


void MultiCoverage::Merge(const MatchInfo &match, double wt) {

    size_t s = match.Start();
    for (size_t i = 0; i < match.Size(); ++i) {
        base_covs_[s+i].all += 1;
        base_covs_[s+i].wt_all += wt;
        const auto& info = match.Get(i);
        if (info.ref == info.base) {
            base_covs_[s+i].matched += 1;
        }
    }
    for (size_t t = 0; t < 3; ++t) {
        size_t half_win_size = 500 * (t+1);
        if (match.Size() >= half_win_size*2+1) {
            auto local_idents = match.LocalIdentity(half_win_size*2+1);
            for (size_t i = 0; i < local_idents.size(); ++i) {
                //for (size_t t = 0; t < quality_threshods_.size(); ++t) {
                    if (local_idents[i] >= quality_threshods_[1]) {
                        assert(half_win_size + i < match.Size());
                        base_covs_[s+i+half_win_size].quals[t] += wt;
                    }
                //}
            }
        }

    }
}

std::vector<MultiCoverage::BaseCov> MultiCoverage::ToCov(const MatchInfo &match, double wt)  const {
    std::vector<BaseCov> covs(match.End() - match.Start());
    
    for (size_t i = 0; i < match.Size(); ++i) {
        covs[i].all += 1;
        covs[i].wt_all += wt;
        const auto& info = match.Get(i);
        if (info.ref == info.base) {
            covs[i].matched += 1;
        }
    }
    for (size_t t = 0; t < 3; ++t) {
        size_t half_win_size = 500 * (t+1);
        if (match.Size() >= half_win_size*2+1) {
            auto local_idents = match.LocalIdentity(half_win_size*2+1);
            for (size_t i = 0; i < local_idents.size(); ++i) {
                //for (size_t t = 0; t < quality_threshods_.size(); ++t) {
                    if (local_idents[i] >= quality_threshods_[1]) {
                        assert(half_win_size + i < match.Size());
                        covs[i+half_win_size].quals[t] += wt;
                    }
                //}
            }
        }

    }
    return covs;
}

void MultiCoverage::Merge(std::vector<BaseCov> &covs, size_t s) {
    for (size_t i = 0; i < covs.size(); ++i) {
        base_covs_[s+i].all += covs[i].all;
        base_covs_[s+i].matched += covs[i].matched;
        base_covs_[s+i].inserts += covs[i].inserts;
        base_covs_[s+i].deletions += covs[i].deletions;
        base_covs_[s+i].wt_all += covs[i].wt_all;
        for (size_t t = 0; t < 3; ++t) {
            base_covs_[s+i].quals[t] += covs[i].quals[t];
        }
    }
}

void MultiCoverage::Dump(std::ofstream& of) {
    for (size_t i = 0; i < base_covs_.size(); ++i) {
        base_covs_[i].Dump(of, i);
    }
}


}