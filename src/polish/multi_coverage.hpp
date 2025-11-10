#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "pol_dataset.hpp"
#include "align/match_info.hpp"

namespace fsa {


class MultiCoverage {
public:

public:
    MultiCoverage(const DnaSeq& ctg) : covs_(5 + 3), win_size_(1000) {
        for (size_t i = 0; i < covs_.size(); ++i) {
            covs_[i].assign(ctg.Size(), 0.0);
        }
    }
    void Merge(const MatchInfo &match, double wt);

    void Dump(std::ofstream& of);

    enum CovType {
        ALL = 0,
        ALL_WT = 1,
        MATCHED = 2,
        MATCHED_WT = 3,
        THRESHOLD = 4,
    };
protected:
    std::vector<std::vector<double>> covs_;
    std::vector<double> locaL_thresholds_ { 0.85, 0.90 };
    size_t win_size_;
};
}