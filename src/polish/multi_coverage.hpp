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
    struct BaseCov {
        double all {0.0};
        double wt_all {0.0};
        double quals[3] {0.0, 0.0, 0.0};  // for different quality thresholds

        void Dump(std::ofstream& of, size_t pos) {
            of << pos << "\t" << all << "\t" << wt_all;
            for (size_t i = 0; i < 3; ++i) {
                of << "\t" << quals[i];
            }
            of << "\n";
        }
    };
public:
    MultiCoverage(const DnaSeq& ctg, double qual_median, double qual_mad) : base_covs_(ctg.Size()), win_size_(1000) {
        assert(quality_threshods_.size() == 3);
        quality_threshods_[0] = qual_median - 6 * 1.4826 * qual_mad;
        quality_threshods_[1] = qual_median - 4 * 1.4826 * qual_mad;
        quality_threshods_[2] = qual_median - 2 * 1.4826 * qual_mad;
    }
    void Merge(const MatchInfo &match, double wt);

    std::vector<BaseCov> ToCov(const MatchInfo &match, double wt) const;
    void Merge(std::vector<BaseCov> &covs, size_t s);

    void Dump(std::ofstream& of);
protected:
    std::vector<BaseCov> base_covs_;
    std::vector<double> quality_threshods_ { 0.85, 0.90, 0.95 };
    size_t win_size_;
};
}