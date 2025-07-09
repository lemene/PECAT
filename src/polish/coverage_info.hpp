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

    
struct BaseCoverage {
    std::string ToString() const {
        char buf[1024];
        sprintf(buf, "%d (all: %zd %.02f %.02f %.02f %.02f %.02f %.02f ) \
(seg: %zd %.02f %.02f %.02f %.02f %.02f %.02f ) \
(hq: %zd %.02f %.02f %.02f %.02f %.02f %.02f ) %.02f %.02f %zd", ref, 
            c[0], bases0[0], bases0[1], bases0[2], bases0[3], bases0[4], bases0[5], 
            c[1], bases1[0], bases1[1], bases1[2], bases1[3], bases1[4], bases1[5], 
            c[2], bases2[0], bases2[1], bases2[2], bases2[3], bases2[4], bases2[5], 
            inssize, clips, top);
        return buf;
    }

    std::array<uint32_t, 3> Coverage() const {
        return {std::accumulate(bases0, bases0+5, 0), std::accumulate(bases1, bases1+5, 0), std::accumulate(bases2, bases2+5, 0)};
    }

    uint8_t ref;
    double bases0[6];
    double bases1[6];
    double bases2[6];
    uint16_t c[3] = {0, 0, 0}; // count of bases in bases1, bases2, and bases0
    double inssize;
    double clips;     
    uint8_t top;   
    uint8_t level;
    
};

class CoverageInfo {
public:
    /** Sliding window information */
    class Base {

    };

public:
    CoverageInfo(const DnaSeq& ctg)
     : base_cov_(ctg.Size(), BaseCoverage()) {
        for (size_t i = 0; i < ctg.Size(); ++i) {
            base_cov_[i].ref = ctg[i];
        }
    }

    size_t Size() const { return base_cov_.size(); }
    const BaseCoverage& Get(size_t i) const { return base_cov_[i]; }
    void Merge(const class MatchInfo &match, size_t offsize, double local_threashold, size_t max_clip, double wt);
    uint8_t GetBestChoice(size_t i) { return base_cov_[i].top; }
    const std::array<double, 3> & AvarageCoverage() const { return average_coverages_; }
    void Scan();
    void Stat();
    void Dump(std::ofstream& of, const std::string& ctg_name);
protected:
    std::vector<BaseCoverage> base_cov_;
    std::array<double, 3> average_coverages_;
};
}