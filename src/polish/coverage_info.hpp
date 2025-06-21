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
        sprintf(buf, "%d (%d %d %d %d %d %d) (%d %d %d %d %d %d)  %d %d %d", ref, 
            abases[0], abases[1], abases[2], abases[3], abases[4], abases[5], 
            cbases[0], cbases[1], cbases[2], cbases[3], cbases[4], cbases[5], 
            inssize, clips, top);
        return buf;
    }

    uint8_t ref;
    uint32_t abases[6];
    uint32_t bbases[6];
    uint32_t cbases[6];
    uint32_t inssize;
    uint16_t clips;     
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
    void Merge(const class MatchInfo &match);
    uint8_t GetBestChoice(size_t i) { return base_cov_[i].top; }
    void Scan();
    void Dump(std::ofstream& of, const std::string& ctg_name);
protected:
    std::vector<BaseCoverage> base_cov_;
};
}