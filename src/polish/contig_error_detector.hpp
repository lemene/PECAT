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
        sprintf(buf, "%d (%d %d %d %d %d %d) (%d %d %d %d %d %d)  %d %d", ref, 
            abases[0], abases[1], abases[2], abases[3], abases[4], abases[5], 
            cbases[0], cbases[1], cbases[2], cbases[3], cbases[4], cbases[5], 
            inssize, clips);
        return buf;
    }

    uint8_t ref;
    uint32_t abases[6];
    uint32_t cbases[6];
    uint32_t inssize;
    uint16_t clips;        
};

class ContigRegion {
public:

};

struct ErrorRegion {
    size_t start;
    size_t end;
    int type;
};

class ContigErrorDetector {
public:
    ContigErrorDetector(Seq::Id tid, const PolDataset& ds);

    void Detect();
    void ComputeCoverage();
    void MergeCoverage(const MatchInfo& match);

    //void EvaluateQuality();

    std::vector<ContigRegion> Split();

    std::vector<ErrorRegion> DetectUnclearRegions(const std::vector<BaseCoverage>& ctg_cov);
    std::vector<ErrorRegion> DetectUncoveredRegions(const std::vector<BaseCoverage>& ctg_cov, size_t min_cov=0);
    std::vector<std::array<uint32_t, 2>> CalculateCoverage();
    void AnalyzeCoverage(const std::vector<std::array<uint32_t, 2>> &covs);
    std::vector<ErrorRegion> MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap=1000);
    
    bool CheckRegion(const ErrorRegion& reg, const std::vector<std::array<uint32_t,2>>& covs);
    void CheckRegions(const std::vector<ErrorRegion>& regs);
    void DumpCoverage(std::ofstream& of);
    void SaveErrors(std::ofstream& of);

protected:
    Seq::Id tid_;
    const PolDataset& dataset_;
    std::vector<BaseCoverage> ctg_cov_;
    std::vector<std::array<uint32_t, 2>> win_cov_;
    std::vector<ErrorRegion> errors_;

    const uint32_t WIN_SIZE = 400;
    const uint32_t STRIDE = 200;
    const uint32_t MIN_CLIP = 500;

    std::vector<MatchInfo> match_;
};

} // namespace fsa {
    
