#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "pol_dataset.hpp"
#include "align/match_info.hpp"

#include "window_slider.hpp"
#include "contig_fragment.hpp"

namespace fsa {

class ContigRegion{};


class ContigAnalyzer {
public:
    struct Segment {
        uint32_t start;
        uint32_t end;
        uint8_t type;
    };
public:
    ContigAnalyzer(Seq::Id tid, const PolDataset& ds);

    void Detect();
    void ComputeCoverage();
    //void EvaluateQuality();

    std::vector<ErrorRegion> MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap=1000);
    std::vector<ErrorRegion> MergeRegions2(const std::vector<ErrorRegion> &regs, size_t max_gap=1000);
    
    bool CheckRegion(const ErrorRegion& reg);
    void DetectErrors();
    void SaveErrors(std::ofstream& of);
    
    const std::string& Name() const { return dataset_.QueryStringById(tid_); }
    std::vector<ContigFragment> Split();

    void DumpCoverage(std::ofstream& of);
    void DumpWindow(std::ofstream& of);
    void DumpMatch(std::ofstream& of);

protected:
    Seq::Id tid_;
    const PolDataset& dataset_;

    const uint32_t WIN_SIZE = 400;
    const uint32_t STRIDE = 200;
    const uint32_t MIN_CLIP = 500;

    std::vector<MatchInfo> match_;
    CoverageInfo cov_info_;
    WindowSlider win_slider_;
    std::vector<ErrorRegion> errors_;
    std::vector<Segment> segs_;
    double max_local_distance_threshold_ {0.0};
};

} // namespace fsa {
    
