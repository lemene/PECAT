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
#include "multi_coverage.hpp"

namespace fsa {

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

    bool CheckRegion(const ErrorRegion& reg);
    void DetectErrors();
    void SaveErrors(std::ofstream& of);

    const std::string& Name() const { return dataset_.QueryStringById(tid_); }
    
    size_t FirstMatch(size_t pos);
    size_t LastMatch(size_t pos);
    std::vector<ContigFragment> Split();
    std::string Polish(size_t s, size_t e);
    std::vector<const MatchInfo*> GetCoverage(size_t pos, int flank);

    /** Save infomations */
    void DumpCoverage(std::ofstream& of);
    void DumpWindow(std::ofstream& of);
    void DumpMatch(std::ofstream& of);
    void DumpMultiCoverage(std::ofstream &of);
    Seq::Id GetId() const { return tid_; }

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

    MultiCoverage multi_cov_;
};

} // namespace fsa {
    
