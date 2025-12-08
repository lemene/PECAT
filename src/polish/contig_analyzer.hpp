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

    void ComputeCoverage(size_t thread_size=1);

    std::vector<ErrorRegion> MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap=1000);

    void SaveErrors(std::ofstream& of);

    const std::string& Name() const { return dataset_.QueryStringById(tid_); }
    

    /** Save infomations */
    void DumpMultiCoverage(std::ofstream &of);
    Seq::Id GetId() const { return tid_; }

protected:
    Seq::Id tid_;
    const PolDataset& dataset_;

    const uint32_t WIN_SIZE = 400;
    const uint32_t STRIDE = 200;
    const uint32_t MIN_CLIP = 500;


    MultiCoverage multi_cov_;
};

} // namespace fsa {
    
