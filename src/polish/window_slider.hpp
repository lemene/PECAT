#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "pol_dataset.hpp"
#include "align/match_info.hpp"

#include "coverage_info.hpp"
namespace fsa {


struct ErrorRegion {
    size_t start;
    size_t end;
    int type;
};


class WinDivider {
public:
    WinDivider(size_t len, size_t win_size, size_t stride) : length_(len), win_size_(win_size), stride_(stride) {
        count_ = (std::max<int>(0, (int)length_ - (int)win_size) + stride_ -1) / stride_ + 1;
    }
    size_t Size() const {  return count_ ; }
    std::array<size_t, 2> Get(size_t i) const { 
        return {stride_*i, std::min<size_t>(length_, stride_*(i) + win_size_)}; 
    }

protected:
    size_t length_;
    size_t win_size_;
    size_t stride_;
    size_t count_;
};
    
        
    

class WindowSlider {
public:
    /** Sliding window information */
    struct WinInfo {
        double c0;
        double c1;
        double c2;
        uint32_t count[3];
        uint32_t min_c[3];
        uint32_t max_c[3];
        uint32_t ins;
        uint32_t inssize;
        uint32_t clips;
        std::string ToString() const {
            char buf[1024];
            sprintf(buf, "%.02f %.02f %.02f (%d %d) %d %d %d  %d %d %d", c0, c1, c2, min_c[1], max_c[1], ins, inssize, clips, count[0], count[1], count[2]);
            return buf;
        }
        uint8_t type;
    };

public:
    WindowSlider(CoverageInfo &cov_info, uint32_t wsize, uint32_t stride)
     : cov_info_(cov_info), win_size_(wsize), stride_(stride) {
    }
    
    std::array<size_t,2> Window2Region(size_t i) { return {i*stride_, std::min(i*stride_+win_size_, cov_info_.Size())}; }
    std::array<size_t,2> Region2Window(const std::array<size_t,2>& r) const { 
        return {r[0] / stride_, (r[1] > win_size_ ? r[1] - win_size_ : 0) / stride_};
    }
    std::array<size_t,2> Region2Window(const ErrorRegion& reg) const { return Region2Window(std::array<size_t,2> ({reg.start, reg.end})); }
    std::vector<ErrorRegion> DetectErrorRegions(size_t max_gap, const std::array<double,3>& ave_covs);
    std::vector<ErrorRegion> DetectErrorRegions2(size_t max_gap, const std::array<double,3>& ave_covs);

    std::vector<ErrorRegion> MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) const ;
    std::vector<ErrorRegion> MergeRegions2(const std::vector<ErrorRegion> &regs)const ;
    void Flush();
    void Dump(std::ofstream &of, const std::string& ctg_name);
    double SurroundingCoverage(const ErrorRegion& reg, size_t inv=10);
    std::array<double,2> ComputeCoverageThresholds(int type);
    std::vector<double> AveCoverages(size_t winnum);
    std::vector<double> StdCoverages(size_t winnum, const std::vector<double>& ave);
    std::vector<std::pair<double,size_t>> MaxCoverages(size_t winnum);
    std::vector<std::pair<double,size_t>> MinCoverages(size_t winnum);

    bool HasBreakpoint(size_t s, size_t e) const;
    bool HasBreakpoint(const std::array<size_t, 2> &w) const { return HasBreakpoint(w[0], w[1]); }
    bool HasAlternate(size_t s, size_t e) const;
    bool HasAlternate(const std::array<size_t,2> &w) const { return HasAlternate(w[0], w[1]); }
    
protected:
    CoverageInfo &cov_info_;
    uint32_t win_size_;
    uint32_t stride_;
    
    std::vector<WinInfo> win_cov_;
};
}