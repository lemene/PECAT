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
        uint32_t c0;
        uint32_t c1;
        uint32_t min_c;
        uint32_t ins;
        uint32_t inssize;
        uint32_t clips;
        std::string ToString() const {
            char buf[1024];
            sprintf(buf, "%d %d %d %d %d", c0, c1, min_c, ins, inssize);
            return buf;
        }
        uint8_t type;
    };

public:
    WindowSlider(CoverageInfo &cov_info, uint32_t wsize, uint32_t stride)
     : cov_info_(cov_info), win_size_(wsize), stride_(stride) {

    }
    
    std::array<size_t,2> Window2Region(size_t i) { return {i*stride_, std::min(i*stride_+win_size_, cov_info_.Size())}; }
    std::array<size_t,2> Region2Window(size_t i) { return {i*stride_, std::min(i*stride_+win_size_, cov_info_.Size())}; }
    std::vector<ErrorRegion> DetectErrorRegions(size_t max_gap);
    std::vector<ErrorRegion> DetectErrorRegions1(size_t max_gap);
    std::vector<ErrorRegion> DetectSimpleRegions() const;
    std::vector<ErrorRegion> MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) const ;
    std::vector<ErrorRegion> MergeRegions2(const std::vector<ErrorRegion> &regs)const ;
    void Flush();
    void Dump(std::ofstream &of, const std::string& ctg_name);
protected:
protected:
    CoverageInfo &cov_info_;
    uint32_t win_size_;
    uint32_t stride_;
    
    std::vector<WinInfo> win_cov_;
};
}