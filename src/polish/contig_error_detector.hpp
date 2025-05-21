#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "pol_dataset.hpp"

namespace fsa {

struct BaseCoverage {
    std::string ToString() const {
        char buf[1024];
        sprintf(buf, "%d (%d %d %d %d) %d %d %d", ref, bases[0], bases[0], bases[0], bases[0], del, ins, inssize);
        return buf;
    }
    uint8_t ref;
    uint8_t bases[4];
    uint8_t del;
    uint8_t ins;
    uint32_t inssize;
};

class ContigErrorDetector {
public:
    ContigErrorDetector(Seq::Id tid, const PolDataset& ds);
    void Detect();
    void ComputeCoverage();
    std::vector<size_t> GetBigInserts(size_t start, size_t end);
protected:
    Seq::Id tid_;
    const PolDataset& dataset_;
    std::vector<BaseCoverage> ctg_cov_;
};

} // namespace fsa {
    
