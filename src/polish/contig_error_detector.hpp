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
        sprintf(buf, "%d (%d %d %d %d) %d %d %d", ref, bases[0], bases[1], bases[2], bases[3], bases[4], bases[5], inssize);
        return buf;
    }
    void Merge(const BaseCoverage &c) {
        assert(ref == c.ref);
        for (size_t i = 0; i < sizeof(bases); ++i) bases[i] += c.bases[i];
        inssize += c.inssize;
        clips += c.clips;
    }

    uint8_t ref;
    uint8_t bases[6];       // A C G T DEL INS
    uint32_t inssize;
    uint16_t clips;         
};

class ContigErrorDetector {
public:
    ContigErrorDetector(Seq::Id tid, const PolDataset& ds);
    void Detect();
    void ComputeCoverage();
    std::vector<BaseCoverage> ComputeCoverage(const Overlap& ol);
    void MergeCoverage(const std::vector<BaseCoverage>& cov, const Overlap &ol);
    void CollectCandidates();
    void VerifyCandidates(const std::vector<std::array<size_t, 2>> &merged);
    std::vector<size_t> GetBigInserts(size_t start, size_t end);
    void EvaluateQuality();

protected:
    Seq::Id tid_;
    const PolDataset& dataset_;
    std::vector<BaseCoverage> ctg_cov_;
};

} // namespace fsa {
    
