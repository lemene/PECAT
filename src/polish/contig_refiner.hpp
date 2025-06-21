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

class ContigRefiner {
public:
    struct Segment {

    };
public:
    ContigRefiner();
    std::string Consensus();

protected:
    
    std::vector<Segment> segs;
};
}