#pragma once
#include <vector>
#include <string>

#include "contig_fragment.hpp"


namespace fsa {

class ContigGraph {
public:
    ContigGraph() = default;
    ContigGraph(const ContigGraph&) = delete;
    ContigGraph& operator=(const ContigGraph&) = delete;

    void AddFragment(const std::vector<ContigFragment>& frags) {
        fragments_.insert(fragments_.end(), frags.begin(), frags.end());
    }

    void BuildGraph();

protected:
    std::vector<ContigFragment> fragments_;

    
    std::vector<std::string> sequences_; // Store sequences of fragments
    std::vector<std::string> qualities_; // Store qualities of fragments
    std::vector<std::pair<size_t, size_t>> edges_; // Store edges between fragments (start, end)
};

}