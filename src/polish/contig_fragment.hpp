#pragma once

#include <string>
#include "sequence.hpp"

namespace fsa {

class ContigAnalyzer;

class ContigFragment {
public:
    ContigFragment(ContigAnalyzer* ctg_analyzer, size_t start, size_t end);

    Seq::Id ContigId() const;
protected:
    ContigAnalyzer *ctg_analyzer_ {nullptr};
    size_t start_ {0};
    size_t end_ {0};
    uint8_t type_ {0}; // 0: normal, 1: error
};

}