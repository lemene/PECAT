#pragma once

#include <string>
#include "sequence.hpp"

namespace fsa {

class ContigAnalyzer;
class WindowSlider;

class ContigFragment {
public:
    ContigFragment(ContigAnalyzer* ctg_analyzer, size_t start, size_t end);

    Seq::Id ContigId() const;
    std::string Polish() const ;
    ContigAnalyzer* Analyzer() const { return ctg_analyzer_; }
    size_t Start() const { return start_; }
    size_t End() const { return end_; }
protected:
    ContigAnalyzer *ctg_analyzer_ {nullptr};
    WindowSlider *win_slider_ {nullptr};

    size_t start_ {0};
    size_t end_ {0};
    uint8_t type_ {0}; // 0: normal, 1: error
};

}