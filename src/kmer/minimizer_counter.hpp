#pragma once

#include "../sequence.hpp"

namespace fsa {


class MinimizerCounter {
public:
    MinimizerCounter(size_t w, size_t k) : w_(w), k_(k) {}

    void Count(const DnaSeq& seq) {
        
    }

protected:
    size_t w_;
    size_t k_;
};
        
}
