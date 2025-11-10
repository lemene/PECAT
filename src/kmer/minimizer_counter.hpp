#pragma once

#include "../sequence.hpp"
#include "kmer_counter.hpp"
#include "ranked_kmers.hpp"

namespace fsa {

struct Minimizer {
    uint64_t dir : 1; // 0: forward, 1: reverse
    uint64_t pos : 31; // position in the sequence
    uint64_t rid : 32; // 31 bits for read id
    uint64_t rank : 4;
    uint64_t hash : 60; // hash value of the minimizer
    uint64_t kmer;
};

class MinimizerCounter {
public:
    MinimizerCounter(size_t k, size_t w) : w_(w), k_(k), kmc_(k) {}

    std::vector<Minimizer> Count(const DnaSeq& seq, const RankedKmers& rkmers);
    std::vector<Minimizer> Count(const DnaSeq& seq, const RankedKmers* rkmers) {
        return rkmers != nullptr ? Count(seq, *rkmers) : Count(seq);
    }
    std::vector<Minimizer> Count(const DnaSeq& seq);

protected:
    size_t w_;
    size_t k_;
    KmerCounter kmc_;
};
        
}
