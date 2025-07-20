#pragma once

#include "kmer.hpp"

namespace fsa {



class KmerCounter { 
public:
    KmerCounter(size_t k) : k_(k) {
        shift1 = 2 * (k - 1);
        mask = (1ULL<<2*k) - 1;
    }

    std::vector<std::array<KmerId, 2>> CountAll(const DnaSeq& seq) {
        std::vector<std::array<KmerId, 2>> kmers;

        if (seq.Size() >= k_) {    
            std::array<KmerId, 2> kmer = {0, 0};
            
            size_t index = 0;
            for (index = 0; index < k_-1; ++index) {
                auto c = seq[index];
                kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
            }
            for (; index < seq.Size(); index++) {
                auto c = seq[index];
                kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
                kmers.push_back(kmer);
            }
        }
        return kmers;
    }

protected:
    size_t k_;
    uint64_t shift1;
    uint64_t mask;
};

}