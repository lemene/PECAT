#pragma once

#include "kmer.hpp"

namespace fsa {



class KmerCounter { 
public:
    KmerCounter(size_t k) : k_(k) {
        shift1 = 2 * (k - 1);
        mask = (1ULL<<2*k) - 1;
    }

    std::vector<KmerId> CountCanon(const DnaSeq& seq) {
        std::vector<KmerId> kmers;
        kmers.reserve(seq.Size() - k_ + 1); // Reserve space for k-mers based on sequence size

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
                kmers.push_back(std::min(kmer[0], kmer[1]));
                assert(k_ % 2 == 0 || kmer[0] != kmer[1]);
            }
        }
        
        return kmers;
    }

    std::vector<std::array<KmerId, 2>> CountAll(const DnaSeq& seq) {
        std::vector<std::array<KmerId, 2>> kmers;
        kmers.reserve(seq.Size() - k_ + 1); // Reserve space for k-mers based on sequence size

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
                assert(kmer[0] != kmer[1]);
            }
        }
        
        return kmers;
    }

    std::vector<KmerId> CountForward(const DnaSeq& seq) {
        std::vector<KmerId> kmers;
        kmers.reserve(seq.Size() - k_ + 1); // Reserve space for k-mers based on sequence size

        if (seq.Size() >= k_) {    
            KmerId kmer = 0;
            
            size_t index = 0;
            for (index = 0; index < k_-1; ++index) {
                auto c = seq[index];
                kmer = (kmer << 2 | c) & mask;           // forward k-mer
            }
            for (; index < seq.Size(); index++) {
                auto c = seq[index];
                kmer = (kmer << 2 | c) & mask;           // forward k-mer
                kmers.push_back(kmer);
            }
        }
        return kmers;
    }

    std::string ToString(KmerId kmer) {
        std::string s(k_, 'A');
        for (size_t i = 0; i < k_; ++i) {
            auto b = (kmer >> (2 * (k_ - i - 1))) & 3;
            s[i] = "ACGT"[b];
        }
        return s;
    }
    size_t k_;

protected:
    uint64_t shift1;
    uint64_t mask;
};

}