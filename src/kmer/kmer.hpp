#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "../sequence.hpp"
namespace fsa {

using KmerId = unsigned long long;
struct KmerItem { KmerId kmer; int size;};
struct KmerSet0 {
    bool Empty() const { return kmers.empty(); }
    bool Find(KmerId kid) const  { return kmers.find(kid) != kmers.end(); }
    int Count(KmerId kid) const {
        auto iter = kmers.find(kid) ;
        return iter == kmers.end() ? 0 : iter->second;
    }
    size_t Size() const { return kmers.size(); }
    size_t k;
    std::unordered_map<KmerId, int> kmers;
};
struct KmerSet1 {
    bool Empty() const { return kmers.empty(); }
    bool Find(KmerId kid) const;
    size_t Size() const { return kmers.size(); }
    void BuildIndex();
    size_t k;
    std::vector<std::array<size_t,2>> index;
    std::vector<KmerItem> kmers;
};
using KmerSet = KmerSet1;


class KmerCount { 
public:
    KmerCount(size_t k) : k_(k) {
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

KmerSet0 LoadKmers0(const std::string &fname);
KmerSet1 LoadKmers1(const std::string &fname);

KmerId KmerStringToId(const std::string &str);
std::string KmerId2String(KmerId id, size_t k);
size_t GetKmerLength(const std::string &fname);
} // namespace fsa
