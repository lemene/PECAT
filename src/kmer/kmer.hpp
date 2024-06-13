#pragma once

#include <unordered_map>
#include <vector>

namespace fsa {

using KmerId = unsigned long long;
struct KmerItem { KmerId kmer; int size;};
struct KmerSet0 {
    bool Empty() const { return kmers.empty(); }
    bool Find(KmerId kid) const  { return kmers.find(kid) != kmers.end(); }
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

} // namespace fsa
