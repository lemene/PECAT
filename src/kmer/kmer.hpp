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
    int Min() const { 
        return std::min_element(kmers.begin(), kmers.end(), 
                    [](const std::pair<KmerId, int>& a, const std::pair<KmerId, int>& b) {
                        return a.second < b.second;
                    })->second;
    }
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



KmerSet0 LoadKmers0(const std::string &fname);
KmerSet1 LoadKmers1(const std::string &fname);

KmerId KmerStringToId(const std::string &str);
std::string KmerId2String(KmerId id, size_t k);
size_t GetKmerLength(const std::string &fname);

static inline uint64_t hash64(uint64_t key)
{
    key = (~key + (key << 21));
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key = key ^ key >> 28;
    key = (key + (key << 31));
    return key;
}


static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

} // namespace fsa
