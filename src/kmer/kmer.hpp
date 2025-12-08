#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../sequence.hpp"
#include "../utils/logger.hpp"

namespace fsa {

using KmerId = unsigned long long;
struct KmerItem { KmerId kmer; int size;};

class KmerStore {
protected:

};

class KmerStoreUsingMap : public KmerStore {
public:
    KmerStoreUsingMap(const std::string& fname, size_t thread_size=1, uint16_t ct=2) : cutoff(ct) {
        Load(fname, thread_size);
    }
    
    void Load(const std::string& fname, size_t thread_size=1);

    bool Empty() const { return kmers.empty(); }
    bool Find(KmerId kid) const  { return kmers.find(kid) != kmers.end(); }
    int Count(KmerId kid) const {
        auto iter = kmers.find(kid) ;
        return iter == kmers.end() ? 0 : iter->second;
    }
    size_t Size() const { return kmers.size(); }
    size_t K() const { return k; }
    int Min() const { 
        return std::min_element(kmers.begin(), kmers.end(), 
                    [](const std::pair<KmerId, int>& a, const std::pair<KmerId, int>& b) {
                        return a.second < b.second;
                    })->second;
    }

protected:
    size_t k;
    std::unordered_map<KmerId, int> kmers;
    uint32_t cutoff {1};
};

class KmerStoreUsingBin : public KmerStore {
public:
    KmerStoreUsingBin(const std::string& fname, size_t bbit=5, size_t thread_size=1)
     : binbit(bbit), kmers(1 << (bbit*2)) {
        LOG(INFO)("KmerStoreUsingBin bbit=%d %zd %zd", kmers.size(), 1 << (bbit*2), bbit);
        Load(fname, thread_size);
    }
    
    void Load(const std::string& fname, size_t thread_size=1);

    bool Empty() const { return kmers.empty(); }
    size_t Size() const { return kmers.size(); }
    size_t K() const { return k; }
protected:
    size_t BinIndex(KmerId kid) const {// 总bit数 = k * 2
        int total_bits = k * 2;
        int middle_bits = binbit * 2;
        
        int skip_bits = (total_bits - middle_bits) / 2;
        uint64_t mask = (1ULL << middle_bits) - 1;
        return (kid >> skip_bits) & mask;
    }

    size_t k;   
    uint8_t binbit{5};
    std::vector<std::unordered_map<KmerId, int>> kmers;
};

class KmerStoreUsingVector : public KmerStore  {
public:
    KmerStoreUsingVector(size_t bbit=5) : binbit(bbit) {
        index.resize(1 << (bbit*2), {0,0});
    }

    KmerStoreUsingVector(const std::string& fname, size_t bbit=5)
     : KmerStoreUsingVector(bbit) {
        Load(fname);
    }

    void Load(const std::string& fname, size_t thread_size=1);
    bool Empty() const { return kmers.empty(); }
    size_t Count(KmerId kid) const ;
    size_t Size() const { return kmers.size(); }
    size_t K() const { return k; }
protected:
    void BuildIndex();
    uint8_t binbit{10};
    size_t k;
    std::vector<std::array<size_t,2>> index;
    std::vector<KmerItem> kmers;
};


class KmerStoreUsingSet : public KmerStore  {
public:
    KmerStoreUsingSet(size_t ks=5) : ksets(ks) {}

    KmerStoreUsingSet(const std::string& fname, size_t ks=5, size_t thread_size=1)
        : KmerStoreUsingSet(ks) {
        Load(fname, thread_size);
    }

    void Load(const std::string& fname, size_t thread_size=1);
    bool Empty() const { return kmers.empty(); }
    size_t Count(KmerId kid) const ;
    size_t Size() const { return kmers.size(); }
    size_t K() const { return k; }
protected:
    size_t k;
    std::vector<std::unordered_set<KmerId>> ksets;
    std::unordered_map<KmerId, int> kmers;
};

    


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
