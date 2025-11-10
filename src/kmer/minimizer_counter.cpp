#include "minimizer_counter.hpp"
#include "../utils/logger.hpp"

namespace fsa {

template<typename T>
struct CircleCache {
public:
    CircleCache(size_t w) : data(w) {}
    size_t Capacity() const { return data.size() ; }
    void Push(const T& d) { 
        data[top] = d;
        Increase(top);
    }
    size_t Size() const { return data.size() - bottom + top; }
    void Pop() {
        Increase(bottom);
    }
    T& operator [](size_t i) { return data[(i + bottom) % data.size()]; }
    const T& operator [](size_t i) const { return data[(i + bottom) % data.size()]; }
    void Increase(size_t &i) {
        i ++;
        if (i >= data.size()) {
            i = 0;
        }
    }
    void Debug_Print() const {
        for (size_t i = 0; i < data.size(); ++i) {
            LOG(INFO)("cache %zd %zd %016llX", i, data[i].pos, data[i].hash);
        }
    }

    size_t top { 0 };
    size_t bottom { 0 }; 
    std::vector<T> data;
};


std::vector<Minimizer> MinimizerCounter::Count(const DnaSeq& seq, const RankedKmers& rkmers) {
    auto kmers = kmc_.CountAll(seq);
    //LOG(INFO)("Counted %zd kmers from sequence of length %zd", kmers.size(), seq.Size());
    // for (auto &k : kmers) {
    //     LOG(INFO)("kmer %s %s %s", kmc_.ToString(k[0]).c_str(), kmc_.ToString(k[1]).c_str(), kmc_.ToString(std::min(k[0], k[1])).c_str());
    // }
    uint64_t mask = (1ULL<<2*k_) - 1;
    std::vector<Minimizer> minimizers;
    minimizers.reserve(kmers.size() / w_); // Reserve space for minimizers based


    auto to_minimizer = [mask, &rkmers, this](const std::array<KmerId, 2>& kmer, size_t pos, size_t rid) -> Minimizer {
        Minimizer m;
        
        if (kmer[0] < kmer[1]) {
            m.dir = 0; // Forward direction
            m.hash = hash64(kmer[0], mask); 
            m.rank = rkmers.GetRank(kmer[0]);
            m.kmer = kmer[0];
        } else {
            if (kmer[0] == kmer[1]) {
                LOG(INFO)("kmer %s %s", kmc_.ToString(kmer[0]).c_str(), kmc_.ToString(kmer[0]).c_str());
            }
            
            assert(kmer[0] > kmer[1]);
            m.dir = 1; // Reverse direction 
            m.hash = hash64(kmer[1], mask); 
            m.rank = rkmers.GetRank(kmer[1]);
            m.kmer = kmer[1];
        }
        if (kmer[0] == 10100601080 || kmer[1] == 10100601080) {
            LOG(INFO)("kmer %s %s %zd %zd", kmc_.ToString(kmer[0]).c_str(), kmc_.ToString(kmer[1]).c_str(), m.rank, m.hash);
        }
        m.pos = pos; // Position in the sequence
        m.rid = rid; // Read ID, can be set based on context
        return m;
    };

    CircleCache<Minimizer> cache(w_);
    Minimizer min;
    for (size_t i = 0; i < w_; ++i) {
        cache.Push(to_minimizer(kmers[i], i, 0));
    }

    auto find_min_in_cache = [](const CircleCache<Minimizer>& cache, std::vector<Minimizer>& minimizers, Minimizer& min) {
        auto mm = std::min_element(cache.data.begin(), cache.data.end(), [](const Minimizer& a, const Minimizer& b) {
            return a.rank < b.rank || (a.rank == b.rank && a.hash < b.hash);
        });
    
        for (size_t i = 0; i != cache.Size(); ++i) {
            if (cache[i].hash == mm->hash) {
                minimizers.push_back(cache[i]);
                min = cache[i];
            }
        }
    };
    
    

    find_min_in_cache(cache, minimizers, min);

    for (size_t i = w_; i < kmers.size(); ++i) {
        auto mm = to_minimizer(kmers[i], i, 0);
        cache.Pop();
        cache.Push(mm);

        if (mm.rank < min.rank || (mm.rank == min.rank &&  mm.hash <= min.hash)) {
            min = mm;
            minimizers.push_back(mm);
        } else if (i >= min.pos + cache.Size()) {
            find_min_in_cache(cache, minimizers, min);
        }
    }

    return minimizers;
}


std::vector<Minimizer> MinimizerCounter::Count(const DnaSeq& seq) {
    return Count(seq, RankedKmers(std::vector<std::string>(), std::vector<double>(), 1));
}
}


