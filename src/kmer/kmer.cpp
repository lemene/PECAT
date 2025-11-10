#include "kmer.hpp"

#include <algorithm>

#include "../utils/logger.hpp"
#include "../sequence.hpp"
#include "../file_io.hpp"
#include "../utility.hpp"
namespace fsa {

void KmerSet::BuildIndex() {
    index.assign(1024, {kmers.size(),0});

    for (size_t i = 0; i<kmers.size(); ++i) {
        size_t idx = kmers[i].kmer >> (k*2 - 10);

        if (index[idx][0] > i) {
            index[idx][0] = i;
        }

        if (index[idx][1] < i+1) {
            index[idx][1] = i+1;
        }
    }

}

bool KmerSet1::Find(KmerId kid) const {
    auto se = index[kid >> (k*2 - 10)];
    size_t s = se[0]; 
    size_t e = se[1];
    //printf("s e %zd %zd\n", s, e);

    while (s < e) {
        size_t m = (s+e) / 2;
        if (kmers[m].kmer == kid) {
            return true;
        } else if (kmers[m].kmer < kid) {
            s = m+1;
        } else {
            e = m;
        }
    }
    return false;
}


auto LoadKmers0(const std::string &fname) -> KmerSet0 {
    KmerSet0 kmers;
    std::mutex mutex_gen;
    std::mutex mutex_comb;

    kmers.k = GetKmerLength(fname);

    std::atomic<size_t> s { 0 };

    const size_t block_size = 1000;
    GzFileReader reader(fname);
    auto generate_func = [&mutex_gen, &reader](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        return reader.GetLines(lines);
    };

    auto combine_func = [&mutex_comb, &kmers](std::unordered_map<KmerId, int>& ks) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        kmers.kmers.insert(ks.begin(), ks.end());
        ks.clear();
    };

    auto work_func = [block_size, generate_func, combine_func, &reader, &s](size_t id) {
        std::vector<std::string> lines(block_size);
        std::unordered_map<KmerId, int> ks;

        size_t sz = generate_func(lines);
        while (sz > 0) {
            
            for (size_t i=0; i<sz; ++i) {
                auto items = SplitStringBySpace(lines[i]);
                ks[KmerStringToId(items[0])] = std::stoi(items[1]);
                s.fetch_add(std::stoi(items[1]));
            }
            
            combine_func(ks);
            sz = generate_func(lines);
        }
    };
    
    int thread_size_ = 10;
    MultiThreadRun(std::min<size_t>(thread_size_, 3), work_func);

    LOG(INFO)("Load %zd kmers, %zd occu (k=%zd) from %s", kmers.kmers.size(), s.fetch_add(0), kmers.k, fname.c_str());
    return kmers;
}

auto LoadKmers1(const std::string &fname) -> KmerSet {
    KmerSet kmers;

    kmers.k = GetKmerLength(fname);

    GzFileReader reader(fname);
    std::string line;
    while (reader.GetLine(line)) {
        auto items = SplitStringBySpace(line);
        kmers.kmers.push_back({KmerStringToId(items[0]),std::stoi(items[1])});
    }

    LOG(INFO)("Load %zd kmers(k=%zd) from %s", kmers.kmers.size(), kmers.k, fname.c_str());
    std::sort(kmers.kmers.begin(), kmers.kmers.end(), [](const KmerItem& a, const KmerItem &b) { return a.kmer < b.kmer; });
    LOG(INFO)("Sort %zd kmers(k=%zd) from %s", kmers.kmers.size(), kmers.k, fname.c_str());
    kmers.BuildIndex();
    return kmers;
}


auto KmerStringToId(const std::string& str) -> KmerId{
    static DnaSerialTable table;
    KmerId id = 0;
    for (auto c : str) {
        id = (id << 2) + table[c];
    }
    return id;    
}

std::string KmerId2String(KmerId kmer, size_t k) {
    std::string s(k, 'A');
    for (size_t i = 0; i < k; ++i) {
        auto b = (kmer >> (2 * (k - i - 1))) & 3;
        s[i] = "ACGT"[b];
    }
    return s;
}


size_t GetKmerLength(const std::string &fname) {
    GzFileReader reader(fname);
    if (reader.Valid()) {
        auto line = reader.GetNoEmptyLine();
        if (!line.empty()) {
            auto items = SplitStringBySpace(line);
            return items[0].size();

        }

    }
    return 0;
}
} // namespace fsa