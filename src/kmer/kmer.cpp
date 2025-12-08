#include "kmer.hpp"

#include <algorithm>

#include "../utils/logger.hpp"
#include "../sequence.hpp"
#include "../file_io.hpp"
#include "../utility.hpp"
namespace fsa {



void KmerStoreUsingMap::Load(const std::string &fname, size_t thread_size) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;

    k = GetKmerLength(fname);
    auto count = CountLinesInFile(fname, std::min<size_t>(thread_size, 20));
    LOG(INFO)("Loading KmerStoreUsingMap from %s: size=%zd, k=%d", fname.c_str(), count, k);

    kmers.reserve(count*1.5);
    std::shared_ptr<Reader> in = fname == "-" ?
        std::shared_ptr<Reader>(new StdioReader()) :
        std::shared_ptr<Reader>(new GzFileReader(fname));

    auto combine_func = [&mutex_comb, this](std::vector<std::pair<KmerId, int>>& ks) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        kmers.insert(ks.begin(), ks.end());
        ks.clear();
    };

    auto work_func = [&mutex_gen, combine_func, in, this](size_t _) {
        LineInBlock line_in_block(*in, 10000000, &mutex_gen);
        const size_t block_size = 10000000;

        std::vector<std::pair<KmerId, int>> ks;

        std::string line;
        for (auto valid = line_in_block.GetLine(line); valid; valid = line_in_block.GetLine(line)) {
            auto items = SplitStringBySpace(line);
            int c =  std::stoi(items[1]);
            if (c >= cutoff) {
                ks.push_back({KmerStringToId(items[0]), c});
            }

            if (ks.size() >= block_size) {
                combine_func(ks);
            }
        }
        LOG(INFO)("Thread %zd loaded %zd kmers", _, ks.size());
        combine_func(ks);
        LOG(INFO)("combine_func %zd loaded %zd kmers", _, ks.size());
    };
    
    MultiThreadRun(thread_size, work_func);
}


void KmerStoreUsingBin::Load(const std::string &fname, size_t thread_size) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;

    k = GetKmerLength(fname);
    
    

    GzFileReader reader(fname);

    auto combine_func = [&mutex_comb, this](std::vector<std::pair<KmerId, int>>& ks) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        for (const auto& p : ks) {
            size_t bidx = BinIndex(p.first);
            kmers[bidx][p.first] = p.second;
        }
        ks.clear();
    };

    auto work_func = [&mutex_gen, combine_func, &reader](size_t _) {
        LineInBlock line_in_block(reader, 10000000, &mutex_gen);
        const size_t block_size = 10000;

        std::vector<std::pair<KmerId, int>> ks;

        std::string line;
        for (auto valid = line_in_block.GetLine(line); valid; valid = line_in_block.GetLine(line)) {
            auto items = SplitStringBySpace(line);
            ks.push_back({KmerStringToId(items[0]), std::stoi(items[1])});

            // if (ks.size() >= block_size) {
            //     combine_func(ks);
            // }
        }
        LOG(INFO)("Thread %zd loaded %zd kmers", _, ks.size());
        combine_func(ks);
        LOG(INFO)("combine_func %zd loaded %zd kmers", _, ks.size());
    };
    
    MultiThreadRun(thread_size, work_func);
}

void KmerStoreUsingVector::Load(const std::string& fname, size_t thread_size) {
    k = GetKmerLength(fname);

    GzFileReader reader(fname);
    std::string line;
    
    while (reader.GetLine(line)) {
        auto items = SplitStringBySpace(line);
        kmers.push_back({KmerStringToId(items[0]),std::stoi(items[1])});
    }
    LOG(INFO)("Load %zd kmers(k=%zd) from %s", kmers.size(), k, fname.c_str());
    std::sort(kmers.begin(), kmers.end(), [](const KmerItem& a, const KmerItem &b) { return a.kmer < b.kmer; });
    LOG(INFO)("Sort %zd kmers(k=%zd) from %s", kmers.size(), k, fname.c_str());

    BuildIndex();
}

void KmerStoreUsingVector::BuildIndex() {
    index.assign(1024, {kmers.size(), 0});

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



size_t KmerStoreUsingVector::Count(KmerId kid) const {
    auto se = index[kid >> (k*2 - 10)];
    size_t s = se[0]; 
    size_t e = se[1];

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

void KmerStoreUsingSet::Load(const std::string &fname, size_t thread_size) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;

    k = GetKmerLength(fname);
    GzFileReader reader(fname);

    auto combine_func = [&mutex_comb, this](std::vector<std::pair<KmerId, int>>& ks) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        for (const auto& p : ks) {
            if (p.second <= ksets.size()) {
                ksets[p.second - 1].insert(p.first);
            } else {
                kmers[p.first] = p.second;
            }
        }
        ks.clear();
    };

    auto work_func = [&mutex_gen, combine_func, &reader](size_t _) {
        LineInBlock line_in_block(reader, 10000000, &mutex_gen);
        const size_t block_size = 1000000;

        std::vector<std::pair<KmerId, int>> ks;

        std::string line;
        for (auto valid = line_in_block.GetLine(line); valid; valid = line_in_block.GetLine(line)) {
            auto items = SplitStringBySpace(line);
            ks.push_back({KmerStringToId(items[0]), std::stoi(items[1])});

            // if (ks.size() >= block_size) {
            //     combine_func(ks);
            // }
        }
        LOG(INFO)("Thread %zd loaded %zd kmers", _, ks.size());
        combine_func(ks);
    };
    
    MultiThreadRun(thread_size, work_func);
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