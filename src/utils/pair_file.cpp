#include "pair_file.hpp"

#include <cassert>
#include <array>

#include "../file_io.hpp"
#include "logger.hpp"
#include "../utility.hpp"

namespace fsa {


void PairFile::Load(const std::string &fname, size_t thread_size) {
    std::mutex mutex_combine;
    std::mutex mutex_generate;

    GzFileReader in(fname);
    std::array<long long, 2> done = {0,0};

    int block_size = 2000;
    auto generate = [&mutex_generate, &in, &done](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_generate);
        auto size = in.GetLines(lines);
        done[0] += size;

        return size;
    };

    auto pair_value = [](int a, int b) {
        return a < b ? ((long long)a << 32) + b : ((long long)b << 32) + a;
    };

    auto combine = [this, &mutex_combine, pair_value](const std::vector<std::vector<StringPool::ID>> &pair0, 
                                                      const std::vector<std::vector<StringPool::ID>> &pair1, size_t sz) {
        std::lock_guard<std::mutex> lock(mutex_combine);
        for (size_t i=0; i<sz; ++i) {
            for (auto p : pair0[i]) {
                auto iter = pairs_.find(p);
                if (iter != pairs_.end()) {
                    iter->second.insert(pair1[i].begin(), pair1[i].end());
                } else {
                    pairs_[p] = std::unordered_set<StringPool::ID>(pair1[i].begin(), pair1[i].end());
                }
            }
        }
    };

    auto work = [&](size_t threadid) {
        StringPool subpool;
        std::vector<std::string> lines(block_size);
        std::vector<std::vector<StringPool::ID>> pair0(block_size);
        std::vector<std::vector<StringPool::ID>> pair1(block_size);
        
        size_t size = generate(lines);
        while (size > 0) {
            for (size_t i=0; i<size; ++i) {
                auto sets = SplitStringByChar(lines[i], ':');
                assert(sets.size() == 2);
                auto set0 = SplitStringByChar(sets[0], ',');
                auto set1 = SplitStringByChar(sets[1], ',');
                
                pair0[i].assign(set0.size(), 0);
                for (size_t j=0; j<set0.size(); ++j) {
                    pair0[i][j] = subpool.GetIdByStringUnsafe(set0[j]);
                }
                pair1[i].assign(set1.size(), 0);
                for (size_t j=0; j<set1.size(); ++j) {
                    pair1[i][j] = subpool.GetIdByStringUnsafe(set1[j]);
                }
            }

            auto mapper = string_pool_.Merge(subpool);
            subpool.Clear();
            for (size_t i=0; i<size; ++i) {
                for (size_t j=0; j<pair0[i].size(); ++j) {
                    pair0[i][j] = mapper[pair0[i][j]];
                }
            }
            for (size_t i=0; i<size; ++i) {
                for (size_t j=0; j<pair1[i].size(); ++j) {
                    pair1[i][j] = mapper[pair1[i][j]];
                }
            }

            combine(pair0, pair1, size);
            size = generate(lines);
        }
    };


    if (in.Valid()) {
        MultiThreadRun(thread_size, work);
    } else {
        LOG(ERROR)("Failed to load file: %s", fname.c_str());
    }
    LOG(INFO)("Load ignored pairs: %zd", pairs_.size());

}

bool PairFile::Contain(const std::string &n1, const std::string &n2) {
    auto id1 = string_pool_.QueryIdByString(n1);
    auto id2 = string_pool_.QueryIdByString(n2);

    if (string_pool_.Valid(id1) && string_pool_.Valid(id2)) {
        auto iter1 = pairs_.find(id1);
        auto iter2 = pairs_.find(id2);
        return (iter1 != pairs_.end() && iter1->second.find(id2) != iter1->second.end()) ||
               (iter2 != pairs_.end() && iter2->second.find(id1) != iter2->second.end()) ;
    } else {
        return false;
    }
}

}