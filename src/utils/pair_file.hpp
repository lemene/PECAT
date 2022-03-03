#pragma once

#include <unordered_set>

#include "string_pool.hpp"

namespace fsa {



class PairFile {
public:
    PairFile() : PairFile(default_string_pool_) {};
    PairFile(StringPool &sp) : string_pool_(sp) {}
    
    void Load(const std::string &fname, size_t thread_size);
    bool Contain(const std::string& n1, const std::string &n2);
    const std::unordered_set<int> &Get(int i) const {
        auto iter = pairs_.find(i);
        return iter != pairs_.end() ? iter->second : empty_;
    }
protected:
    long long int PairValue(int a, int b) { return a < b ? ((long long)a << 32) + b : ((long long)b << 32) + a; }

protected:
    StringPool default_string_pool_;
    StringPool& string_pool_;

    std::unordered_map<int, std::unordered_set<int>> pairs_;
    std::unordered_set<int> empty_;
};

}