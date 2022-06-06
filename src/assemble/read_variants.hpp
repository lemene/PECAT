#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "utils/string_pool.hpp"

namespace fsa {

class Overlap;


class ReadVariants {
public:
    struct Variants {
        StringPool::ID contig;
        char d;
        int offset;     
        std::unordered_map<int, int> vars;
    };

public:
    ReadVariants(StringPool &sp) : string_pool_(sp) {}
    void Load(const std::string &fname);

    std::array<int, 2> GetClosestSnps(const Overlap &ol) const;
    std::array<int, 2> GetSnps(const Variants &a, const Variants &b) const;
    std::array<int, 2> Test(const Overlap &ol, bool bad=false) const;
    std::array<int, 2> Test(int a, int b) const ;
    template<typename T>
    std::array<int, 2> Test(const T& as, int b) const  {
        std::array<int, 2> result = {0, 0};
        for (auto &a : as) {
            auto r = Test(a, b);
            result[0] += r[0];
            result[1] += r[1];
        }
        return result;
    }

    bool IsCompatible(const Variants &a, const Variants &b, const Overlap &ol) const;
    const std::vector<Variants>* GetVariants(int a) {
        auto iter = reads.find(a); 
        return iter != reads.end() ? &iter->second : (std::vector<Variants>*)nullptr;
    }

    std::unordered_map<StringPool::ID, std::vector<Variants>> reads;
    StringPool &string_pool_;
};

}