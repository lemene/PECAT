#pragma once

#include <unordered_map>
#include <unordered_set>

#include "../utils/string_pool.hpp"
namespace fsa {


class SnpFile {
public:
    SnpFile(const StringPool& sp) : sp_(sp) {}

    void Load(const std::string &fname);

    const std::unordered_set<int>& QuerySnps(StringPool::ID id) const;

protected:
    std::unordered_set<int>& GetSnpSet(StringPool::ID id); 
    std::unordered_map<StringPool::ID, std::unordered_set<int>> snps_;
    const StringPool& sp_;
    std::unordered_set<int> empty_;
};

} // namespace fsa {


