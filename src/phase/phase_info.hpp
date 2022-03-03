#pragma once

#include <unordered_set>
#include "../sequence.hpp"

#include "../utils/string_pool.hpp"
#include "../overlap.hpp"

namespace fsa {
    

class PhaseInfoFile {
    struct Item {
        int strand : 8;
        int offset : 24;
    };


public:
    PhaseInfoFile() : PhaseInfoFile(default_string_pool_) {};
    PhaseInfoFile(StringPool &sp) : string_pool_(sp) {}
    
    void Load(const std::string &fname);
    void Load(const std::string &fname, size_t thread_size);
    bool Contain(const Overlap& ol, int threshold=1000) const;
    bool IsRemoved(int id) const { return removed_.find(id) != removed_.end() && phased_.find(id) == phased_.end(); }
    bool Contain(const std::string &name0, const std::string &name1);
    bool Contain(int id0, int id1, bool d, int off_0to1, int off_1to0, int threshold) const;
    bool Contain(int id0, int id1) const;
    const StringPool& GetStringPool() const { return string_pool_; }
    std::unordered_set<int> Get(int id) const;
protected:
    StringPool default_string_pool_;
    StringPool& string_pool_;

    std::unordered_map<int, std::unordered_map<int, std::vector<Item>>> phased_;
    std::unordered_set<int> removed_;
};


} // namespace fsa
