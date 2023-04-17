#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <mutex>
#include <algorithm>
#include "./utils/logger.hpp"
#include "sequence.hpp"
#include "seq_io.hpp"
#include "utils/string_pool.hpp"

namespace fsa {

class BaseStore {
public:
    void Add(const char* str, size_t len);
    size_t Size() const { return size_; }
protected:
    uint8_t& Locate(size_t i)  { 
        size_t ii = i/4;
        return data_[ii/block_size_][ii%block_size_];
    }
protected:
    std::vector<std::vector<uint8_t>> data_;
    const int block_size_ = 1000000000;
    size_t size_;
};

class SequenceStore {
public:
    SequenceStore(StringPool &sp) : string_pool_(sp) {}
    SequenceStore() : SequenceStore(default_string_pool_) {}


    struct Sequence {
        size_t offset;
        size_t len;
    };

    Seq::Id QueryIdByName(const std::string &name) const;
    const std::string& QueryNameById(Seq::Id id) const;

    StringPool& GetStringPool() { return string_pool_; }
    const StringPool& GetStringPool() const { return string_pool_; }

    const Sequence& Get(size_t i) {
        assert(i >= 0 && i < items_.size());
        return items_[i];
    } 
    
    size_t Size() const { return items_.size(); }

    void Load(const std::string &fname, const std::string &type="");
    void LoadFasta(const std::string &fname);
    void LoadFastq(const std::string &fname);
    void LoadFofn(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadTxt(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) { LoadFofn(fname, all, seqids); }

    // @return start and end poistion of head, subhead, seq and quality. 
    std::vector<std::array<size_t, 2>> GetFastqSegment(const std::vector<char> &block, size_t bsize, size_t start=0);
    
    static std::string DetectFileType(const std::string &fname);

protected:
    StringPool default_string_pool_;
    StringPool &string_pool_;   
    std::vector<Sequence> items_;

    std::unordered_map<Seq::Id, size_t> id_2_index_;
    BaseStore data_;
};

///
size_t FindLastFastqInBlock(const std::vector<char> &block, size_t bsize);

} // namespace fsa {

