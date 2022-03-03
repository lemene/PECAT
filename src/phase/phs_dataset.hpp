#pragma once

#include "../overlap_store.hpp"
#include "../read_store.hpp"

namespace fsa {

class PhsOptions;

class PhsDataset {
public:
    PhsDataset(const PhsOptions &opts) : opts_(opts) {}

    void Load();
    std::vector<Seq::Id> GetSortedContigs() const;
    const std::string& QueryNameById(Seq::Id i) const { return rd_store_.QueryNameById(i); }
    size_t GetOverlapSize(Seq::Id i) const {
        auto iter = read_ids2_.find(i);
        return iter != read_ids2_.end() ? iter->second : 0;
    }
protected:
    
    // Load twice to avoid loading unwanted reads
    void PreloadReads();
    void LoadReads();
    void LoadOverlaps(const std::string &fname);
    
public:
    const PhsOptions &opts_;

    ReadStore rd_store_;
    OverlapStore ol_store_ { rd_store_.GetStringPool() };
    
    std::unordered_set<std::string> contig_names_;
    std::unordered_set<Seq::Id> contig_ids_;
    std::unordered_map<Seq::Id, size_t> read_ids2_;

};

} // namespace fsa
