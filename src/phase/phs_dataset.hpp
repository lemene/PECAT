#pragma once

#include "../overlap_store.hpp"
#include "../read_store.hpp"

namespace fsa {

class PhsOptions;

class PhsDataset {
public:
    PhsDataset(PhsOptions &opts) : opts_(opts) {}

    void Load();
    std::vector<Seq::Id> GetSortedContigs() const;
    const std::string& QueryNameById(Seq::Id i) const { return rd_store_.QueryNameById(i); }
    size_t GetOverlapSize(Seq::Id i) const {
        auto iter = read_ids2_.find(i);
        return iter != read_ids2_.end() ? iter->second : 0;
    }

    
    std::unordered_set<int> QueryGroup(int id);
    bool QueryAva(int r0, int r1) const {
        auto i0 = ava_groups_.find(r0);
        if (i0 != ava_groups_.end()) {
            return i0->second.find(r1) != i0->second.end();
        } else {
            return false;
        }
    }
    bool HasAva() const { return ol_ava_.Size() > 0; }
protected:
    
    // Load twice to avoid loading unwanted reads
    void PreloadReads();
    void LoadReads();
    void LoadOverlaps(const std::string &fname);
    void LoadAva(const std::string &fname);
    
    void LoadSnpFromVcf(const std::string &fname, StringPool& string_pool);
public:
    PhsOptions &opts_;

    ReadStore rd_store_;
    OverlapStore ol_store_ { rd_store_.GetStringPool() };

    OverlapStore ol_ava_ { rd_store_.GetStringPool() };
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> ava_groups_;
    
    std::unordered_set<std::string> contig_names_;
    std::unordered_set<Seq::Id> contig_ids_;
    std::unordered_map<Seq::Id, size_t> read_ids2_;

    std::unordered_map<Seq::Id, std::unordered_map<size_t, std::array<uint8_t, 2>>> snps_;
};

} // namespace fsa
