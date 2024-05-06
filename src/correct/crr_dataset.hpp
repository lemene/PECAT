#pragma once

#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "corrector.hpp"

namespace fsa {

class CrrOptions;

struct CrrDataset {
public:
    CrrDataset(CrrOptions &opt) : opts_(opt) {}

    void Load();
    const std::string& QueryStringById(Seq::Id id) const { return string_pool_.QueryStringById(id); }    

protected:
    void LoadReadIds();
    void LoadOverlaps();
    void LoadReads();
    void GroupReadIds();
    void EstimateParameters();
public:
    CrrOptions& opts_;
    
    StringPool string_pool_;
    ReadStore read_store_ {string_pool_};
    OverlapStore ol_store_{string_pool_ };
    
    OverlapGrouper grouper_ { ol_store_ };
    
    std::vector<Seq::Id> read_ids_;
    //std::vector<Seq::Id> grouped_ids_;
    //std::vector<size_t> group_ticks;o
    std::vector<std::vector<Seq::Id>> clu_ids_;
};

} // namespace fsa
