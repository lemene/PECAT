#pragma once

#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "corrector.hpp"
#include "test_future.hpp"

namespace fsa {

class CrrOptions;

struct CrrDataset {
public:
    CrrDataset(const CrrOptions &opt) : opts_(opt) {}
    const CrrOptions& opts_;

    void Load(const class StringPool &sp);    
    std::shared_ptr<SnpFile> variants;

    
    StringPool string_pool_;
    ReadStore read_store_ {string_pool_};
    OverlapStore ol_store_{string_pool_ };
    
    OverlapGrouper grouper_ { ol_store_ };
    
    std::vector<Seq::Id> read_ids_;
    std::vector<Seq::Id> grouped_ids_;
    std::vector<size_t> group_ticks;
    size_t  group_size  = 1000;
};

} // namespace fsa
