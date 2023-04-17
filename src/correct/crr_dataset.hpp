#pragma once

#include <atomic>

#include "test_future.hpp"

namespace fsa {

class CrrOptions;

struct CrrDataset {
public:
    CrrDataset(const CrrOptions &opt) : opts_(opt) {}
    const CrrOptions& opts_;

    void Load(const class StringPool &sp);    
    std::shared_ptr<SnpFile> variants;
};

} // namespace fsa
