#pragma once

#include <string>

#include "utils/string_pool.hpp"

namespace fsa {
    
class VariantStore {
public:
    VariantStore(class StringPool &sp) : string_pool_(sp) {}

    void LoadFromVcf(const std::string &fname);

protected:
    class StringPool &string_pool_;
};


}