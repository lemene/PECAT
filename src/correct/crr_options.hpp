#pragma once

#include "../utils/argument_parser.hpp"

#include <atomic>
namespace fsa {


struct CrrOptions {
public:
    void SetArguments(ArgumentParser &ap);
    void CheckArguments();


    int thread_size { 4 };
    int min_coverage { 4 };

    bool use_cache { false };
    bool debug { false };
    int debug_flag { 0 };  
    std::string variants;
};

} // namespace fsa
