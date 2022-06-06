#pragma once

#include "../utils/argument_parser.hpp"
#include "../overlap.hpp"

namespace fsa {

struct AsmOptions {
    void SetArguments(ArgumentParser &ap);
    void CheckArguments() {
        filter0.From(filter0_opts);
        filter0_opts = filter0.ToString();
    }

    std::string OutputPath(const std::string &fname) const { return output_directory +"/"+fname; }

    static bool ParamToGenomeSize(const std::string& str, long long *v);

    void UpdateByGenomeSize(size_t gsize);

    std::string filter0_opts {"l=5000:al=2500:i=10:oh=100"};               //!< overlap filtering options when loading 
    Overlap::Filter filter0 { filter0_opts };
    int thread_size { 4 };           //!< 
    std::string output_directory {"."};
    int min_coverage { 2 }; 
    double min_identity { 0.95 };
    int max_spur_length { 1000000 };
    int max_spur_nodesize { 20 } ; 
    int max_trivial_length { -1 };

    double max_offset_rate {0.05};

    std::string read_file {""};
    std::string debug_name {""};

    std::string ifname;                        
    long long genome_size {0};
    int coverage { 60 };
    std::string coverage_opts { "2|1000|500" };

    std::string reduction0 {"transitive:fuzz=500|quality|spur|best|spur|phase|unreliable|spur|bridge|"};
    std::string reduction1 {"spur|duplicate|cross|bridge|bubble:s=1|bridge|bubble:s=0|spur|loop|semi"};

    std::string reducer0 {""};
    std::string reducer1 {""};
    
    int max_unreliable_length { 10000 };
    double max_unreliable_rate { 0.66 };

    int min_contig_length { 500 };
    // int max_bubble_identity { 96 };
    // int max_bubble_coverage { 97 };
    int max_bubble_length { 50000000 };
    std::string select_branch { "no" }; 
    std::string phased;
    std::string variants;

    int diploid_count { 3 };
    double diploid_rate { 0.30 };

    int dump { 0 };
    bool skip_purge { false };
};

}