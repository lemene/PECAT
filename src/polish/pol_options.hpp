#pragma once

#include "../utils/argument_parser.hpp"

#include "overlap.hpp"

#include <atomic>
namespace fsa {


struct PolOptions {
public:
    
    void SetArguments(ArgumentParser &ap);
    void CheckArguments();

    std::string OutputPath(const std::string &fname) { return output_directory+"/"+fname; }

    std::string output_directory {"."};

    std::string filter0_opts_ {"l=2000:al=2000:alr=0.50"};
    std::string filter1_opts_ {"l=2000:al=3000:alr=0.50:aal=6000:oh=2000:ohr=0.2"};

    Overlap::Filter filter0_;
    Overlap::Filter filter1_; 

    std::string aligner_ { "diff" };
    std::string score_ { "weight" };

    std::string read_name_ {""};
    std::string read_name_fname_ { "" };
    
    
    int min_coverage_ { 4 };
    double min_identity_ { 60 };
    double min_local_identity_ { 50 };
    int local_window_size_ { 1000 };
    bool check_local_identity_ { false };

    int coverage_ { 50 };
    int window_size_ { 50000 };
    int overlap_size_ { 500 };

    std::string ctg_fname_;
    std::string overlap_fname_;
    std::string rread_fname_;
    std::string cread_fname_;

    
    std::string graph_fname_ {""};
    std::string infos_fname_ {""}; 
    std::string rd_2_ref_fname_ { "" };

    int thread_size { 4 };
    int min_coverage { 4 };
    bool skip_branch_check { false };
    bool use_cache { false };
    bool debug { false };
    std::string variants;
    double secondary_to_primary_ratio { 0.80 };
};

} // namespace fsa
