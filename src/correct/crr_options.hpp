#pragma once

#include "../utils/argument_parser.hpp"

#include "overlap.hpp"

#include <atomic>
namespace fsa {


struct CrrOptions {
public:
    // the options for selecting candiate overlaps
    struct CandidateOptions {
        CandidateOptions(const std::string &str) { From(str); }
        void From(const std::string &str);
        std::string ToString() const;

        bool IsEndCondition(const std::vector<int> &cov, size_t number, size_t fail) const {
            return (int)fail >= failures || IsEnough(cov);
        }
        bool IsEnough(const std::vector<int> &cov) const ;
        bool IsEnough(size_t number) const { return false; }

        double percent { 0.95 };             // p Percentage of filled matrix
        double overhang_weight   { 0.0 };                // w overhang的比重
        int failures { 10 };               // f 连续失败次数
        int max_number { 200 };             // 
        int coverage { 80 };                    // 需要多少层数据
    };
    
    void SetArguments(ArgumentParser &ap);
    void CheckArguments();

    std::string OutputPath(const std::string &fname) { return output_directory+"/"+fname; }

    std::string output_directory {"."};

    std::string filter0_opts_ {"l=2000:al=2000:alr=0.50"};
    std::string filter1_opts_ {"l=2000:al=3000:alr=0.50:aal=6000:oh=2000:ohr=0.2"};
    std::string cands_opts_str_ { "c=80:f=10:p=0.95:ohwt=0.1"};


    Overlap::Filter filter0_;
    Overlap::Filter filter1_; 
    CandidateOptions cands_opts_ { cands_opts_str_ };

    std::string aligner_ { "diff" };
    std::string score_ { "weight" };

    std::string read_name_ {""};
    std::string read_name_fname_ { "" };
    
    double min_identity_ { 60 };
    double min_local_identity_ { 50 };
    bool check_local_identity_ { false };

    std::string overlap_fname_;
    std::string rread_fname_;
    std::string cread_fname_;
    std::string graph_fname_ {""};
    std::string infos_fname_ {""}; 

    int thread_size { 4 };
    int min_coverage { 4 };
    bool skip_branch_check { false };
    bool use_cache { false };
    bool debug { false };
    std::string variants;
};

} // namespace fsa
