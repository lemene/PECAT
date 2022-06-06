#pragma once

#include "../utils/argument_parser.hpp"

#include "utility.hpp"
#include <atomic>
namespace fsa {

#define VARIANT_TYPE_M 1
#define VARIANT_TYPE_D 2
#define VARIANT_TYPE_I 4

struct PhsOptions {
public:

    struct PhaserOptions : public StringOptions {
        PhaserOptions() { }
        PhaserOptions(const std::string &str) { From(str); }
        void From(const std::string &str);
        std::string ToString() const;

        // group reads
        int     min_count { 2 }; // min_consist_count / 2
        int     min_coverage_count { 6 };
        double  min_coverage_rate { 0.2 };

        // kmean distribute reads
        int     min_similar_count { 2 };
        double  min_similar_rate { 0.5 };
        double  min_similar_coverage { 0.3 };
        double  min_similar_diff { 0.01 };

        // kmean sse

        int     number_of_iteration { 3 };
        int     min_support_count { 6 };
        double  min_support_rate { 0.66 };

        int min_link_support_count { 2 };
        double min_link_support_rate { 0.33 };
        int min_link_valid_count { 2 };
        double min_link_valid_rate { 0.33 };
        
        // consist    
        double min_consist_rate = 0.7;
        int min_consist_count = 4;
        double min_consist_coverage = { 0.66 };

        // inconsist
        int min_inconsist_count = 4;
        double min_inconsist_rate = 0.05;
        std::string debug;
    };

    void SetArguments(ArgumentParser &ap);
    void CheckArguments();

    std::string OutputPath(const std::string &fname) const { return output_directory_ +"/"+fname; }

    std::string Consistent() const { return OutputPath("consistent"); }
    std::string Inconsistent() const { return OutputPath("inconsistent"); }
    std::string Varaints() const { return OutputPath("variants"); }
    std::string Readinfos() const { return OutputPath("readinfos"); }

    std::string ol_fname_;
    std::string rd_fname_;
    std::string ctg_fname_;
    std::string output_directory_ {"."};
    std::string ctgname_fname_ {""};
    std::string rd2rd {""};

    std::string filter_opts_str_ {"l=3000:al=2500:alr=0.6:i=90:aal=5000:aalr=0.5:oh=1000:ohr=0.1:ilid=100"};               //!< overlap filtering options when loading 
    Overlap::Filter filter_opts_ { filter_opts_str_ };
    std::string cov_opts_str_ {""};
    CoverageOptions cov_opts_ {cov_opts_str_};
    std::string phase_opts_str_ {"sc=6:sr=0.66"};
    PhaserOptions phase_opts_ { phase_opts_str_ };
    int thread_size_ { 4 };
    int variant_type_ { 1 };
    int snp_match_length_ { 3 };
    int max_count { 6 };
    std::atomic<size_t> curr_thread_size { 0 };
    int loglevel { 0 };

};

} // namespace fsa
