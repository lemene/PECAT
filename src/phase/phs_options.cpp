#include "phs_options.hpp"

#include "../utils/logger.hpp"
#include "../utility.hpp"

namespace fsa {
void PhsOptions::SetArguments(ArgumentParser &ap) {
    
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    ap.AddNamedOption(ol_fname_, "ol_fname", "alignment filename");
    ap.AddNamedOption(rd_fname_, "rd_fname", "Read filename");
    ap.AddNamedOption(ctg_fname_, "ctg_fname", "Contig filename");
    ap.AddNamedOption(output_directory_, "output_directory", "");
    ap.AddNamedOption(ctgname_fname_, "ctgname_fname", "contig names that are processed");
    ap.AddNamedOption(filter_opts_str_, "filter", "overlap filtering options");
    ap.AddNamedOption(cov_opts_str_, "coverage", "coverage options to identify candidate SNPs, \"l=40:h=200:s=0.2\" means that [l,h] to coverage range to identity SNP and s is alternative base freqency cutoff");
    ap.AddNamedOption(phase_opts_str_, "phase_options", "Phasing options");
    ap.AddNamedOption(snp_match_length_, "snp_match_length", "");
    ap.AddNamedOption(rd2rd, "rd2rd", "overlaps between reads");
    ap.AddNamedOption(loglevel, "loglevel", "the larger the value(0-4), the more detailed the log. ");
    
}

void PhsOptions::CheckArguments() {

    filter_opts_.From(filter_opts_str_);
    cov_opts_.From(cov_opts_str_);
    phase_opts_.From(phase_opts_str_);
    
    filter_opts_str_ = filter_opts_.ToString();
    cov_opts_str_ = cov_opts_.ToString();
    phase_opts_str_ = phase_opts_.ToString();
}

void PhsOptions::PhaserOptions::From(const std::string &str) {
    auto items = SplitStringByChar(str, ':');

    for (auto &i : items) {
        auto kv = SplitStringByChar(i, '=');
        if (kv[0] == "sc") {
            min_support_count = std::stoi(kv[1]);
        } else if (kv[0] == "sr") {
            min_support_rate = std::stod(kv[1]);
        } else if (kv[0] == "icc") {
            min_inconsist_count = std::stoi(kv[1]);
        } else if (kv[0] == "icr") {
            min_inconsist_rate = std::stod(kv[1]);
        } else if (kv[0] == "ic") {
            min_consist_count = std::stoi(kv[1]);
        } else if (kv[0] == "cr") {
            min_consist_rate = std::stod(kv[1]);
        } else if (kv[0] == "debug") {
            debug = kv[1];
        } else {
            LOG(ERROR)("Unrecoginze phasing option %s", kv[0].c_str());
        }
    }
}

std::string PhsOptions::PhaserOptions::ToString() const {
    Output oput;
    oput.Add("sc", min_support_count);
    oput.Add("sr", min_support_rate);
    oput.Add("icc", min_inconsist_count);
    oput.Add("icr", min_inconsist_rate);
    oput.Add("cc", min_consist_count);
    oput.Add("cr", min_consist_rate);
    oput.Add("debug", debug);
    return oput.ToString();
}

}   // namespace fsa