#include "pol_options.hpp"

#include "../utils/logger.hpp"
#include "../utility.hpp"
#include <sstream>

namespace fsa {
void PolOptions::SetArguments(ArgumentParser &ap) {

    ap.AddPositionOption(overlap_fname_, "ol_fname", "overlap file name");
    ap.AddPositionOption(ctg_fname_, "ctg_fname", "contig file name");
    ap.AddPositionOption(rread_fname_, "rr_fname", "raw read file name");
    ap.AddPositionOption(cread_fname_, "cr_fname", "corrected read file name");

    ap.AddNamedOption(filter0_opts_, "filter0", "overlap filtering options using at loading step", "");
    ap.AddNamedOption(filter1_opts_, "filter1", "overlap filtering options using at loading step", "");
    ap.AddNamedOption(rd_2_ref_fname_, "rd_2_ref", "alignments between reads and reference", "");
    
    ap.AddNamedOption(read_name_, "read_name", "read name for correcting");
    ap.AddNamedOption(read_name_fname_, "read_name_fname", "Set read name for correcting");

    ap.AddNamedOption(graph_fname_, "graph_fname", "The file recording graph");
    ap.AddNamedOption(infos_fname_, "infos_fname", "The file recording score infos");

    ap.AddNamedOption(min_identity_, "min_identity", "");
    ap.AddNamedOption(min_local_identity_, "min_local_identity", "");
    ap.AddNamedOption(check_local_identity_, "check_local_identity", "");
    
    ap.AddNamedOption(aligner_, "aligner", "method for local alignment, diff|edlib.");
    ap.AddNamedOption(score_, "score", "");
    
    ap.AddNamedOption(thread_size, "thread_size", "thread size");
    ap.AddNamedOption(min_coverage, "min_coverage", "");

    ap.AddNamedOption(variants, "variants", "");
    ap.AddNamedOption(use_cache, "use_cache", "use cache for alignment");
    ap.AddNamedOption(debug, "debug", "Output debugging information");
    ap.AddNamedOption(skip_branch_check, "skip_branch_check", "skip branch check");
    ap.AddNamedOption(output_directory, "output_directory", "The directory for temporary files");
}

void PolOptions::CheckArguments() {
    
    filter0_.From(filter0_opts_);
    filter1_.From(filter1_opts_);
    
    filter0_opts_ = filter0_.ToString();
    filter1_opts_ = filter1_.ToString();

    if (debug) SetDebug();
}


}   // namespace fsa