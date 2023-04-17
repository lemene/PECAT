#include "crr_options.hpp"

#include "../utils/logger.hpp"


namespace fsa {
void CrrOptions::SetArguments(ArgumentParser &ap) {

    ap.AddNamedOption(thread_size, "thread_size", "thread size");
    ap.AddNamedOption(min_coverage, "min_coverage", "");

    ap.AddNamedOption(variants, "variants", "");
    ap.AddNamedOption(use_cache, "use_cache", "use cache for alignment");
    ap.AddNamedOption(debug, "debug", "Output debugging information");
    ap.AddNamedOption(debug_flag, "debug_flag", "Output debugging information");
}

void CrrOptions::CheckArguments() {
}


}   // namespace fsa