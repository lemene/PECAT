#include "contig_fragment.hpp"
#include "contig_analyzer.hpp"

namespace fsa {
ContigFragment::ContigFragment(ContigAnalyzer* ctg_analyzer, size_t start, size_t end)
 : ctg_analyzer_(ctg_analyzer), start_(start), end_(end) {
    // Initialize other members if needed   
}
    

Seq::Id ContigFragment::ContigId() const {
    return ctg_analyzer_->GetId(); 
}
}