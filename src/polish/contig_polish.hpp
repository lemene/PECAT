#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <atomic>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "../correct/alignment_graph.hpp"
#include "utils/program.hpp"
#include "pol_dataset.hpp"
#include "pol_options.hpp"
#include "contig_analyzer.hpp"
#include "contig_graph.hpp"

namespace fsa {
using ArrayGraph = AlignmentGraph;

class ContigPolish : public Program {
public:
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();
    virtual void CheckArguments() { opts_.CheckArguments(); }
    
protected:
    void AnalyzeContigs();
    void AnalyzeContigs2();
protected:
    PolOptions opts_;
    PolDataset dataset_ { opts_ };
    std::vector<std::shared_ptr<ContigAnalyzer>> ctg_analyzers_;
    ContigGraph graph_;
};

} // namespace fsa {
    
