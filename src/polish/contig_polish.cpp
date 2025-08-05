#include "contig_polish.hpp"

#include <edlib.h>
#include <iostream>
#include "./utils/logger.hpp"
#include "utility.hpp"
#include "contig_graph.hpp"

namespace fsa {


ArgumentParser ContigPolish::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}

void ContigPolish::Running() {
    SET_LOG_LEVEL(DEBUG);
    dataset_.Load();

    LOG(INFO)("Start detecting misassemblies");
    AnalyzeContigs();

    LOG(INFO)("Build Contig graph");
    BuildGraph();
}

void ContigPolish::AnalyzeContigs() {
    std::ofstream of_cov("cov_all");
    std::ofstream of_win("cov_win");
    std::ofstream of_mis("mis.bed");
    std::ofstream of_match("match_info");
    std::ofstream of_pol(opts_.cread_fname_);
    std::mutex mutex;

    auto dump = [&mutex, &of_cov, &of_mis, &of_win, &of_pol, &of_match](std::shared_ptr<ContigAnalyzer> ctg_dtr) {
        std::lock_guard<std::mutex> locked(mutex);

        ctg_dtr->DumpCoverage(of_cov);
        ctg_dtr->DumpWindow(of_win);
        ctg_dtr->DumpMatch(of_match);
        ctg_dtr->SaveErrors(of_mis);

    };

    ContigGraph graph;

    auto make_ctg_analyser = [&mutex, this](Seq::Id tid, PolDataset &dataset) -> std::shared_ptr<ContigAnalyzer> {
        std::lock_guard<std::mutex> locked(mutex);
        ctg_analyzers_.emplace_back(new ContigAnalyzer(tid, dataset_));
        return ctg_analyzers_.back();
    };

    std::atomic<size_t> index {0};
    auto work_func = [dump, &index, this, &graph, make_ctg_analyser](size_t _) {
        
        for (size_t i = index.fetch_add(1); i < dataset_.ctg_ids_.size(); i = index.fetch_add(1)) {
            auto tid = dataset_.ctg_ids_[i];    

            auto ctg_analyser = make_ctg_analyser(tid, dataset_);
            ctg_analyser->Detect();
            auto frgs = ctg_analyser->Split();
            graph.AddFragment(frgs);
            dump(ctg_analyser);
        }
    };

    MultiThreadRun((size_t)opts_.thread_size, work_func);

}

void ContigPolish::BuildGraph() {
    for (const auto& ctg_alzr : ctg_analyzers_) {
        graph_.AddFragment(ctg_alzr->Split());
    }
    graph_.Build();
}

void ContigPolish::PolishContigs() {
    auto chains = graph_.GetChains();
    for (auto& chain : chains) {
        chain.Polish();
    }
}


} // namespace fsa {
