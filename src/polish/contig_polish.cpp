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
    AnalyzeContigs2();
}

// void ContigPolish::AnalyzeContigs() {
//     std::ofstream of_mcov("multi_covs");
//     std::mutex mutex;

//     auto dump = [&mutex, &of_mcov](std::shared_ptr<ContigAnalyzer> ctg_dtr) {
//         std::lock_guard<std::mutex> locked(mutex);
//         ctg_dtr->DumpMultiCoverage(of_mcov);

//     };

//     ContigGraph graph;

//     auto make_ctg_analyser = [&mutex, this](Seq::Id tid, PolDataset &dataset) -> std::shared_ptr<ContigAnalyzer> {
//         std::lock_guard<std::mutex> locked(mutex);
//         ctg_analyzers_.emplace_back(new ContigAnalyzer(tid, dataset_));
//         return ctg_analyzers_.back();
//     };

//     std::atomic<size_t> index {0};
//     auto work_func = [dump, &index, this, &graph, make_ctg_analyser](size_t _) {
        
//         for (size_t i = index.fetch_add(1); i < dataset_.ctg_ids_.size(); i = index.fetch_add(1)) {
//             auto tid = dataset_.ctg_ids_[i];    

//             auto ctg_analyser = make_ctg_analyser(tid, dataset_);
//             ctg_analyser->Detect();
//             dump(ctg_analyser);
//         }
//     };

//     MultiThreadRun((size_t)opts_.thread_size, work_func);

// }

void ContigPolish::AnalyzeContigs2() {

        
    for (size_t i = 0; i < dataset_.ctg_ids_.size(); i++) {
        auto tid = dataset_.ctg_ids_[i];    
        std::ofstream of_mcov(opts_.OutputPath(dataset_.QueryStringById(tid) + ".mcov"));
        ContigAnalyzer ctg_analyser(tid, dataset_);
        ctg_analyser.ComputeCoverage(opts_.thread_size);
        ctg_analyser.DumpMultiCoverage(of_mcov);
    }


}

} // namespace fsa {
