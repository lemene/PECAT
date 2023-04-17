#include "overlap_assemble.hpp"

#include <cassert>
#include <algorithm>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include <sstream>
#include <iostream>
#include <atomic>
#include <edlib.h>

#include "../file_io.hpp"
#include "../utility.hpp"


#include "../utils/logger.hpp"
#include "../phase/hic_read_infos.hpp"

namespace fsa {

OverlapAssemble::OverlapAssemble() : string_graph_(dataset_), path_graph_(dataset_) {
}

OverlapAssemble::~OverlapAssemble() {
}


ArgumentParser OverlapAssemble::GetArgumentParser() {
    ArgumentParser ap("fsa_ol_assemble", "Assemble overlaps", "1.0");
    opts_.SetArguments(ap);
    
    return ap;
}

void OverlapAssemble::CheckArguments() {
    opts_.CheckArguments();
    //DUMPER.SetLevel(opts_.dump);
    DUMPER.SetLevel(1);
    DUMPER.SetDirectory(opts_.output_directory);
}


void OverlapAssemble::Running() {
    // if read_file is provided, overlap-loading step can be accelerated

    if (!opts_.skip_purge) {
        dataset_.Load();

        dataset_.Purge();

    } else {
        dataset_.LoadPurged();

    }

    CreateStringGraph();

    CreatePathGraph();

    SaveGraph();

    if (dataset_.GetReadStore().Size() > 0) {
        SaveContigs();
    }
}

void OverlapAssemble::CreateStringGraph() {
    LOG(INFO)("Create string graph");
    
    string_graph_.AddOverlaps(dataset_.GetOverlapStore());

    LOG(INFO)("Simplify string graph");
    string_graph_.Simplify(opts_.reduction0, opts_.reducer0);
    
    LOG(INFO)("Identify paths from string graph");
    string_graph_.IdentifySimplePaths();
}

void OverlapAssemble::CreatePathGraph() {
    LOG(INFO)("Create path graph");
    path_graph_.BuildFrom(string_graph_);
    
    LOG(INFO)("Simplify path graph");
    path_graph_.Simplify(opts_.reduction1, opts_.reducer1);
    
    LOG(INFO)("Identify paths from path graph");
    path_graph_.IdentifyPaths(opts_.select_branch);
}

void OverlapAssemble::SaveContigs() {
    LOG(INFO)("Save Contigs");
    ContigGenerator generator(opts_, dataset_, string_graph_, path_graph_);
    generator.Save();
}


void OverlapAssemble::SaveGraph() {
    
    LOG(INFO)("Save Graph");
    string_graph_.SaveEdges(OutputPath("graph_edges.gz"));
    path_graph_.SaveEdges(OutputPath("graph_paths.gz"));
    path_graph_.SavePaths(OutputPath("contig_paths.gz"));
}


// double OverlapAssemble::ComputeSequenceSimilarity(const std::string &qseq, const std::string &tseq) {
//     auto r = edlibAlign(qseq.c_str(), qseq.size(), tseq.c_str(), tseq.size(),edlibDefaultAlignConfig());
//     double identity = r.status == EDLIB_STATUS_OK ? 1 - r.editDistance *1.0 / qseq.size() : 0;

//     edlibFreeAlignResult(r);
//     return identity;
// }


} // namespace fsa {
