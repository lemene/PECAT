#include "asm_options.hpp"

#include <regex>

namespace fsa {
void AsmOptions::SetArguments(ArgumentParser &ap) {
    //ap.AddNamedOption(bestn, "bestn", " output best n overlaps on 5' or 3' end for each read");
    ap.AddNamedOption(min_identity, "min_identity", "minimum identity of edges in the string graph");
    ap.AddNamedOption(output_directory, "output_directory", "directory for output files");
    ap.AddNamedOption(thread_size, "thread_size", "number of threads");
    ap.AddNamedOption(filter0_opts, "filter0", "filter options using in loading overlaps");
    ap.AddNamedOption(min_coverage, "min_coverage", "minimum base coverage, negative number = determined by the program.");
    ap.AddNamedOption(genome_size, "genome_size", "genome size. It determines the maximum length of reads with coverage together", "INT", ParamToGenomeSize);

    ap.AddNamedOption(read_file, "read_file", "The file are used to increase the speed of loading overlaps.");

    ap.AddNamedOption(max_offset_rate, "max_offset_rate", "");
    ap.AddNamedOption(debug_name, "debug_name", "");
    ap.AddNamedOption(reduction0, "reduction0", "reduction method0 for the graph");
    ap.AddNamedOption(reduction1, "reduction1", "reduction method1 for the graph");
    ap.AddNamedOption(reducer0, "reducer0", "parameters of reducer0");
    ap.AddNamedOption(reducer1, "reducer1", "parameters of reducer1");

    ap.AddNamedOption(min_contig_length, "min_contig_length", "minimum length of contigs");
    ap.AddNamedOption(select_branch, "select_branch", "selecting method when encountering branches in the graph, \"no\" = do not select any branch, \"best\" = select the most probable branch", "\"no|best\"");
    ap.AddNamedOption(phased, "phased", "phased read names from the previous steps");
    ap.AddNamedOption(variants, "variants", "variants in reads");
    ap.AddNamedOption(max_bubble_length, "max_bubble_length", "");
    ap.AddNamedOption(skip_purge, "skip_purge", "it is for debuging");
    ap.AddNamedOption(dump, "dump", "dump intermediate data, 0: dump nothing, 4: dump all data");

    ap.AddPositionOption(ifname, "overlaps", "input filename");
}



bool AsmOptions::ParamToGenomeSize(const std::string& str, long long *v) {
    
    std::regex pattern("(\\d+)([gGmMkK]?)");
    std::smatch m;
    bool r = std::regex_match(str, m, pattern);
    if (r) {
       *v = atoll(m.str(1).c_str());
       *v *= m.str(2) == "k" || m.str(2) == "K" ? 1024 :
             m.str(2) == "m" || m.str(2) == "M" ? 1024*1024 :
             m.str(2) == "g" || m.str(2) == "G" ? 1024*1024*1024 : 1;
    }
    return r;
}
}   // namespace fsa