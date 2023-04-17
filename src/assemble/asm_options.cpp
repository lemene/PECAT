#include "asm_options.hpp"

#include <regex>

namespace fsa {
void AsmOptions::SetArguments(ArgumentParser &ap) {
    ap.AddNamedOption(min_identity, "min_identity", "minimum identity of edges in the string graph");
    ap.AddNamedOption(output_directory, "output_directory", "directory for output files");
    ap.AddNamedOption(thread_size, "thread_size", "number of threads");
    ap.AddNamedOption(filter0_opts, "filter0", "filter options using in loading overlaps");
    ap.AddNamedOption(min_coverage, "min_coverage", "minimum base coverage, negative number = determined by the program.");
    ap.AddNamedOption(coverage, "coverage", "mincov|maxcov|maxdiff, mincov: minimum base coverage, maxcov: maximum base coverage, maxdiff: maximum coverage difference. negative number means the threshold value is determined by the program.");
    ap.AddNamedOption(genome_size, "genome_size", "genome size. It determines the maximum length of reads with coverage together", "INT", ParamToGenomeSize);

    ap.AddNamedOption(read_file, "read_file", "The file are used to increase the speed of loading overlaps.");

    ap.AddNamedOption(max_offset_rate, "max_offset_rate", "");
    ap.AddNamedOption(reduction0, "reduction0", "reduction method0 for the graph");
    ap.AddNamedOption(reduction1, "reduction1", "reduction method1 for the graph");
    ap.AddNamedOption(reducer0, "reducer0", "parameters of reducer0");
    ap.AddNamedOption(reducer1, "reducer1", "parameters of reducer1");

    ap.AddNamedOption(min_contig_length, "min_contig_length", "minimum length of contigs");
    ap.AddNamedOption(contig_format, "contig_format", "\"prialt\" = primary/alterate format, \"dual\" = dual format");
    ap.AddNamedOption(select_branch, "select_branch", "selecting method when encountering branches in the graph, \"no\" = do not select any branch, \"best\" = select the most probable branch", "\"no|best\"");
    ap.AddNamedOption(phased, "phased", "phased read names from the previous steps");
    ap.AddNamedOption(variants, "variants", "variants in reads");
    ap.AddNamedOption(max_bubble_length, "max_bubble_length", "");
    ap.AddNamedOption(max_trivial_length, "max_trivial_length", "maximun trivial branch length, the branch will be removed if its length its length is less than this value");
    ap.AddNamedOption(diploid_rate, "contig_dup_rate", "conting is marked as a duplicate of another conting. If the proportion of its reads inconsistent with other contings exceeds this value.");
    ap.AddNamedOption(hic_info, "hic_info", "information of alignment between hic sequences and the first round of assembly");
    
    // for debuging the program
    ap.AddNamedOption(skip_purge, "skip_purge", "it is for debuging");
    ap.AddNamedOption(debug_name, "debug_name", "");
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

void AsmOptions::UpdateByGenomeSize(size_t gsize) {
    if (genome_size == 0) {
        genome_size = gsize;
    }

    if (max_trivial_length < 0) {
        if (genome_size < 100000000) {
            max_trivial_length = 1000;
        } else if (genome_size < 500000000) {
            max_trivial_length = 10000;
        } else if (genome_size < 1000000000) {
            max_trivial_length = 100000;            
        } else {
            max_trivial_length = 500000;      
        }
    }

    if (max_spur_length < 0) {

    }
    
}

}   // namespace fsa