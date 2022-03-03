#include "read_haplotype.hpp"

#include <unordered_set>
#include <iostream>
#include <atomic>
#include "../utils/logger.hpp"
#include "overlap_store.hpp"
#include "utility.hpp"
#include "contig_phaser.hpp"
//#include "local_phaser.hpp"

namespace fsa {
    
bool ReadHaplotype::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}


void ReadHaplotype::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser ReadHaplotype::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}


void ReadHaplotype::CheckArguments() {
    opts_.CheckArguments();
}


void ReadHaplotype::Running() {
    dataset_.Load();
    
    FindVariants();
}

void ReadHaplotype::FindVariants() {

    std::mutex mutex;

    std::ofstream of_vars(opts_.Varaints());
    std::ofstream of_rdinfos(opts_.Readinfos());
    std::ofstream of_phased(opts_.Inconsistent());
    std::ofstream of_consistent(opts_.Consistent());
    
    auto dump_func = [&](std::ostringstream &oss_vars, std::ostringstream &oss_rdinfos, std::ostringstream &oss_phased, std::ostringstream &oss_consistent) {
        std::lock_guard<std::mutex> lock(mutex);

        of_vars << oss_vars.str();
        oss_vars.str("");

        of_rdinfos << oss_rdinfos.str();
        oss_rdinfos.str("");

        of_phased << oss_phased.str();
        oss_phased.str("");

        of_consistent << oss_consistent.str();
        oss_consistent.str("");
    };

    std::vector<Seq::Id> contigs_list = dataset_.GetSortedContigs();

    std::atomic<size_t> index {0};

    auto work_func = [&](size_t i) {
        opts_.curr_thread_size.fetch_add(1);
        std::ostringstream oss_vars;
        std::ostringstream oss_rdinfos;
        std::ostringstream oss_phased;
        std::ostringstream oss_consistent;

        size_t curr = index.fetch_add(1);
        while (curr < contigs_list.size()) {
            Seq::Id c = contigs_list[curr];
        
            LOG(INFO)("Begin Contig %s, %zd/%zd", dataset_.QueryNameById(c).c_str(), curr, contigs_list.size());

            ContigPhaser phaser(c, dataset_, opts_);
            
            phaser.DumpReadInfos(oss_rdinfos);
            phaser.DumpInconsistent(oss_phased);
            phaser.DumpVariants(oss_vars);
            phaser.DumpConsistent(oss_consistent);
            LOG(INFO)("End Contig %s", dataset_.QueryNameById(c).c_str());
            curr = index.fetch_add(1);
            dump_func(oss_vars, oss_rdinfos, oss_phased, oss_consistent);
        }
        opts_.curr_thread_size.fetch_add(-1);
    };

    MultiThreadRun((size_t)opts_.thread_size_, work_func);
    //work_func(0);
}



} // namespace fsa {
