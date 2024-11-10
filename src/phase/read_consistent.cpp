#include "read_consistent.hpp"

#include <unordered_set>
#include <iostream>
#include <atomic>
#include "../utils/logger.hpp"
#include "overlap_store.hpp"
#include "utility.hpp"
#include "contig_phaser.hpp"
//#include "local_phaser.hpp"

namespace fsa {
    
bool ReadConsistent::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}


void ReadConsistent::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser ReadConsistent::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}


void ReadConsistent::CheckArguments() {
    opts_.CheckArguments();
}


void ReadConsistent::Running() {
    dataset_.Load();
    
    FindVariants();
}

void ReadConsistent::FindVariants() {

    std::mutex mutex;

    std::ofstream of_rdinfos(opts_.Readinfos());
    std::ofstream of_phased(opts_.Inconsistent());
    std::ofstream of_consistent(opts_.Consistent());
    
    auto dump_func = [&](std::ostringstream &oss_rdinfos, std::ostringstream &oss_phased, std::ostringstream &oss_consistent) {
        std::lock_guard<std::mutex> lock(mutex);


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
        std::ostringstream oss_rdinfos;
        std::ostringstream oss_phased;
        std::ostringstream oss_consistent;

        size_t curr = index.fetch_add(1);
        while (curr < contigs_list.size()) {
            Seq::Id c = contigs_list[curr];
        
            LOG(INFO)("Begin Contig %s, %zd/%zd", dataset_.QueryNameById(c).c_str(), curr, contigs_list.size());

            ContigPhaser phaser(c, dataset_, opts_);
            
            phaser.DetectConsistent();
            phaser.DumpReadInfos(oss_rdinfos);
            phaser.DumpInconsistent(oss_phased);
            phaser.DumpConsistent(oss_consistent);
            LOG(INFO)("End Contig %s", dataset_.QueryNameById(c).c_str());
            curr = index.fetch_add(1);
            dump_func(oss_rdinfos, oss_phased, oss_consistent);
        }
        opts_.curr_thread_size.fetch_add(-1);
    };

    MultiThreadRun((size_t)opts_.thread_size_, work_func);
    //work_func(0);
}



} // namespace fsa {
