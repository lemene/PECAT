#include "contig_call.hpp"

#include <unordered_set>
#include <iostream>
#include <atomic>
#include "../utils/logger.hpp"
#include "overlap_store.hpp"
#include "utility.hpp"
#include "contig_phaser.hpp"

namespace fsa {
    
bool ContigCall::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}


void ContigCall::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser ContigCall::GetArgumentParser() {
    ArgumentParser ap;
    opts_.SetArguments(ap);
    return ap;
}


void ContigCall::CheckArguments() {
    opts_.CheckArguments();
}


void ContigCall::Running() {
    dataset_.Load();
    
    FindVariants();
}

void ContigCall::FindVariants() {

    std::mutex mutex;

    std::ofstream of_vars(opts_.Varaints());
    
    auto dump_func = [&](std::ostringstream &oss_vars) {
        std::lock_guard<std::mutex> lock(mutex);

        of_vars << oss_vars.str();
        oss_vars.str("");
    };

    std::vector<Seq::Id> contigs_list = dataset_.GetSortedContigs();

    std::atomic<size_t> index {0};

    auto work_func = [&](size_t i) {
        opts_.curr_thread_size.fetch_add(1);
        std::ostringstream oss_vars;

        size_t curr = index.fetch_add(1);
        while (curr < contigs_list.size()) {
            Seq::Id c = contigs_list[curr];
        
            LOG(INFO)("Begin Contig %s, %zd/%zd", dataset_.QueryNameById(c).c_str(), curr, contigs_list.size());

            ContigPhaser phaser(c, dataset_, opts_);
            
            phaser.Call();
            phaser.SaveVariantInVcf(oss_vars);
            LOG(INFO)("End Contig %s", dataset_.QueryNameById(c).c_str());
            curr = index.fetch_add(1);
            dump_func(oss_vars);
        }
        opts_.curr_thread_size.fetch_add(-1);
    };

    MultiThreadRun((size_t)opts_.thread_size_, work_func);
    //work_func(0);
}



} // namespace fsa {
