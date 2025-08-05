#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "utils/program.hpp"
#include "utility.hpp"

#include "kmer_counter.hpp"

namespace fsa {

class KmerBin : public Program {
public:

    
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();

    std::array<size_t,3> CountKmers(size_t k, const std::string& seq, 
        const KmerSet& patkmers, const KmerSet& matkmers,const KmerSet& offkmers);

    std::string OutputPath(const std::string &fname) const { return output_directory_+"/"+fname; }

    size_t CheckKmerSet(const KmerSet& patkmers, const KmerSet& matkmers, const KmerSet& offkmers) const;

protected:
    std::string paternal_;
    std::string maternal_;
    std::string offspring_;
    std::string ifname_;
    std::string ofname_;
    std::string output_directory_ { "." };
    int thread_size_ { 4 };
    int th_count_ { 0 };
    double th_rate_ { 0.0 };
};



} // namespace fsa

