#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "utils/program.hpp"
#include "utility.hpp"

namespace fsa {

class KmerBin : public Program {
public:
    using KmerId = unsigned long long;
    struct KmerItem { KmerId kmer; int size;};
    struct KmerSet0 {
        bool Empty() const { return kmers.empty(); }
        bool Find(KmerId kid) const  { return kmers.find(kid) != kmers.end(); }
        size_t Size() const { return kmers.size(); }
        size_t k;
        std::unordered_map<KmerId, int> kmers;
    };
    struct KmerSet1 {
        bool Empty() const { return kmers.empty(); }
        bool Find(KmerId kid) const;
        size_t Size() const { return kmers.size(); }
        void BuildIndex();
        size_t k;
        std::vector<std::array<size_t,2>> index;
        std::vector<KmerItem> kmers;
    };
    using KmerSet = KmerSet1;
    
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();

    KmerSet0 LoadKmers0(const std::string &fname);
    KmerSet1 LoadKmers1(const std::string &fname);

    std::array<size_t,3> CountKmers(size_t k, const std::string& seq, 
        const KmerSet& patkmers, const KmerSet& matkmers,const KmerSet& offkmers);

    std::string OutputPath(const std::string &fname) const { return output_directory_+"/"+fname; }

    size_t CheckKmerSet(const KmerSet& patkmers, const KmerSet& matkmers, const KmerSet& offkmers) const;

    KmerId KmerStringToId(const std::string &str);
    size_t GetKmerLength(const std::string &fname);
    size_t CountLines(const std::string &fname);
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

