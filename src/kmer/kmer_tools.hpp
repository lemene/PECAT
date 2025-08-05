#pragma once

#include <string>
#include <unordered_set>

#include "utils/program.hpp"
#include "utility.hpp"
#include "utils/string_pool.hpp"

namespace fsa {


class Program_Count : public Program {
public:
    Program_Count() {
        name_ = "check";
        desc_ = "check sequences";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        return ap;
    }

    virtual void Running();
protected:
    std::string ifname_;
};

class Program_Test : public Program {
public:
    Program_Test() {
        name_ = "test";
        desc_ = "test";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(ranked_kmer_fnames_, "rank_kmers", "");
        return ap;
    }
 
    virtual void Running(); 
    void Test_CountingKmer(); 
    void Test_CountingMinimizer();
    std::shared_ptr<class RankedKmers> BuildRankedKmers();

protected:
    std::string ifname_;
    uint8_t k_ { 19 };
    uint8_t w_ { 10 };
    int thread_size_ { 4 };
    std::string ranked_kmer_fnames_;
    std::shared_ptr<class RankedKmers> ranked_kmers_;
};

class Program_Gap : public Program {
public:
    Program_Gap() {
        name_ = "gap";
        desc_ = "gap between err kmer";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(ranked_kmer_fnames_, "rank_kmers", "");
        return ap;
    }
 
    virtual void Running(); 
protected:
    std::string ifname_;
    uint8_t k_ { 19 };
    uint8_t w_ { 50 };
    int thread_size_ { 4 };
    std::string ranked_kmer_fnames_;
};

class KmerTools : public MultiProgram {
public:
    KmerTools() {
        Add(new Program_Count());
        Add(new Program_Test());
        Add(new Program_Gap());
    }
}; 
    
} // namespace fsa

