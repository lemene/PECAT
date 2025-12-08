#pragma once

#include <string>
#include <unordered_set>

#include "utils/program.hpp"
#include "utility.hpp"
#include "utils/string_pool.hpp"

#include "kmer_counter.hpp"

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

class Program_Bin : public Program {
public:
    Program_Bin() {
        name_ = "bin";
        desc_ = "";
    }
    virtual ArgumentParser GetArgumentParser() { 
        ArgumentParser ap("fsa_kmer_tools", "tools about reads", "1.0");
        ap.AddNamedOption(specific_, "specific", "specific kmers of each group, separated by semicolons");
        ap.AddNamedOption(thread_size_, "thread_size", "number of threads");
        ap.AddNamedOption(ifname_, "ifname", "input file name");
        ap.AddNamedOption(ofname_, "ofname", "output file name");
        ap.AddNamedOption(atol_, "atol", "absolute tolerance allowed between two groups");
        ap.AddNamedOption(rtol_, "rtol", "relative tolerance allowed between two groups");
        return ap;
    }
    virtual void Running();

    std::vector<size_t> CountKmers(size_t k, const std::string& seq, const std::vector<KmerStoreUsingVector>& kmers);

protected:
    std::string specific_;
    std::string ifname_;
    std::string ofname_;
    int thread_size_ { 4 };
    int atol_ { 0 };
    double rtol_ { 0.0 };
};


class Program_Graph : public Program {
public:
    Program_Graph() {
        name_ = "graph";
        desc_ = "Build graph from minimizers";
    }
    virtual ArgumentParser GetArgumentParser()
    {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(ranked_kmer_fnames_, "rank_kmers", "");
        ap.AddNamedOption(k_, "k", "k-mer size");
        ap.AddNamedOption(w_, "w", "window size");
        return ap;
    }

    virtual void Running();
    void CountMinimizers();
    std::shared_ptr<class RankedKmers> BuildRankedKmers();

protected:
    std::string ifname_;
    int k_{19};
    int w_{10};
    int thread_size_{4};
    std::string ranked_kmer_fnames_;
    std::shared_ptr<class RankedKmers> ranked_kmers_;
};

class Program_Histo : public Program {
public:
    Program_Histo() {
        name_ = "histo";
        desc_ = "kmer histogram";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddPositionOption(ofname_, "ofname", "output file name");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(freq_fname_, "freq", "k-mer frequency");
        return ap;
    }

    virtual void Running();
protected:
    std::string ifname_;
    std::string ofname_;
    int thread_size_{4};
    std::string freq_fname_;

};


class Program_FreqFreq : public Program {
public:
    Program_FreqFreq() {
        name_ = "freqfreq";
        desc_ = "frequncy of read kmer frequency";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(global_fname_, "global", "read_global_kmer_freq_distribution");
        ap.AddNamedOption(local_fname_, "local", "read_local_kmer_freq_distribution");
        ap.AddNamedOption(freq_fname_, "kmer_freq", "");
        ap.AddNamedOption(k_, "k", "kmer length");
        return ap;
    }

    virtual void Running();
protected:
    int k_ { 17 };
    std::string ifname_;
    std::string ofname_;
    int thread_size_{4};
    std::string freq_fname_;
    std::string global_fname_;
    std::string local_fname_;
    
};

class Program_SegFreq : public Program {
public:
    Program_SegFreq() {
        name_ = "segfreq";
        desc_ = "segment frequency";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "freq", "k-mer frequency file");
        ap.AddPositionOption(ofname_, "ofname", "output file name");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(start_, "start", "start position");
        ap.AddNamedOption(len_, "end", "end position");
        return ap;
    }   
    virtual void Running();
protected:
    std::string ifname_;
    std::string ofname_;
    int thread_size_{4};
    int start_{0};
    int len_ {5};
};

class Program_Verify : public Program {
public:
    Program_Verify() {
        name_ = "verify";
        desc_ = "";
    }
    virtual ArgumentParser GetArgumentParser() { 
        ArgumentParser ap("fsa_kmer_tools", "tools about reads", "1.0");
        ap.AddNamedOption(thread_size_, "thread_size", "number of threads");
        ap.AddPositionOption(freq_fname0_, "ifname0", "input file name");
        ap.AddPositionOption(freq_fname1_, "ifname1", "input file name");
        ap.AddPositionOption(ofname_, "ofname", "output file name");
        return ap;
    }
    virtual void Running();

protected:
    std::string freq_fname0_;
    std::string freq_fname1_;
    std::string ofname_;
    int thread_size_ { 4 };

};


class Program_Test : public Program
{
public:
    Program_Test() {
        name_ = "test";
        desc_ = "test";
    }
    virtual ArgumentParser GetArgumentParser()
    {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddNamedOption(thread_size_, "thread_size", "sequnece file");
        ap.AddNamedOption(ranked_kmer_fnames_, "rank_kmers", "");
        ap.AddNamedOption(k_, "k", "k-mer size");
        ap.AddNamedOption(w_, "w", "window size");
        return ap;
    }

    virtual void Running();
    void Test_CountingKmer();
    void Test_CountingMinimizer();
    void CountMinimizers();
    std::shared_ptr<class RankedKmers> BuildRankedKmers();

protected:
    std::string ifname_;
    int k_{19};
    int w_{10};
    int thread_size_{4};
    std::string ranked_kmer_fnames_;
    std::shared_ptr<class RankedKmers> ranked_kmers_;
};

class KmerTools : public MultiProgram {
public:
    KmerTools() {
        Add(new Program_Count());
        Add(new Program_Gap());
        Add(new Program_Bin());
        Add(new Program_Graph());
        Add(new Program_Histo());
        Add(new Program_FreqFreq());
        Add(new Program_SegFreq());
        Add(new Program_Verify());
        Add(new Program_Test());
    }
}; 
    
} // namespace fsa

