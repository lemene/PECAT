#pragma once

#include <string>
#include <unordered_set>

#include "utils/program.hpp"
#include "utility.hpp"
#include "utils/string_pool.hpp"

namespace fsa {


class Program_Check : public Program {
public:
    Program_Check() {
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

class Program_N50 : public Program {
public:
    Program_N50() {
        name_ = "n50";
        desc_ = "stat N50 of sequences";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddNamedOption(genome_size_, "genome_size", "genome size");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        return ap;
    }

    virtual void Running();
protected:
    std::string ifname_;
    long long genome_size_ { 0 };
};

class Program_Split : public Program {
public:
    Program_Split() {
        name_ = "split";
        desc_ = "split sequence file to smaller files";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "sequnece file");
        ap.AddPositionOption(opattern_, "opattern", "output filename pattern, {} indicates which file it is");
        ap.AddNamedOption(block_size_, "block_size", "size of split files");
        ap.AddNamedOption(base_size_, "base_size", "split longest sequeces");
        return ap;
    }

    virtual void Running();
protected:
    std::string ifname_;
    std::string opattern_;
    long long block_size_ { 0 };
    long long base_size_ { 0 };
    int min_length_ { 0 }; 
};


class Program_SplitName : public Program {
public:
    Program_SplitName() {
        name_ = "split_name";
        desc_ = "split sequence file to smaller files";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(fn_reads_, "reads", "sequnece file");
        ap.AddPositionOption(fn_sub_reads_, "sub_reads", "output filename pattern, {} indicates which file it is");
        ap.AddNamedOption(block_size_, "block_size", "size of split files");
        ap.AddNamedOption(base_size_, "base_size", "split longest sequeces");
        ap.AddNamedOption(thread_size_, "thread_size", "number of threads");
        ap.AddNamedOption(method_, "method", "method to split reads");
        ap.AddNamedOption(fn_ols_, "overlaps", "overlap file");
        ap.AddNamedOption(fn_sub_ols_, "sub_overlaps", "sub files that overlaps are split into");
        return ap;
    }

    virtual void Running();
    void LoadReadnames(const std::string &fname);   
    void GroupReadsPlainly();
    void GroupReadsByOverlaps();
    void Group(const std::vector<std::vector<int>>& groups, const class ReadStore &rd_store);
    void SaveReadnames(const std::string &opattern);
    void SaveOverlaps(const std::string &fn_ols, const std::string &fn_sub_ols);
protected:
    std::string fn_reads_;
    std::string fn_sub_reads_;
    long long block_size_ { 0 };
    long long base_size_ { 0 };
    int min_length_ { 0 }; 
    int thread_size_ { 1 };
    std::string method_ { "plain" };
    std::string fn_ols_;
    std::string fn_sub_ols_;
    StringPool string_pool_;
    std::vector<int> lengths_;
    std::unordered_map<Seq::Id, int> read_ids_index_;
    std::vector<Seq::Id> read_ids_;
    std::vector<int> groups_;

};

class Program_Longest : public Program {
public:
    Program_Longest() {
        name_ = "longest";
        desc_ = "select the longest sequences";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddNamedOption(base_size_, "base_size", "total number of selected bases");
        ap.AddNamedOption(min_length_, "min_length", "mininum length of reads");
        ap.AddNamedOption(id2name_, "id2name", "if it is set, the names of sequences will be replaced with serial numbers");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofname_, "ofname", "output file");
        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string ofname_;
    int min_length_ { 0 }; 
    long long base_size_ { 0 };
    std::string id2name_;
};


class Program_Random : public Program {
public:
    Program_Random() {
        name_ = "random";
        desc_ = "randomly select sequences";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddNamedOption(base_size_, "base_size", "total number of selected bases");
        ap.AddNamedOption(min_length_, "min_length", "mininum length of reads");
        ap.AddNamedOption(id2name_, "id2name", "if it is set, the names of sequences will be replaced with serial numbers");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofname_, "ofname", "output file");
        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string ofname_;
    int min_length_ { 0 }; 
    long long base_size_ { 0 };
    std::string id2name_;
};


class Program_Fasta2Fastq : public Program {
public:
    Program_Fasta2Fastq() {
        name_ = "fa2fq";
        desc_ = "convert fasta to fastq";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(fasta_, "fasta", "sequnece file");
        ap.AddPositionOption(fastq_, "fastq", "sequnece file");
        ap.AddNamedOption(quality_, "quality", "total number of selected bases");
        
        return ap;
    }
protected:
    std::string fasta_;
    std::string fastq_;
    std::string quality_ { "g" };
};


class Program_Sub : public Program {
public:
    Program_Sub() {
        name_ = "sub";
        desc_ = "extract the reads in names";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofname_, "ofname", "output file");
        ap.AddNamedOption(names_, "names", "read names, splitted by a comma");
        ap.AddNamedOption(names_fname_, "names_fname", "read names, splitted by a newline");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string names_;
    std::string names_fname_;
    std::string ofname_;
};


class Program_Test : public Program {
public:
    Program_Test() {
        name_ = "test";
        desc_ = "test some functions";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname0_, "ifname0", "input file");
        ap.AddPositionOption(ifname1_, "ifname1", "output file");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname0_;
    std::string ifname1_;
};

template<typename C>
void FilterReadFile(const std::string &ifname, const std::string &ofname, const std::string &id2name, C check);

class ReadTools : public MultiProgram {
public:
    ReadTools() {
        Add(new Program_Check());
        Add(new Program_N50());
        Add(new Program_Split());
        Add(new Program_SplitName());
        Add(new Program_Longest());
        Add(new Program_Random());
        Add(new Program_Sub());
        Add(new Program_Test());
    }
};

} // namespace fsa

