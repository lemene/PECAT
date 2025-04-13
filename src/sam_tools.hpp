#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>

#include "utils/program.hpp"
#include "utility.hpp"
#include "file_io.hpp"
#include "./phase/phase_info.hpp"


namespace fsa {

class Program_GroupRead : public Program {
public:
    Program_GroupRead() {
        name_ = "group";
        desc_ = "";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "bam file");
        ap.AddPositionOption(opattern_, "opattern", "read name");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");
        ap.AddNamedOption(block_size_, "block_size", "block size");

        return ap;
    }
    virtual void Running();
protected:
    void Save(int id, const std::vector<std::string> &names);
protected:
    std::string ifname_;
    std::string opattern_;
    int thread_size_ { 4 };
    long long block_size_ { 4000000000L };

};

class Program_N50 : public Program {

    public:
    Program_N50() {
        name_ = "n50";
        desc_ = "for testing";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "paf file");
        ap.AddNamedOption(genome_size_, "genome_size", "genome size");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string rname_;
    long long genome_size_ { 0 };
    int thread_size_ { 4 };
};


class Program_Stat : public Program {
    struct Position {
        uint16_t match {0};
        uint16_t mismatch[4] = {0,0,0,0};
        uint16_t deletion {0};
        uint16_t insertion {0};
        uint32_t inssize {0};
    };
public:
    Program_Stat() {
        name_ = "stat";
        desc_ = "for testing";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "paf file");
        ap.AddPositionOption(ref_fname_, "reference", "paf file");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
protected:
    void InitializeTables();
    void PrintAccuracy();
protected:
    std::string ifname_;
    std::string ref_fname_;
    long long genome_size_ { 0 };
    int thread_size_ { 4 };
    std::unordered_map<std::string, std::vector<Position>> tables_;
};

class Program_Test : public Program {
public:
    Program_Test() {
        name_ = "test";
        desc_ = "for testing";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "paf file");
        ap.AddPositionOption(rname_, "rname", "read name");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string rname_;
    int thread_size_ { 4 };
};


class SamTools : public MultiProgram {
public:
    SamTools() {
        Add(new Program_GroupRead());
        Add(new Program_N50());
        Add(new Program_Stat());
        Add(new Program_Test());
    }
};


} // namespace fsa

