#pragma once

#include <string>
#include <unordered_set>

#include "utils/program.hpp"
#include "utility.hpp"
#include "file_io.hpp"
#include "./phase/phase_info.hpp"

namespace fsa {

class Program_Filter : public Program {
public:
    Program_Filter() {
        name_ = "filter";
        desc_ = "filter out inconsistent overlaps";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofname_, "ofname", "output file");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");
        ap.AddNamedOption(inconsistent_, "inconsistent", "inconsistent overlaps");
        ap.AddNamedOption(consistent_, "consistent", "consistent overlaps");
        ap.AddNamedOption(threshold_, "threshold", "");
        ap.AddNamedOption(filtered_, "filtered", "file to record filtered overlaps");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string ofname_;
    int thread_size_ { 4 };
    std::string inconsistent_;
    std::string consistent_;
    int threshold_ { 500 };
    std::string filtered_;
};

class OverlapTools : public MultiProgram {
public:
    OverlapTools() {
        Add(new Program_Filter());
    }
};


} // namespace fsa

