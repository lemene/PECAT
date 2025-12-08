#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>

#include "utils/program.hpp"
#include "utility.hpp"
#include "file_io.hpp"
#include "./phase/phase_info.hpp"

#include "prog_coverage_profile.hpp"

namespace fsa {

class Program_NextBaseAccuracy : public Program {
public:
    Program_NextBaseAccuracy() {
        name_ = "next_base_accuracy";
        desc_ = "calculate accuracy of the (n+1)th base given the first n bases are correct";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddNamedOption(thread_size_, "thread_size", "number of threads");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    int thread_size_ { 4 };
};

class Program_MultiCoverage : public Program {
public:
    Program_MultiCoverage() {
        name_ = "mcov";
        desc_ = "calculate multi-coverage from mapping file";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
};

class MapTools : public MultiProgram {
public:
    MapTools() {
        Add(new Program_NextBaseAccuracy());
        Add(new Program_CoverageProfile());
    }
};


} // namespace fsa

