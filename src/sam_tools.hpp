#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>

#include "utils/program.hpp"
#include "utility.hpp"
#include "file_io.hpp"
#include "./phase/phase_info.hpp"

namespace fsa {

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
        Add(new Program_Test());
    }
};


} // namespace fsa

