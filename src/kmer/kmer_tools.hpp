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


class KmerTools : public MultiProgram {
public:
    KmerTools() {
        Add(new Program_Count());
    }
};
    
} // namespace fsa

