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

class Program_Split : public Program {
public:
    Program_Split() {
        name_ = "split";
        desc_ = "split overlaps to files according to the names";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofnames_, "ofnames", "output file");
        ap.AddPositionOption(namesets_, "namesets", "names by which the overlaps are split");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
    std::unordered_set<std::string> LoadNameset(const std::string &fname);
protected:
    std::string ifname_;
    std::string ofnames_;
    std::string namesets_;
    int thread_size_ { 4 };
};

class Program_Purge : public Program {
public:
    Program_Purge() {
        name_ = "purge";
        desc_ = "purge overlaps";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofnames_, "ofnames", "output file");
        ap.AddPositionOption(tiles_, "tiles", "tiles of contigs");
        ap.AddPositionOption(readinfos_, "readinfos", "tiles of contigs");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
    std::unordered_set<std::string> LoadNameset(const std::string &fname);
protected:
    std::string ifname_;
    std::string ofnames_;
    std::string tiles_;
    std::string readinfos_;
    int thread_size_ { 4 };
};

class OverlapTools : public MultiProgram {
public:
    OverlapTools() {
        Add(new Program_Filter());
        Add(new Program_Split());
        Add(new Program_Purge());
    }
};


} // namespace fsa

