#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>

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
        ap.AddNamedOption(rate_, "rate_", "");
        ap.AddNamedOption(strict_, "strict", "");
        ap.AddNamedOption(filtered_, "filtered", "file to record filtered overlaps");
        ap.AddNamedOption(range_fn_, "range", "corrected reads which record its range in raw reads");

        return ap;
    }
    virtual void Running();
    void LoadRanges(const std::string& fname, StringPool& sp);
protected:
    std::string ifname_;
    std::string ofname_;
    int thread_size_ { 4 };
    std::string inconsistent_;
    std::string consistent_;
    int threshold_ { 500 };
    double rate_ { 0.004 };
    std::string filtered_;
    bool strict_ { false };
    std::string range_fn_ { "" };
    std::unordered_map<Seq::Id, std::array<size_t, 2>> ranges_;
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

class Program_Sub : public Program {
public:
    Program_Sub() {
        name_ = "sub";
        desc_ = "purge overlaps";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofname_, "ofnames", "output file");
        ap.AddNamedOption(names_, "names", "");
        ap.AddNamedOption(names_fname_, "names_fname", "");
        ap.AddNamedOption(thread_size_, "thread_size", "number of threads");
        return ap;
    }
    virtual void Running();

protected:
    std::unordered_set<std::string> LoadNames() const;

protected:
    std::string ifname_;
    std::string ofname_;
    std::string names_;
    std::string names_fname_;
    int thread_size_ { 4 };
};

class Program_Accuracy : public Program {
public:
    Program_Accuracy() {
        name_ = "accuracy";
        desc_ = "calculate accuracy of paf";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "paf file");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    int thread_size_ { 4 };
};

class Program_Accuracy2 : public Program {
public:
    Program_Accuracy2() {
        name_ = "accuracy2";
        desc_ = "calculate accuracy of paf";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "paf file");
        ap.AddPositionOption(bed_fname_, "bed", "bed file");
        ap.AddNamedOption(thread_size_, "thread_size", "thread size");

        return ap;
    }
    virtual void Running();
protected:
    std::string ifname_;
    std::string bed_fname_;
    int thread_size_ { 4 };
};

class OverlapTools : public MultiProgram {
public:
    OverlapTools() {
        Add(new Program_Filter());
        Add(new Program_Split());
        Add(new Program_Purge());
        Add(new Program_Sub());
        Add(new Program_Accuracy());
        Add(new Program_Accuracy2());
    }
};


} // namespace fsa

