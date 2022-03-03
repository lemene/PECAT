#pragma once

#include <iostream>
#include <cassert>

#include "argument_parser.hpp"
#include "logger.hpp"


namespace fsa {
    
class Program {
public:
    virtual ~Program() {}
    virtual bool ParseArgument(int argc, char *const argv[]) {
        return  GetArgumentParser().ParseArgument(argc, argv);
    }

    virtual void Run() {
        LOG(INFO)("Start");
        CheckArguments();
        PrintArguments();
        Running();
        LOG(INFO)("End");
    };
    
    virtual void Usage() {
        std::cout << GetArgumentParser().Usage();
    }

    virtual void PrintArguments() {
        LOG(INFO)("Arguments: \n%s", GetArgumentParser().PrintOptions().c_str());
    }

    virtual void CheckArguments() {};
    virtual ArgumentParser GetArgumentParser() { return ArgumentParser(Name(), Description(), ""); }
    virtual void Running() {};

    const std::string& Name() const { return name_; }
    const std::string& Description() const { return desc_; }
protected:
    std::string name_;
    std::string desc_;
};

#define PROGRAM_INSTANCE(cls) \
    int main(int argc, char *argv[]) { \
        fsa::cls runner;   \
        if (runner.ParseArgument(argc, argv)) {    \
            runner.Run();   \
        } else {    \
            runner.Usage(); \
        }   \
        return 0;   \
    }

class MultiProgram : public Program {
public:
    virtual ~MultiProgram() { for (auto p : progs_) delete p; }

    virtual bool ParseArgument(int argc, char *const argv[]) {
        if (argc > 1) {
            for (auto p : progs_) {
                if (p->Name() == argv[1]) {
                    curr_ = p;
                    break;
                }
            }
        }

        if (curr_ != nullptr) {
            return curr_->ParseArgument(argc-1, argv+1);
        } else {
            return false;
        }
    }

    virtual void Run() {
        assert(curr_ != nullptr);
        curr_->Run();
    }

    virtual void Usage() {
        if (curr_ != nullptr) {
            curr_->Usage();
        } else {
            for (auto p : progs_) {
                std::cout << p->Name() << "\t\t" << p->Description() << "\n";
            }
        }
    }

    void Add(Program* p) { progs_.push_back(p); }
protected:
    Program* curr_ {nullptr};
    std::vector<Program*> progs_;
};

} // namespace fsa {

