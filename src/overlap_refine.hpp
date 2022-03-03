#pragma once

#include <string>
#include <memory>

#include "utils/program.hpp"
#include "sequence.hpp"
#include "overlap_store.hpp"
#include "file_io.hpp"
#include "read_store.hpp"

namespace fsa {

class OverlapRefine : public Program {
public:
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();
    virtual void CheckArguments();

    Reader* GetReader(const std::string &fname);
    Writer* GetWriter(const std::string &fname);
    bool SplitPafLine(const std::string &line, Overlap& ol);

    bool DetectFileType();

    void TaskExtend(Overlap& o, const StringPool::NameId &ni);
    void TaskIdentity(Overlap& o, const StringPool::NameId &ni);
    void TaskIdle(Overlap& o, const StringPool::NameId &ni) { return; }
    void (OverlapRefine::* GetTask(const std::string& task))(Overlap& o, const StringPool::NameId &ni) ;
protected:

    std::string ifname_;
    std::string ofname_;
    std::string itype_;
    std::string otype_;
    std::string read_fname_;
    int thread_size_ { 1 };
    std::string task_ { "" };

    std::string filter0_opts_;
    std::string filter1_opts_;
    Overlap::Filter filter0_;
    Overlap::Filter filter1_;
    std::string aligner_opts_ {"diff:s=500:e=0.1"};

    ReadStore rd_store_;
    bool (*fromLine)(const std::string& line, Overlap& ol, StringPool::NameId &ni, int &replen);
    std::string (*toLine)(const Overlap &ol, const StringPool::NameId& ni);
};

} // namespace fsa {


