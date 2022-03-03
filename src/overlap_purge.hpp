#ifndef FSA_OVERLAP_PURGE_HPP
#define FSA_OVERLAP_PURGE_HPP


#include <string>

#include "utils/program.hpp"
#include "sequence.hpp"
#include "overlap_store.hpp"
#include "file_io.hpp"
#include "read_store.hpp"

namespace fsa {

class OverlapPurge : public Program {
public:
    virtual ArgumentParser GetArgumentParser();
    virtual void Running();

    void LoadFiles();

    void SaveFiles();
    void GroupReadIds();

    void SaveName(const std::string &opattern, size_t group_size);
protected:

    std::string ifname_;
    std::string ofname_;
    int thread_size_ { 1 };
    std::string read_file_;

    ReadStore rd_store_;
    OverlapStore ol_store_{rd_store_.GetStringPool()};
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> groups_;
 
    std::vector<Seq::Id> read_ids_;
    std::vector<Seq::Id> grouped_ids_;
    std::vector<size_t> group_ticks;
    std::array<size_t, 2> group_size {{50000, 55000}};
};

} // namespace fsa {

#endif // FSA_OVERLAP_PURGE_HPP

