#pragma once

#include <string>
#include <unordered_map>

#include "sequence.hpp"
#include "utils/program.hpp"

#include "phs_options.hpp"
#include "phs_dataset.hpp"
#include "local_phaser.hpp"

namespace fsa {


class ReadStore;

  
class ContigCall : public Program {
public:
    bool ParseArgument(int argc, const char* const argv[]);
    virtual void Running();
    void Usage();
protected:
    ArgumentParser GetArgumentParser();
    void CheckArguments();

    void FindVariants();
    
    void DumpVariants(Seq::Id c, std::ofstream& of, const std::vector<Variant>& vars);
    void DumpReadInfos(Seq::Id c, std::ofstream& of, const std::unordered_map<ReadOffset, ReadInfo>& rdinfos, const std::vector<Variant>& vars);
    void DumpIgnored(std::ofstream& of, const std::unordered_map<int, std::vector<PhaseItem>> &phased);
    void DumpPositions(const std::string& fname,  const std::vector<Variant>& vars);
    void DumpReadInfos(const std::string& fname, const std::unordered_map<Seq::Id, ReadInfo>& read_infos);
    void SaveVcfHead(std::ofstream& of);
protected:
    PhsOptions opts_;
    PhsDataset dataset_ { opts_ };
};

} // namespace fsa {     
