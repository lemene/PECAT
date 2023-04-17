#pragma once

#include <climits>
#include <sstream>
#include <array>
#include <list>

#include "../utils/program.hpp"

#include "../overlap_store.hpp"
#include "../read_store.hpp"

#include "asm_options.hpp"
#include "string_graph.hpp"
#include "phase/phase_info.hpp"
#include "asm_dataset.hpp"
#include "contig_generator.hpp"

namespace fsa {

class OverlapAssemble : public Program {
public:
    OverlapAssemble();
    virtual ~OverlapAssemble();

    virtual ArgumentParser GetArgumentParser();
    virtual void CheckArguments();
    virtual void Running();

    // misc
    std::string OutputPath(const std::string &fname) const { return opts_.OutputPath(fname); }
public:
    AsmOptions opts_;
    AsmDataset dataset_ { opts_ };
    StringGraph string_graph_;
    PathGraph path_graph_;
public:

    void CreateStringGraph();
    void CreatePathGraph();
    void SaveGraph();
    void SaveContigs();

    std::vector<std::vector<Seq::Tile>> StringEdgesToTiles(const std::list<BaseEdge*> &path);


    void ClassifyContigs();



};

} // namespace fsa {

