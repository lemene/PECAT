#pragma once

#include <climits>
#include <sstream>
#include <array>

#include "../utils/program.hpp"

#include "../overlap_store.hpp"
#include "../read_store.hpp"

#include "asm_options.hpp"
#include "string_graph.hpp"
#include "path_graph.hpp"
#include "phase/phase_info.hpp"
#include "asm_dataset.hpp"

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
protected:
    AsmOptions opts_;
    AsmDataset dataset_ { opts_ };
    StringGraph string_graph_;
    PathGraph path_graph_;
public:

    void CreateStringGraph();
    void CreatePathGraph();
    void SaveGraph();
    void SaveContigs();
    void SaveContigs1();

    std::vector<std::vector<Seq::Tile>> StringEdgesToTiles(const std::list<BaseEdge*> &path);


    void ClassifyContigs();

protected:
    std::vector<Seq::Tile> EdgesToTiles(const std::vector<BaseEdge*> &path);

    std::string ConstructContigS(const std::list<BaseEdge*> &contig);
    std::string ConstructContigStraight(const std::list<BaseEdge*> &contig);
    std::string ConstructContig(const std::list<BaseEdge*> &contig);
    double ComputeSequenceSimilarity(const std::string &qseq, const std::string &tseq);
    std::string EdgeToSeq(const BaseEdge *e);
    Seq::Tile EdgeToTile(const BaseEdge *e);


protected:
    struct Contig {
    public:
        Contig(size_t id, const std::list<PathEdge*> &path, StringGraph &sg);
        bool IsDiploid(const Contig& a, const AsmOptions& opts);
        void PhasedReads(PhaseInfoFile *phased);
        void Save(std::ostream& os, OverlapAssemble& ass, int min_contig_length);
        void SaveTiles(std::ostream& os, const StringGraph& sg);
        void SaveBubbles(std::ostream& os, std::ostream& ftile, OverlapAssemble& ass);
        int Length() const { return std::accumulate(pcontig.begin(), pcontig.end(), 0, [](int a, const BaseEdge*b) {
            return a + b->Length();
            }) ;
        }
        void SetHomo(const Contig& ctg) { homo_ = &ctg; }
        bool IsPrimary() const { return homo_ == nullptr; }
        bool IsCircular() const;
        std::string MainName() const;
        std::string SubName(size_t ibubble, size_t ipath) const;
        std::string Description(size_t len) const;

        size_t id_;
        std::list<BaseEdge*> pcontig;
        std::list<std::pair<PathEdge*, std::vector<std::list<BaseEdge*>>>> acontigs;
        std::unordered_set<int> reads;
        std::unordered_set<int> vreads;
        const Contig* homo_ { nullptr };
    };
    
};

} // namespace fsa {

