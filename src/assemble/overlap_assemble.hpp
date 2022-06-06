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

namespace fsa {

class OverlapAssemble : public Program {
public:
    OverlapAssemble();
    virtual ~OverlapAssemble();

    virtual ArgumentParser GetArgumentParser();
    virtual void CheckArguments();
    virtual void Running();

    PhaseInfoFile *GetPhasedInfo() { return dataset_.GetInconsistentOverlaps(); }
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
        bool IsDiploid(const Contig& a, const AsmOptions& opts) const;
        bool IsOverlapped(const Contig& a, const AsmOptions& opts) const;
        void PhasedReads(PhaseInfoFile *phased);
        void Save(std::ostream& os, OverlapAssemble& ass, int min_contig_length);
        void SaveTiles(std::ostream& os, const StringGraph& sg);
        void SaveBubbles(std::ostream& fctg, std::ostream& ftile, OverlapAssemble& ass);
        void SaveBubbles(std::ostream& fctg0, std::ostream& ftile0, std::ostream& fctg1, std::ostream& ftile1, OverlapAssemble& ass);
        void SaveBubbles(std::ostream& fctg, std::ostream& ftile, OverlapAssemble& ass, int id, const std::vector<std::list<BaseEdge*>>& paths);
        int Length() const { return std::accumulate(pcontig.begin(), pcontig.end(), 0, [](int a, const BaseEdge*b) {
            return a + b->Length();
            }) ;
        }
        int Weight() const { return weight_; }
        void SetPrimary(const Contig& ctg) { pri_ = &ctg; }
        const Contig* GetPrimary() const { return pri_; }
        void AddAlternate(const Contig& ctg) { alts_.push_back(&ctg); }
        bool IsPrimary() const { return pri_ == nullptr; }
        bool IsCircular() const;
        bool IsCovered(const std::vector<std::list<BaseEdge*>>& paths) const;
        bool IsTrivial() const { return trivial_; }
        void SetTrivial() { trivial_ = true; }
        std::string MainName() const;
        std::string SubName(size_t ibubble, size_t ipath) const;
        std::string Description(size_t len) const;

        size_t id_;
        std::list<BaseEdge*> pcontig;
        std::list<std::pair<PathEdge*, std::vector<std::list<BaseEdge*>>>> acontigs;
        std::unordered_set<int> reads;
        std::unordered_set<int> vreads;
        const Contig* pri_ { nullptr };
        std::vector<const Contig*> alts_;
        size_t weight_ = 0;
        bool trivial_ = false;
        const std::list<PathEdge*> *path_;
    };
    
};

} // namespace fsa {

