/**
 * @brief generate contigs from assembly graph.
 * 
 */
#pragma once
#include <list>
#include <vector>

#include "../read_store.hpp"

#include "asm_options.hpp"
#include "string_graph.hpp"
#include "phase/phase_info.hpp"
#include "asm_dataset.hpp"

namespace fsa {

class ContigGenerator {
public:
    struct Contig {
    public:
        Contig(ContigGenerator* ower, size_t id, const std::list<PathEdge*> &path, StringGraph &sg);
        bool IsDiploid(const Contig& a, const AsmOptions& opts) const;
        double DiploidScore(const Contig& a) const; 
        bool IsOverlapped(const Contig& a, const AsmOptions& opts) const;
        void PhasedReads(PhaseInfoFile *phased);
        void SavePrimary(std::ostream& fctg, std::ostream& ftile, uint8_t hap=0);
        void SaveBubbles(std::ostream& fctg, std::ostream& ftile);
        void SaveBubbles(std::ostream& fctg0, std::ostream& ftile0, std::ostream& fctg1, std::ostream& ftile1, const std::string& format);
        void SaveBubbles(std::ostream& fctg, std::ostream& ftile, int id, const std::vector<std::list<BaseEdge*>>& paths);
        size_t Length() const { 
            return std::accumulate(pcontig.begin(), pcontig.end(), (size_t)0, [](size_t a, const BaseEdge*b) {
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

        void PhaseBubbles(const class HicReadInfos &infos, const ReadVariants& rvs);
        void PhaseBubbles();
        void PhaseHeadTail();
        int PhaseCross(BaseEdge* in_edge_0, BaseEdge *in_edge_1, BaseEdge *out_edge_0, BaseEdge *out_edge_1);

        typedef std::unordered_map<int, std::unordered_map<int, std::vector<std::array<int, 2>>>> Variants;
        std::array<Variants,2> ExtendOutSnps(BaseNode* in_node_0, BaseNode* in_node_1, BaseNode* out_node_0, BaseNode* out_node_1);
        std::array<Variants,2> ExtendInSnps(BaseNode* in_node_0, BaseNode* in_node_1, BaseNode* out_node_0, BaseNode* out_node_1);
        std::unordered_set<Seq::Id> CollectOutSnps(Variants& vars, BaseNode* node, BaseNode* altnode);
        std::unordered_set<Seq::Id> CollectOutSnps(Variants& vars, BaseNode* node);
        std::unordered_set<Seq::Id> CollectInSnps(Variants& vars, BaseNode* node, BaseNode* altnode);
        std::unordered_set<Seq::Id> CollectInSnps(Variants& vars, BaseNode* node);
        //std::array<size_t,2> Similarity(const Variants &vars0, const Variants &vars1);
        std::list<BaseEdge*> ConstructPrimaryPath(uint8_t hap);
        std::string MainName() const;
        std::string SubName(size_t ibubble, size_t ipath) const;
        std::string Description(size_t len) const;

        size_t id_;
        std::list<BaseEdge*> pcontig;
        std::vector<std::vector<std::list<BaseEdge*>>> segs_;
        std::list<std::pair<PathEdge*, std::vector<std::list<BaseEdge*>>>> acontigs;
        std::unordered_set<int> reads;
        std::unordered_set<int> vreads;
        const Contig* pri_ { nullptr };
        std::vector<const Contig*> alts_;
        size_t weight_ = 0;
        bool trivial_ = false;
        const std::list<PathEdge*> *path_;
        std::vector<uint8_t> phasing_;
        ContigGenerator* owner_;
        int hap_ = { 0 };
    };

    struct Group {  // Chromosome 
        std::vector<Contig*> ctgs;
    };
public:
    ContigGenerator(AsmOptions &opts, AsmDataset &dataset, StringGraph &string_graph, PathGraph &path_graph)
     : opts_(opts), dataset_(dataset), string_graph_(string_graph), path_graph_(path_graph)
    {}

    void IdentifyContigs();
    void FindInconsistentOverlap(const Contig &a, const Contig &b);
    void GroupContigs();
    bool IsTrivialContig(const Contig& contig);
    void Save();
    void CollectContigFromGraph();
    void PhaseBubbles(std::vector<Contig> &contigs);
    void SaveContigs(std::vector<Contig>& contigs, const std::string& format);

    std::string OutputPath(const std::string &fname) const { return opts_.OutputPath(fname); }

    std::vector<Seq::Tile> EdgesToTiles(const std::vector<BaseEdge*> &path);

    std::string ConstructContigS(const std::list<BaseEdge*> &contig);
    std::string ConstructContigStraight(const std::list<BaseEdge*> &contig);
    std::string ConstructContig(const std::list<BaseEdge*> &contig);
    std::string EdgeToSeq(const BaseEdge *e);
    Seq::Tile EdgeToTile(const BaseEdge *e);


protected:
public:
    AsmOptions &opts_;
    AsmDataset &dataset_;
    StringGraph &string_graph_;
    PathGraph &path_graph_;

    std::vector<Contig> contigs;
    std::vector<Group>  groups_;

};

} // namespace fsa {

