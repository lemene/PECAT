#pragma once

#include <array>
#include <set>
#include <unordered_map>
#include <vector>

#include "../overlap_store.hpp"
#include "../read_store.hpp"
#include "asm_options.hpp"
#include "../graph.hpp"
#include "stat_readinfo.hpp"
#include "read_variants.hpp"
#include "phase/phase_info.hpp"

namespace fsa {

class AsmDataset {

public:
    AsmDataset(AsmOptions& opts) 
     : opts_(opts)
     , sread_info_(read_infos_) {

    }

    const StringPool& GetStringPool() const { return rd_store_.GetStringPool(); }
    ReadStore& GetReadStore() { return rd_store_; }
    const ReadStore& GetReadStore() const { return rd_store_; }
    OverlapStore& GetOverlapStore() { return ol_store_; }
    const OverlapStore& GetOverlapStore() const { return ol_store_; }

    void Load();
    void LoadPurged();

    void LoadOverlaps(const std::string &fname);
    void LoadOverlaps2(const std::string &fname);
    void LoadOverlapsWithoutLowQuality(const std::string &fname);

    void Purge();
    void FilterLowQuality();
    void FilterLowQuality(int id, const std::unordered_map<int, const Overlap*> &group, std::unordered_set<const Overlap*> &ignored);
    double GetOverlapQuality(const Overlap &ol);
    void GroupAndFilterDuplicate();
    void CheckOverlapEnd();
    void CheckOverlapEnd(int id, const std::unordered_map<int, const Overlap*> &group, std::unordered_set<const Overlap*> &ignored, std::unordered_set<int>& ignReads);
    void ExtendOverlapToEnd();
    void ExtendOverlap();

    void FilterContained();
    bool IsContained(const Overlap& o, std::array<int, 2> &rel);
    void FilterCoverage();

    void FilterConsistency();
    void CalcConsistency(int id, const std::unordered_map<int, const Overlap*> &group, std::unordered_set<const Overlap*> &best, std::unordered_set<const Overlap*> &ignored);
    MatrixGraph CalcConsistencyGraph(int id, const std::vector<const Overlap*> ols);
    void SelectBestExtends(Seq::Id id, const MatrixGraph& graph, const std::vector<std::set<int>>& clusters, const std::vector<const Overlap*>& ols, 
                           std::unordered_set<const Overlap*> &best, std::unordered_set<const Overlap*> &ignored);

    void Debug_PrintGraph(const std::string &name, const MatrixGraph& graph, const std::vector<std::set<int>>& clusters, const std::vector<const Overlap*>& ols);

    double CalcLocalIdentityThreshold(std::vector<std::array<double,2>> &idents, int base_lower_limit, int base_upper_limit);

    void ModifyEnd(const Overlap &o, int maxoh);

    /** \return (type mincov maxcov) */
    std::array<int, 3> AnalyzeCoverage(int id, const std::unordered_map<int, const Overlap*>& group);
    int AnalyzeCoverageType(const std::vector<int>& cov, bool log=false);

    void CoverageConfidencePoints(const std::vector<int>& cov, bool log=false);
    std::array<int,2> CoverageConfidencePoints1(const std::vector<int>& cov, bool log=false);

    int CalcMinCoverage() const;
    std::array<int, 3> CalcCoverageThreshold() const;
  
    bool IsReserved(const Overlap &o) const { return GetOlReason(o).type == OlReason::RS_OK; }
    bool IsReserved(Seq::Id id)       const { auto it = read_infos_.find(id); return it->second.filtered.IsOk(); }
    size_t ReservedSize() const {
        size_t sz = 0;
        for (size_t i=0; i < ol_store_.Size(); ++i) {
            if (IsReserved(ol_store_.Get(i))) sz++;
        }
        return sz;
    }

    void UpdateFilteredRead();
    std::string QueryNameById(Seq::Id id) const { return rd_store_.QueryNameById(id); }
    std::unordered_set<Seq::Id> GetNearbyReads(Seq::Id);
    std::unordered_set<const Overlap*> GetExtendOverlaps(Seq::Id tid, int end) const ;
    std::unordered_set<const Overlap*> GetBackOverlaps(Seq::Id tid, int end) const ;
    const Overlap* QueryOverlap(Seq::Id a, Seq::Id b) const {
        auto iter0  = groups_.find(a);
        if (iter0 != groups_.end()) {
            auto iter1 = iter0->second.find(b);
            return iter1 != iter0->second.end() ? iter1->second : (Overlap*)nullptr;
        }
        return nullptr;
    }

    std::unordered_set<Seq::Id> ReservedReads();
    void ExtendReservedReads();

    /** Record internal state and variables */
    void Dump() const;      
    void DumpOverlaps(const std::string &fname) const;
    void DumpFilteredOverlaps(const std::string &fname) const;
    void DumpReadInfos(const std::string &fname, const std::unordered_map<int, ReadStatInfo> &readInfos) const;

    static bool BetterAlignedLength(const Overlap &a, const Overlap &b) { return a.AlignedLength()*a.identity_ > b.AlignedLength()*b.identity_; }
    static void SetOlReason(const Overlap &o, OlReason rs);
    static OlReason GetOlReason(const Overlap &o);


    ReadVariants* GetReadVariants() { 
        if (read_variants_ == nullptr && !opts_.variants.empty()) {
            read_variants_.reset(new class ReadVariants(string_pool_));
            read_variants_->Load(opts_.variants);
        }
        return read_variants_.get(); 
    }
    PhaseInfoFile* GetInconsistentOverlaps() {
        if (phased_reads_ == nullptr && !opts_.variants.empty()) {
            phased_reads_.reset(new PhaseInfoFile(string_pool_));
            phased_reads_->Load(opts_.phased);
        }
        return phased_reads_.get(); 
    }
    // misc
    std::string OutputPath(const std::string &fname) const { return opts_.OutputPath(fname); }

    const ReadStatInfo& GetReadInfo(Seq::Id id) const {
        auto iter = read_infos_.find(id);
        assert(iter != read_infos_.end());
        return iter->second;
    }
public:
    static int Percentile(const std::vector<int>& data, double percent);
    static int FirstTrough(const std::vector<int>& data, size_t last, size_t k);

    void ScanOverlapsToSelectParams(const std::string &fname);    
    double CalcLocalOverhangThreshold(std::vector<std::array<double,2>> &overhang);

    void EstimateGenomeSize();
   
    AsmOptions& opts_;
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> groups_;

    std::unordered_map<Seq::Id, ReadStatInfo> read_infos_;

    StatReadInfo sread_info_;

    StringPool string_pool_;
    ReadStore rd_store_ { string_pool_ };
    OverlapStore ol_store_ { string_pool_ };
    
    std::shared_ptr<ReadVariants> read_variants_;
    std::shared_ptr<PhaseInfoFile> phased_reads_;
};

}