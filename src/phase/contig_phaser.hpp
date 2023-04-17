#pragma once

#include <unordered_map>
#include <unordered_set>

#include "../sequence.hpp"

#include "utility.hpp"
#include "phs_dataset.hpp"

namespace fsa {

class ReadStore;
class OverlapStore;
class PhsOptions;
class PhsDataset;

class ContigPhaser {
public:
    ContigPhaser(Seq::Id ctg, const PhsDataset& dataset, PhsOptions &opts);

    void Phase();

    std::unordered_map<ReadOffset, ReadInfo> CollectReads(Seq::Id ctg);
    void FindVariantsInContig(Seq::Id c, const std::unordered_map<ReadOffset, ReadInfo>& read_infos, std::vector<Variant> &vars);
    void ConfirmVariants(std::vector<Variant> &vars);
    void GetVariantsFromSnps(const std::unordered_map<ReadOffset, ReadInfo>& read_infos, std::vector<Variant> &vars);
    void FindVariantsInReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos, const std::vector<Variant> &vars);
    // return ignored overlaps {iid, jid}
    std::unordered_map<ReadOffset, std::unordered_set<ReadOffset>> GroupReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos, int shared);
    void ClassifyReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos);   
    bool IsValidPair(const Overlap &irio, const Overlap& jrio);

    void ScanContig();

    void DumpVariants(std::ostream& of) const;
    void DumpReadInfos(std::ostream& of) const;
    void DumpInconsistent(std::ostream& of) const;
    void DumpConsistent(std::ostream& of) const;
protected:
    const ReadStore& GetReadStore() const { return dataset_.rd_store_; }

protected:
    Seq::Id ctg_;
    const PhsDataset &dataset_;
    const ReadStore &rd_store_;
    const OverlapStore &ol_store_;
    PhsOptions &opts_;

    std::vector<Variant> variants_;
    std::unordered_map<ReadOffset, ReadInfo> read_infos_;
    std::unordered_map<ReadOffset, std::vector<PhaseItem>> consistent_;
    std::unordered_map<ReadOffset, std::vector<PhaseItem>> inconsistent_;
};

} // namespace fsa