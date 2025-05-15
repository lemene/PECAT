#pragma once

#include <atomic>
#include <memory>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "overlap/mapping.hpp"

namespace fsa {

class PolOptions;
class Dispatcher;

struct PolDataset {
public:
    PolDataset(PolOptions &opt) : opts_(opt) {}

    void Load();
    const std::string& QueryStringById(Seq::Id id) const { return string_pool_.QueryStringById(id); }    
    const StringPool& GetStringPool() const { return string_pool_; }
protected:
    /** reads in  */
    std::vector<Bed> CollectBedFromBam(const std::vector<Seq::Id>& read_ids);
    void LoadOverlaps(const std::string &fname);
    void LoadMappings();
    
    void LoadReadIds();
    void CalcCoverage();
public:
    PolOptions& opts_;

    StringPool string_pool_;
    ReadStore read_store_ {string_pool_};
    OverlapStore ol_store_{string_pool_ };
    std::unordered_map<int, std::unordered_map<int, std::vector<const Overlap*>>> groups_;
    
    OverlapStore rd_2_ref_ {string_pool_ };
    Mapping mapping_ { rd_2_ref_ };

    
    std::vector<Seq::Id> read_ids_;
};


} // namespace fsa
