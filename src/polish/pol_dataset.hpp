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
    double GetOverlapQualityThreshold() const { return overlap_quality_threshold_; }
    double GetLocalQualityThreshold() const { return local_quality_threshold_; }
    size_t CountReadMap(Seq::Id id) const ;
    size_t MaxReadLength() const { return max_read_length_; }
protected:
    void Stat();
    /** reads in  */
    std::vector<Bed> CollectBedFromBam(const std::vector<Seq::Id>& read_ids);
    void LoadOverlaps(const std::string &fname);
    void LoadMappings();
    
    std::vector<Seq::Id> LoadContigIds(const std::string& name, const std::string &fname, const ReadStore& store);
    void SelectBestMapping();
public:
    PolOptions& opts_;

    StringPool string_pool_;
    ReadStore seq_store_ {string_pool_};
    std::array<size_t,2> rd_ids_;
    
    OverlapStore rd_2_ctg_ {string_pool_ };
    OverlapGrouper grouper_ { rd_2_ctg_ };
    
    
    std::vector<Seq::Id> ctg_ids_;
    double overlap_quality_threshold_ {0.0};    // TODO Move it to pol_options.hpp
    double local_quality_threshold_ {0.0};
    size_t max_read_length_ {0};
    size_t ave_read_length_ {0};
};


} // namespace fsa
