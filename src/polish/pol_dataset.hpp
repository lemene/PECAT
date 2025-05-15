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
public:
    PolOptions& opts_;
    
    std::string overlap_fname_;
    std::string rread_fname_;
    std::string ctg_fname_;
    std::string cread_fname_;

    StringPool string_pool_;
    ReadStore read_store_ {string_pool_};
    OverlapStore ol_store_{string_pool_ };
    
    OverlapStore rd_2_ref_ {string_pool_ };
    Mapping mapping_ { rd_2_ref_ };

    
    std::vector<Seq::Id> read_ids_;
};


} // namespace fsa
