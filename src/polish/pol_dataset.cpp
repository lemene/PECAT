#include "pol_dataset.hpp"

#include <random>

#include "pol_options.hpp"
#include "../correct/align/tool_aligner.hpp"

namespace fsa {

void PolDataset::Load() {
    read_store_.Load(opts_.ctg_fname_, "", true);

    LoadOverlaps();

    read_store_.Load(opts_.rread_fname_, "", true);


    if (opts_.debug) {
        read_store_.SaveIdToName("id_2_name");
    }
    
}
    

void PolDataset::LoadOverlaps() {
    const std::string& fname = opts_.overlap_fname_;
    
    ol_store_.Load(fname, "", (size_t)opts_.thread_size, [this](Overlap &o) {
        return opts_.filter0_.Valid(o);
    });

    LOG(INFO)("Load %zd overlaps from file %s", ol_store_.Size(), fname.c_str());
   
}

void PolDataset::LoadMappings() {
    rd_2_ref_.LoadFileBam(opts_.rd_2_ref_fname_, [](const Overlap &o){return true;}, opts_.thread_size);
    assert(!opts_.rd_2_ref_fname_.empty());

    if (!opts_.read_name_.empty() || !opts_.read_name_fname_.empty()) {
        LOG(INFO)("Start collecting bed");
        
        auto beds = CollectBedFromBam(read_ids_);
        LOG(INFO)("BED: %zd", beds.size());
        rd_2_ref_.LoadFileBam(opts_.rd_2_ref_fname_, [](const Overlap &o){return true;}, opts_.thread_size, beds);
    } else {
    }

    mapping_.BuildIndex();
    LOG(INFO)("Load rd_2_ref size = %zd", rd_2_ref_.Size());
    

}



}   // namespace fsa