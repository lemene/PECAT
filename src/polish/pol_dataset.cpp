#include "pol_dataset.hpp"

#include <random>

#include "pol_options.hpp"
#include "../correct/align/tool_aligner.hpp"

namespace fsa {

void PolDataset::Load() {
    
    read_store_.Load(opts_.ctg_fname_, "", true);
    LoadReadIds();
    read_store_.Load(opts_.rread_fname_, "", true);

    LoadOverlaps(overlap_fname_);

    

    ol_store_.GroupTarget(groups_, opts_.thread_size);


    if (opts_.debug) {
        read_store_.SaveIdToName("id_2_name");
    }
    
}
    


void PolDataset::LoadReadIds() {
    if (!opts_.read_name_.empty()) {
        read_ids_.push_back(read_store_.QueryIdByName(opts_.read_name_));
    } else if (!opts_.read_name_fname_.empty()) {
        std::ifstream file(opts_.read_name_fname_);
        std::string line;
        while (std::getline(file, line)) {
            read_ids_.push_back(read_store_.QueryIdByName(line));
        }

    } else {
        std::array<size_t, 2> range = read_store_.GetIdRange();
        for (int i=range[0]; i<range[1]; ++i) {
            read_ids_.push_back(i);
        }
    }
}


void PolDataset::LoadOverlaps(const std::string &fname) {
    std::unordered_set<Seq::Id> ids(read_ids_.begin(), read_ids_.end());

    auto filter = [this, &ids](Overlap &o) {
        bool rel = ids.find(o.b_.id) != ids.end();
        return rel;// && filter0_.ValidQuery(o);

    };
    
    ol_store_.LoadFast(fname, "", (size_t)opts_.thread_size, filter);

    LOG(INFO)("Load %zd overlaps from file %s", ol_store_.Size(), fname.c_str());

}

// void PolDataset::LoadMappings() {
//     rd_2_ref_.LoadFileBam(opts_.rd_2_ref_fname_, [](const Overlap &o){return true;}, opts_.thread_size);
//     assert(!opts_.rd_2_ref_fname_.empty());

//     if (!opts_.read_name_.empty() || !opts_.read_name_fname_.empty()) {
//         LOG(INFO)("Start collecting bed");
        
//         auto beds = CollectBedFromBam(read_ids_);
//         LOG(INFO)("BED: %zd", beds.size());
//         rd_2_ref_.LoadFileBam(opts_.rd_2_ref_fname_, [](const Overlap &o){return true;}, opts_.thread_size, beds);
//     } else {
//     }

//     mapping_.BuildIndex();
//     LOG(INFO)("Load rd_2_ref size = %zd", rd_2_ref_.Size());
    

// }



}   // namespace fsa