#include "pol_dataset.hpp"

#include <random>

#include "pol_options.hpp"
#include "../correct/align/tool_aligner.hpp"

namespace fsa {

void PolDataset::Load() {
    
    seq_store_.Load(opts_.ctg_fname_, "", true);
    ctg_ids_ = LoadContigIds(opts_.read_name_, opts_.read_name_fname_, seq_store_);

    seq_store_.Load(opts_.rread_fname_, "", true);

    LoadOverlaps(opts_.overlap_fname_);

    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>());

    if (opts_.debug) {
        seq_store_.SaveIdToName("id_2_name");
    }
    
}
    


std::vector<Seq::Id> PolDataset::LoadContigIds(const std::string& name, const std::string &fname, const ReadStore& seq_store) {
    std::vector<Seq::Id> ctg_ids;
    if (!name.empty()) {
        ctg_ids.push_back(seq_store_.QueryIdByName(opts_.read_name_));
    } else if (!fname.empty()) {
        std::ifstream file(opts_.read_name_fname_);
        std::string line;
        while (std::getline(file, line)) {
            ctg_ids.push_back(seq_store_.QueryIdByName(line));
        }

    } else {
        std::array<size_t, 2> range = seq_store.GetIdRange();
        for (int i=range[0]; i<range[1]; ++i) {
            ctg_ids.push_back(i);
        }
    }
    LOG(INFO)("Contig size: %zd", ctg_ids.size());
    return ctg_ids;
}


void PolDataset::LoadOverlaps(const std::string &fname) {
    std::unordered_set<Seq::Id> ids(ctg_ids_.begin(), ctg_ids_.end());

    auto filter = [this, &ids](Overlap &o) {
        bool rel = ids.find(o.b_.id) != ids.end();
        return rel;// && filter0_.ValidQuery(o);

    };
    
    rd_2_ctg_.LoadFast(fname, "", (size_t)opts_.thread_size, filter);

    LOG(INFO)("Load %zd overlaps from file %s", rd_2_ctg_.Size(), fname.c_str());

}


}   // namespace fsa