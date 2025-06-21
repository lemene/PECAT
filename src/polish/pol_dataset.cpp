#include "pol_dataset.hpp"

#include <random>

#include "pol_options.hpp"
#include "../correct/align/tool_aligner.hpp"

namespace fsa {

void PolDataset::Load() {
    
    seq_store_.Load(opts_.ctg_fname_, "", true);
    size_t id_start = seq_store_.GetIdUp();
    ctg_ids_ = LoadContigIds(opts_.read_name_, opts_.read_name_fname_, seq_store_);

    seq_store_.Load(opts_.rread_fname_, "", true);
    rd_ids_ = {id_start, seq_store_.GetIdUp()};

    seq_store_.SaveIdToName("id2name");
    LoadOverlaps(opts_.overlap_fname_);

    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>());

    if (opts_.debug) {
        seq_store_.SaveIdToName("id_2_name");
    }
    SelectBestMapping();
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
        return rel && opts_.filter0_.ValidQuery(o);

    };
    
    rd_2_ctg_.LoadFast(fname, "", (size_t)opts_.thread_size, filter);

    LOG(INFO)("Load %zd overlaps from file %s", rd_2_ctg_.Size(), fname.c_str());

}

void PolDataset::SelectBestMapping() {
    struct Item {
        const Overlap* ol;
        double wt;
        uint32_t rank;
    };
    for (size_t ri = rd_ids_[0]; ri < rd_ids_[1]; ++ri) {
        auto gp = grouper_.Get(ri);
        std::vector<Item> items;
        for (size_t i = 0; i < gp.Size(); i++) {
            for (size_t j = 0; j < gp.Size(i); j++) {
                auto ol = gp.Get(i,j);
                //items.push_back({ol, ol->Identity()*ol->QueryLength(), 0});
                items.push_back({ol, ol->Identity()*ol->AlignedLength(), 0});
            }
        }
        if (items.size() > 0) {
            std::vector<int> ranks(items.size());
            std::iota(ranks.begin(), ranks.end(), 0);
            std::random_shuffle(ranks.begin(), ranks.end());
            for (size_t i = 0; i < items.size(); ++i) {
                items[i].rank = ranks[i];
            }
            auto mx = std::max_element(items.begin(), items.end(), [](const Item& a, const Item& b) {
                // if (a.wt + 500 < b.wt) {
                //     return true;
                // } else if (a.wt > b.wt + 500) {
                //     return false;
                // } else {
                //     return a.rank > b.rank;
                // }
                return a.wt < b.wt;
            });
            mx->ol->attached = 1;

            for (auto &it : items) {
                if (it.wt + 500 > mx->wt && it.wt < mx->wt + 500) {
                    mx->ol->attached = 1;
                }

            }

            if (items.size() > 0) {
                LOG(INFO)("SEL %s %.02f %zd", OverlapStore::ToPafLine(*mx->ol, StringPool::TempNameId2(string_pool_)).c_str(), 
                    mx->wt, mx->rank);   
                for (auto it : items) {
                    LOG(INFO)("SEL_all %s %.02f %zd", OverlapStore::ToPafLine(*it.ol, StringPool::TempNameId2(string_pool_)).c_str(), 
                    it.wt, it.rank);   

                } 
            }
        }
    }
}

}   // namespace fsa