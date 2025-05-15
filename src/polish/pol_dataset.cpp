#include "pol_dataset.hpp"

#include <random>

#include "pol_options.hpp"
#include "../correct/align/tool_aligner.hpp"

namespace fsa {

void PolDataset::Load() {
    
    read_store_.Load(opts_.ctg_fname_, "", true);
    LoadReadIds();
    read_store_.Load(opts_.rread_fname_, "", true);

    LoadOverlaps(opts_.overlap_fname_);
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
    LOG(INFO)("Contig size: %zd", read_ids_.size());
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


void PolDataset::CalcCoverage() {
    for (auto ctgid : read_ids_) {
        const auto ctglen = read_store_.GetSeqLength(ctgid);
        std::vector<int> cov(ctglen+1, 0);

        auto ctg_ols = groups_.find(ctgid);
        if (ctg_ols != groups_.end()) {
            for (auto ols : ctg_ols->second) {
                for (auto ol : ols.second) {
                    if (ol->a_.start < ol->a_.len*0.1 && (ol->a_.len - ol->a_.end) < ol->a_.len*0.1) {
                        cov[ol->b_.start]++;
                        cov[ol->b_.end]--;
                    }
                }
            }
            for (size_t i=1; i < cov.size(); ++i) {
                cov[i] += cov[i-1];
            }
            int start = -1;
            for (size_t i = 0; i + 1 < cov.size(); ++i) {
                if (cov[i] < 5) {
                    //printf("cov\t%s\t%zd\t%d\n", read_store_.QueryNameById(ctgid).c_str(), i, cov[i]);
                    if (start == -1) {
                        start = i;
                    }
                } else {
                    if (start != -1) {
                        printf("lowcov\t%s\t%d\t%zd\n", read_store_.QueryNameById(ctgid).c_str(), start, i);
                        start = -1;
                    }
                }
            }
        }
    }

    std::unordered_set<int> done;
    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        auto ol = ol_store_.Get(i);
        if (ol.a_.start < ol.a_.len*0.1 && (ol.a_.len - ol.a_.end) < ol.a_.len*0.1) {
            done.insert(ol.a_.id);
        }
    }

    for (size_t i = 0; i < read_store_.GetIdRange()[1]; ++i) {
        if (done.find(i) == done.end()) {
            printf("abn %s\n", read_store_.QueryNameById(i).c_str());
        }
    }
}


}   // namespace fsa