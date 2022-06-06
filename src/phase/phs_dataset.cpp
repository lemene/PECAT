
#include "phs_dataset.hpp"

#include "../utils/logger.hpp"

#include "phs_options.hpp"


namespace fsa {

void PhsDataset::Load() {

    LOG(INFO)("Preload reads and contigs");
    PreloadReads();

    LOG(INFO)("Load overlaps");
    LoadOverlaps(opts_.ol_fname_);

    if (!opts_.rd2rd.empty()) {
        LoadAva(opts_.rd2rd);
    }

    LOG(INFO)("Load reads and contigs");
    LoadReads();
    
}   


void PhsDataset::LoadOverlaps(const std::string &fname) {
    auto filter_simple = [&](Overlap& o) {
        if (!opts_.filter_opts_.Valid(o)) {
            return false;
        }

        if (contig_ids_.find(o.b_.id) == contig_ids_.end()) {
            return false;
        }

//        // TODO for debug
//        if (o.b_.start < 66229153 - 1000000 || o.b_.start > 66332061 + 1000000) return false;

        return true;
    };

    ol_store_.LoadFast(fname, "", (size_t)opts_.thread_size_, filter_simple);

    if (ol_store_.Size() > 0) {
        LOG(INFO)("Overlap size: %zd", ol_store_.Size());
    } else {
        LOG(ERROR)("No overlap was loaded");
    }

    //ol_store_.Save(opts_.OutputPath("load.paf"));
}


void PhsDataset::LoadAva(const std::string &fname) {
 

    ol_ava_.LoadFast(fname, "", (size_t)opts_.thread_size_);

    if (ol_ava_.Size() > 0) {
        LOG(INFO)("Overlap size: %zd", ol_ava_.Size());
    } else {
        LOG(ERROR)("No overlap was loaded");
    }
    ol_ava_.Group(ava_groups_, opts_.thread_size_);
}

std::unordered_set<int> PhsDataset::QueryGroup(int id) {
    auto iter = ava_groups_.find(id);
    if (iter != ava_groups_.end()) {
        std::unordered_set<int> sss;
        for (auto ii : iter->second) {
            sss.insert(ii.first);
        }
        return sss;
    } else {
        return {};
    }
}

void PhsDataset::LoadReads() {
    std::unordered_set<Seq::Id> ids;

    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        ids.insert(ol_store_.Get(i).a_.id);
        read_ids2_[ol_store_.Get(i).a_.id]++;
    }

    rd_store_.Load(opts_.rd_fname_, "", true, ids);
    rd_store_.Load(opts_.ctg_fname_, "", true, contig_ids_);
}

void PhsDataset::PreloadReads() {
    rd_store_.Load(opts_.rd_fname_, "", false, std::unordered_set<Seq::Id>());
    rd_store_.Load(opts_.ctg_fname_, "", false, std::unordered_set<Seq::Id>());
    
    if (!opts_.ctgname_fname_.empty()) {
        GzFileReader reader(opts_.ctgname_fname_);
        if (reader.Valid()) {
            std::string name = reader.GetNoEmptyLine();
            while (!name.empty()) {
                contig_names_.insert(name);
                contig_ids_.insert(rd_store_.QueryIdByName(name));
                name = reader.GetNoEmptyLine();
            }
        } else {
            LOG(ERROR)("Failed to load file: %s", opts_.ctgname_fname_.c_str());
        }
    } else {
        contig_ids_ = rd_store_.IdsInFile(opts_.ctg_fname_);
        for (auto id : contig_ids_) {
            contig_names_.insert(rd_store_.QueryNameById(id));
        }
    }

    if (opts_.loglevel >= 4) rd_store_.SaveIdToName(opts_.OutputPath("id2name.gz"));

}

std::vector<Seq::Id> PhsDataset::GetSortedContigs() const {
    std::vector<Seq::Id> sorted(contig_ids_.begin(), contig_ids_.end());

    std::unordered_map<Seq::Id, size_t> counts;

    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);
        counts[o.b_.id] += 1;
    }

    for (auto &c : counts) {
        c.second = c.second * c.second  / std::max<size_t>(1, (rd_store_.GetSeqLength(c.first) / 1000000)); 
    }

    std::sort(sorted.begin(), sorted.end(), [&](Seq::Id a, Seq::Id b) {
        return counts[a] > counts[b];
    });

    return sorted;
}

}