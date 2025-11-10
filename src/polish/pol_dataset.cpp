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

    Stat();
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
    
    std::vector<std::array<double, 2>> quals;
    quals.reserve(rd_ids_[1] - rd_ids_[0]);

    struct Item {
        const Overlap* ol;
        double wt;
    };
    for (size_t ri = rd_ids_[0]; ri < rd_ids_[1]; ++ri) {
        auto gp = grouper_.Get(ri);
        std::vector<Item> items;
        for (size_t i = 0; i < gp.Size(); i++) {
            for (size_t j = 0; j < gp.Size(i); j++) {
                auto ol = gp.Get(i,j);
                items.push_back({ol, ol->Identity()*ol->AlignedLength() }); // TODO use a better weight function
                //items.push_back({ol, ol->Identity() }); // TODO use a better weight function
            }
        }
        if (items.size() > 0) {
            std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
                return a.wt > b.wt;
            });
            quals.push_back({items[0].ol->Identity(), items[0].ol->AlignedLength()});

            std::vector<int16_t> covs(items[0].ol->a_.len+1, 0);
            covs[items[0].ol->a_.start] ++;
            covs[items[0].ol->a_.end] --;
            size_t count = 1;
            for (; count < items.size(); ++count) {

                if (items[count].ol->AlignedLength()*2 > items[count].ol->a_.len) {
                    covs[items[count].ol->a_.start] ++;
                    covs[items[count].ol->a_.end] --;
                    break;
                } else {
                    break;
                }
            }

            for (size_t i = 1; i < covs.size(); ++i) {
                covs[i] += covs[i-1];
            }
            assert(covs.back() == 0);

            // if (count >= 2) {
            //     LOG(INFO)("data_sec: %s", QueryStringById(items[0].ol->a_.id).c_str());
            // }

            for (size_t i = 0; i < count; ++i) {
                double c = std::accumulate(covs.begin() + items[i].ol->a_.start, covs.begin() + items[i].ol->a_.end, 0) / 
                    (items[i].ol->a_.end - items[i].ol->a_.start);
                items[i].ol->attached = std::round(c);
            }
        }
    }

    std::sort(quals.begin(), quals.end(), [](const std::array<double,2> &a, const std::array<double,2> &b) {
        return a[0] > b[0];
    });

    double median = 0;
    double mad = 0;

    ComputeMedianAbsoluteDeviation(std::vector<std::array<double, 2>>(quals.begin(), quals.begin()+quals.size()*3/4), median, mad);
    overlap_quality_threshold_ = median-6*1.4826*mad;
    LOG(INFO)("Median = %.02f, MAD = %.02f, threshold=%.02f", median, mad, overlap_quality_threshold_);

}


    
size_t PolDataset::CountReadMap(Seq::Id id) const {
    size_t count = 0;
    auto gp = grouper_.Get(id);
    for (size_t i = 0; i < gp.Size(); ++i) {
        for (size_t j = 0; j < gp.Size(i); ++j) {
            const Overlap& ol = *gp.Get(i, j);
            if (ol.attached == 1) {
                count ++;
            }
        }
       
    }
    return count;
}

void PolDataset::Stat() {

    // calculate max read length
    size_t ave_length = 0;
    size_t count = 0;
    size_t max_read_length_ = 0;
    for (size_t i = rd_ids_[0]; i < rd_ids_[1]; ++i) {
        max_read_length_ = std::max(max_read_length_, seq_store_.GetSeqLength(i));
        ave_length += seq_store_.GetSeqLength(i);
        count++;
    }
    ave_read_length_ = ave_length / count;
    LOG(INFO)("Max read length: %zd, average read length: %zd", max_read_length_, ave_read_length_);
}

}   // namespace fsa