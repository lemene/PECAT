#include "prog_coverage_profile.hpp"

#include "utils/logger.hpp"
#include "../polish/contig_analyzer.hpp"

namespace fsa {

void Program_CoverageProfile::Running() {
    SET_LOG_LEVEL(DEBUG);
    dataset_.Load();
    AnalyzeContigs();
}


void Program_CoverageProfile::AnalyzeContigs() {

    for (size_t i = 0; i < dataset_.ctg_ids_.size(); i++) {
        auto tid = dataset_.ctg_ids_[i];    
        std::ofstream of_mcov(opts_.OutputPath(dataset_.QueryStringById(tid) + ".mcov"));
        Worker worker(tid, dataset_);
        worker.ComputeCoverage(opts_.thread_size);
        worker.DumpMultiCoverage(of_mcov);
    }
}


void Program_CoverageProfile::Worker::ComputeCoverage(size_t thread_size) {
    // short name
    const ReadStore& seq_store = dataset_.seq_store_;
    auto  ol_group = dataset_.grouper_.Get(tid_);
    // 
    const DnaSeq& target = seq_store.GetSeq(tid_);

    
    std::mutex mutex;
    auto combine_func = [&mutex, this](std::vector<std::vector<MultiCoverage::BaseCov>> &covs, std::vector<size_t>& pos) {
        std::lock_guard<std::mutex> lock(mutex);
        assert(covs.size() == pos.size());
        for (size_t i = 0; i < covs.size(); ++i) {
            multi_cov_.Merge(covs[i], pos[i]);
        }
        covs.clear();
        pos.clear();
    };
    
    std::atomic<size_t> index {0};
    auto work_func = [&](size_t _) {
        thread_local std::vector<std::vector<MultiCoverage::BaseCov>> local_covs;
        thread_local std::vector<size_t> local_pos;
        for (size_t i = index.fetch_add(1); i < ol_group.Size(); i = index.fetch_add(1)) {
            for (size_t j = 0; j < ol_group.Size(i); j++) {
                const auto &ol = *ol_group.Get(i, j);
                const DnaSeq& query = seq_store.GetSeq(ol.a_.id);
                if (query.Size() < 2000) continue;`
                // assert(query.Size() >= 2000);
                if (ol.attached > 0) {
                    auto mi = MatchInfo(&ol, query, target);
                    local_covs.push_back(multi_cov_.ToCov(mi, 1.0 / ol.attached));
                    local_pos.push_back(mi.Start());
                } else {
                    assert(ol.attached == 0);
                }
            }

            if (local_covs.size() > 1000) {
                combine_func(local_covs, local_pos);
            }
        }
        combine_func(local_covs, local_pos);
    };
    MultiThreadRun(thread_size, work_func);
}


void Program_CoverageProfile::Worker::DumpMultiCoverage(std::ofstream& of) {
    const auto &ctg_name = 
    of << ">" << dataset_.QueryStringById(tid_) << "\n";
    multi_cov_.Dump(of);
}


void Program_CoverageProfile::Dataset::Load() {
    
    seq_store_.Load(opts_.ctg_fname_, "", true);
    size_t id_start = seq_store_.GetIdUp();
    ctg_ids_ = LoadContigIds(seq_store_);

    seq_store_.Load(opts_.read_fname_, "", true);
    rd_ids_ = {id_start, seq_store_.GetIdUp()};

    rd_2_ctg_.LoadFast(opts_.mapping_fname_, "", (size_t)opts_.thread_size);
    LOG(INFO)("Load %zd overlaps from file %s", rd_2_ctg_.Size(), opts_.mapping_fname_.c_str());

    grouper_.BuildIndex(opts_.thread_size, std::unordered_set<int>());

    if (opts_.debug) {
        seq_store_.SaveIdToName("id_2_name");
    }
    SelectBestMapping();

    Stat();
}
    


std::vector<Seq::Id> Program_CoverageProfile::Dataset::LoadContigIds(const ReadStore& seq_store) {
    std::vector<Seq::Id> ctg_ids;
    std::array<size_t, 2> range = seq_store.GetIdRange();
    for (int i=range[0]; i<range[1]; ++i) {
        ctg_ids.push_back(i);
    }
    LOG(INFO)("Contig size: %zd", ctg_ids.size());
    return ctg_ids;
}

void Program_CoverageProfile::Dataset::SelectBestMapping() {
    
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

    ComputeMedianAbsoluteDeviation(std::vector<std::array<double, 2>>(quals.begin(), quals.begin()+quals.size()*3/4), overlap_quality_median_, overlap_quality_mad_);
    LOG(INFO)("Median = %.02f, MAD = %.02f", overlap_quality_median_, overlap_quality_mad_);

}

void Program_CoverageProfile::Dataset::Stat() {
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


} // namespace fsa
