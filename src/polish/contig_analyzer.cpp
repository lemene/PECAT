#include "contig_analyzer.hpp"

#include "align/match_info.hpp"

namespace fsa {

ContigAnalyzer::ContigAnalyzer(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds), multi_cov_(ds.seq_store_.GetSeq(tid), ds.overlap_quality_median_, ds.overlap_quality_mad_) {
}



void ContigAnalyzer::ComputeCoverage(size_t thread_size) {
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
                assert(query.Size() >= 2000);
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


std::vector<ErrorRegion> ContigAnalyzer::MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) {

    std::vector<ErrorRegion> merged;
    if (regs.size() > 0) {
        merged.push_back(regs[0]);

        for (size_t i = 1; i < regs.size(); ++i) {
            if (regs[i].start <= merged.back().end + max_gap) {
                assert(merged.back().end <= regs[i].end);
                merged.back().end = regs[i].end;
            } else {
                merged.push_back(regs[i]);
            }
        }
    }
    return merged;
}


void ContigAnalyzer::DumpMultiCoverage(std::ofstream& of) {
    const auto &ctg_name = 
    of << ">" << dataset_.QueryStringById(tid_) << "\n";
    multi_cov_.Dump(of);
}

}