#include "contig_analyzer.hpp"

#include "align/match_info.hpp"

namespace fsa {

ContigAnalyzer::ContigAnalyzer(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds), win_slider_(cov_info_, WIN_SIZE, STRIDE), cov_info_(ds.seq_store_.GetSeq(tid)) {
}

void ContigAnalyzer::Detect() {

    ComputeCoverage();
    cov_info_.Scan();
    cov_info_.Stat();
 
    LOG(INFO)("ComputeCoverage");
    win_slider_.Flush();
    DetectErrors();

}

void ContigAnalyzer::ComputeCoverage() {
    // short name
    const ReadStore& seq_store = dataset_.seq_store_;
    auto  ol_group = dataset_.grouper_.Get(tid_);
    // 
    const DnaSeq& target = seq_store.GetSeq(tid_);

    for (size_t i = 0; i < ol_group.Size(); i++) {
        for (size_t j = 0; j < ol_group.Size(i); j++) {
            const auto &ol = *ol_group.Get(i, j);
            const DnaSeq& query = seq_store.GetSeq(ol.a_.id);
            assert(query.Size() >= 2000);
            if (ol.attached > 0) {
                match_.push_back(MatchInfo(&ol, query, target));
            } else {
                assert(ol.attached == 0);
            }
        }
    }

    std::vector<double> max_local_distances;
    for (const auto& m : match_) {
        max_local_distances.push_back(m.MaxLocalDistance());
    }

    auto mm = ComputeMedianAbsoluteDeviation(max_local_distances); // median, mad
    max_local_distance_threshold_ = mm[0] + 3*1.4826 * mm[1];
    LOG(INFO)("max_local_distance_threshold: %.02f = %.02f + 3*1.4826 * %.02f", max_local_distance_threshold_, mm[0], mm[1]);


    for (auto& m : match_) {
        if ( m.MatchedIdentity() >= dataset_.GetOverlapQualityThreshold() && (m.LClip() < MIN_CLIP || m.RClip() < MIN_CLIP )) {
            LOG(INFO)("add_seq %d: %s", m.GetOverlap()->a_.id, dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
            assert(m.GetOverlap()->attached > 0);
            cov_info_.Merge(m, MIN_CLIP*2, max_local_distance_threshold_, MIN_CLIP, 1.0 / m.GetOverlap()->attached);
            //cov_info_.Merge(m, MIN_CLIP*2, 0.50, MIN_CLIP, 1.0 / m.GetOverlap()->attached);
        }

    }
}


void ContigAnalyzer::SaveErrors(std::ofstream& of) {

    auto ctg_name = dataset_.QueryStringById(tid_);
    for (auto w : errors_) {
        of << ctg_name << " " << w.start << " " << w.end << " " << w.type << "\n";    
    }
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



bool ContigAnalyzer::CheckRegion(const ErrorRegion& reg) {
    double count = 0.0;
    auto s = reg.start / STRIDE;
    auto e = (reg.end - WIN_SIZE + 1) / STRIDE;

    for (auto& m : match_) {
        if ((m.Start() + 1000 < reg.start || m.Start() < 100) && 
            (m.End() > reg.end + 1000 || m.Len() - m.End() < 100) && 
            m.MatchedIdentity() >= dataset_.GetLocalQualityThreshold() &&
            m.MaxLocalDistance() < max_local_distance_threshold_ && m.LClip() < MIN_CLIP && m.RClip() < MIN_CLIP &&
            m.GetOverlap()->attached > 0 ) {
            
            count += 1.0 / m.GetOverlap()->attached;
            LOG(INFO)("checkreg: sup %zd %zd %s", m.Start(), m.End(), dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
        }
    }
    LOG(INFO)("checkreg: %s:%zd-%zd %.02f < %0.2f", Name().c_str(), reg.start, reg.end, count, this->win_slider_.SurroundingCoverage(reg) * 0.2);
    return count == 0 || count < this->win_slider_.SurroundingCoverage(reg) * 0.2;
}


void ContigAnalyzer::DetectErrors() {
    LOG(INFO)("DetectErrors");
    errors_ = win_slider_.DetectErrorRegions(1000, cov_info_.AvarageCoverage());
    LOG(INFO)("DetectErrors: %zd", errors_.size());
    errors_ = MergeRegions(errors_, 1000);
    LOG(INFO)("DetectErrors: merged %zd", errors_.size());
    
    errors_.erase(std::remove_if(errors_.begin(), errors_.end(), [this](const ErrorRegion& e) {
        bool torf = CheckRegion(e);
        if (torf) {
            LOG(INFO)("DetectErrors: add %s:%zd-%zd", Name().c_str(), e.start, e.end);

        } else {
            LOG(INFO)("DetectErrors: skip %s:%zd-%zd", Name().c_str(), e.start, e.end);
        }
        return !CheckRegion(e);
    }), errors_.end());
}


std::vector<ContigFragment> ContigAnalyzer::Split() {
    std::vector<ContigFragment> frgs;

    size_t idx = 0;
    for (size_t i = 0; i < errors_.size(); ++i) {
        const auto &e = errors_[i];
        if (e.start > idx) {
            frgs.push_back(ContigFragment(this, idx, e.start));
        }
        if (e.end > e.start) {
            frgs.push_back(ContigFragment(this, e.start, e.end));
        }
        idx = e.end;
    }
    if (idx < cov_info_.Size()) {
        frgs.push_back(ContigFragment(this, idx, cov_info_.Size()));
    }
    LOG(INFO)("Split: %zd fragments", frgs.size());
    return frgs;
}


void ContigAnalyzer::DumpCoverage(std::ofstream& of) {
    const auto &ctg_name = dataset_.QueryStringById(tid_);
    cov_info_.Dump(of, ctg_name);
}

void ContigAnalyzer::DumpWindow(std::ofstream& of) {
    const auto &ctg_name = dataset_.QueryStringById(tid_);
    win_slider_.Dump(of, ctg_name);

}


void ContigAnalyzer::DumpMatch(std::ofstream& of) {
    
    for (size_t i = 0; i < match_.size(); ++i) {
        const auto &m = match_[i];
        const auto ol = m.GetOverlap();
       const auto &rd_name = dataset_.QueryStringById(ol->a_.id);
       const auto &ctg_name = dataset_.QueryStringById(ol->b_.id);
        of << ctg_name << ":" << ol->b_.start << '-' << ol->b_.end << " " 
           << rd_name << " " << ol->a_.start << " " << ol->a_.end << " " << ol->a_.len << " " << ol->identity_ << "\n";
    }

}

}