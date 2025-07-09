#include "contig_error_detector.hpp"

#include "align/match_info.hpp"

namespace fsa {

ContigErrorDetector::ContigErrorDetector(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds), win_slider_(cov_info_, WIN_SIZE, STRIDE), cov_info_(ds.seq_store_.GetSeq(tid)) {
}

void ContigErrorDetector::Detect() {

    ComputeCoverage();
    cov_info_.Scan();
    cov_info_.Stat();
    //EvaluateQuality();
 
    LOG(INFO)("ComputeCoverage");
    win_slider_.Flush();
    LOG(INFO)("Start_Detect_Error_Regions %d", tid_);
    auto errs = win_slider_.DetectErrorRegions(1000, cov_info_.AvarageCoverage());
    LOG(INFO)("ComputeCoverage0");
    for (auto e : errs) {
        if (CheckRegion(e)) {
            errors_.push_back(e);
        }
    }
    LOG(INFO)("SplitSegments");
    SplitSegments();
}

void ContigErrorDetector::ComputeCoverage() {
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

    double median = 0.0;
    double mad = 0.0;
    ComputeMedianAbsoluteDeviation(max_local_distances, median, mad);
    LOG(INFO)("local distance median: %.2f, mad: %.2f", median, mad);
    max_local_distance_threshold_ = median + 3 * mad;
    LOG(INFO)("max_local_distance_threshold: %.2f", max_local_distance_threshold_);


    for (auto& m : match_) {
        if ( m.MatchedIdentity() >= dataset_.GetOverlapQualityThreshold() && (m.LClip() < MIN_CLIP || m.RClip() < MIN_CLIP )) {
            LOG(INFO)("add_seq %d: %s", m.GetOverlap()->a_.id, dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
            assert(m.GetOverlap()->attached > 0);
            //cov_info_.Merge(m, MIN_CLIP*2, max_local_distance_threshold_, MIN_CLIP, 1.0 / m.GetOverlap()->attached);
            cov_info_.Merge(m, MIN_CLIP*2, 0.50, MIN_CLIP, 1.0 / m.GetOverlap()->attached);
        }

    }
}


void ContigErrorDetector::SaveErrors(std::ofstream& of) {

    auto ctg_name = dataset_.QueryStringById(tid_);
    for (auto w : errors_) {
        of << ctg_name << " " << w.start << " " << w.end << " " << w.type << "\n";    
    }
}


std::vector<ErrorRegion> ContigErrorDetector::MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) {


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



bool ContigErrorDetector::CheckRegion(const ErrorRegion& reg) {
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
    LOG(INFO)("checkreg: %s:%zd-%zd %.02f %0.2f", Name().c_str(), reg.start, reg.end, count, this->win_slider_.SurroundingCoverage(reg) * 0.2);
    return count == 0 || count < this->win_slider_.SurroundingCoverage(reg) * 0.2;
    return count == 0 || count < std::min(this->win_slider_.SurroundingCoverage(reg) * 0.2, cov_info_.AvarageCoverage()[1]/2.0);
}


void ContigErrorDetector::Correct() {

}
std::string ContigErrorDetector::Consensus() {
    std::string seq;

    for (const auto& seg : segs_) {
        if (seg.type == 0) {
            seq += ConsensusSimple(seg);
        } else {
            seq += ConsensusComplex(seg);
        }
    }
    return seq;
}

std::string ContigErrorDetector::ConsensusSimple(const Segment& seg) {
    std::string seq;
    for (size_t i = seg.start; i < seg.end; ++i) {
        auto c = cov_info_.GetBestChoice(i);
        assert (c < 6);
        if (c <= 3) {
            seq += "ACGT"[c];
        }

    }
    return seq;
}

std::string ContigErrorDetector::ConsensusComplex(const Segment& seg) {
    std::string seq;
    for (size_t i = seg.start; i < seg.end; ++i) {
        auto c = cov_info_.GetBestChoice(i);
        if (c <= 3) {
            seq += "ACGT"[c];
        }

    }
    return seq;
}

void ContigErrorDetector::SplitSegments() {
    
    size_t i = 0;
    for (const auto &r : win_slider_.DetectSimpleRegions()) {

        if ( i < r.start) {
            segs_.push_back({i, r.start, 1});
            i = r.start;
        }
        segs_.push_back({r.start, r.end, 0});
        i = r.end;
    }

    if (i < cov_info_.Size()) {
        
        segs_.push_back({i, cov_info_.Size(), 1});
    }
    LOG(INFO)("segs: %zd", segs_.size());
    for (auto r : segs_) {
        LOG(INFO)("segs: %zd-%zd %d", r.start, r.end, r.type);
    }
}

void ContigErrorDetector::SaveContig(std::ofstream& of) {
    const auto &ctg_name = dataset_.QueryStringById(tid_);
    auto seq = Consensus();
    of << ">" << ctg_name << "\n" << seq << "\n";
}


void ContigErrorDetector::DumpCoverage(std::ofstream& of) {
    const auto &ctg_name = dataset_.QueryStringById(tid_);
    cov_info_.Dump(of, ctg_name);
}

void ContigErrorDetector::DumpWindow(std::ofstream& of) {
    const auto &ctg_name = dataset_.QueryStringById(tid_);
    win_slider_.Dump(of, ctg_name);

}


void ContigErrorDetector::DumpMatch(std::ofstream& of) {
    
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