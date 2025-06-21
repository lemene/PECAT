#include "contig_error_detector.hpp"

#include "align/match_info.hpp"

namespace fsa {

ContigErrorDetector::ContigErrorDetector(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds), win_slider_(cov_info_, WIN_SIZE, STRIDE), cov_info_(ds.seq_store_.GetSeq(tid)) {
}

void ContigErrorDetector::Detect() {

    ComputeCoverage();
    cov_info_.Scan();
    //EvaluateQuality();
 
    LOG(INFO)("ComputeCoverage");
    win_slider_.Flush();
    LOG(INFO)("Start_Detect_Error_Regions %d", tid_);
    auto errs = win_slider_.DetectErrorRegions1(1000);
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
            if (ol.attached == 1) {
                match_.push_back(MatchInfo(&ol, query, target));
                //if ( match_.back().Identity() >= 0.85 && (match_.back().LClip() < MIN_CLIP || match_.back().RClip() < MIN_CLIP )) {
                if ( match_.back().MatchedIdentity() >= 0.90 && (match_.back().LClip() < MIN_CLIP || match_.back().RClip() < MIN_CLIP )) {
                    LOG(INFO)("add_seq %d: %s", match_.back().GetOverlap()->a_.id, dataset_.QueryStringById(match_.back().GetOverlap()->a_.id).c_str());
                    cov_info_.Merge(match_.back());
                }
                //break;
            }
            
            //break;  // 
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
    size_t count = 0;
    auto s = reg.start / STRIDE;
    auto e = (reg.end - WIN_SIZE + 1) / STRIDE;

    for (auto& m : match_) {
        if (m.Start() + 1000 < reg.start && m.End() > reg.end + 1000 && 
            m.MaxLocalDistance(1000) < 0.2 && m.LClip() < 1000 && m.RClip() < 1000) {
            count++;
            LOG(INFO)("checkreg: sup %zd %zd %s", m.Start(), m.End(), dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
        }
    }
    LOG(INFO)("checkreg: %s:%zd-%zd %zd", Name().c_str(), reg.start, reg.end, count);
    return count < 2;
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
        assert (c < 5);
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

// void ContigErrorDetector::EvaluateQuality() {
//     size_t ref_mis = 0;
//     size_t rd_mat = 0;
//     size_t rd_total = 0;
//     //for (const auto& c : ctg_cov_) {
//     for (size_t i = 0; i < ctg_cov_.size(); ++i) {
//         const auto& c = ctg_cov_[i];
//         auto m = std::max_element(c.bases, c.bases+6);
//         //LOG(INFO)("S(%zd): %s | %zd %zd", i, c.ToString().c_str(), m-c.bases, *m);

//         if ((m-c.bases) != c.ref) {
//             ref_mis ++;
//         }

//         rd_total += c.bases[0] + c.bases[1] + c.bases[2] + c.bases[3] + c.bases[5];

//         switch ((m-c.bases)) {
//         case 0:
//         case 1:
//         case 2:
//         case 3:
//             rd_mat += *m;
//             break;
//         case 4:
//             // 
//             break;
//         case 5:
//             rd_mat += c.inssize;
//             break;
//         default:
//             LOG(INFO)("XXX %zd", (m-c.bases));
//             assert(!"never coming here");
//         }
    
//     }
//     LOG(INFO)("%zd, %zd, %zd, %zd", ref_mis, ctg_cov_.size(), rd_mat, rd_total);
//     LOG(INFO)("%.02f, %.02f", ref_mis*1.0/ctg_cov_.size(), rd_mat*1.0 / rd_total);
// }

}