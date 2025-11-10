#include "contig_analyzer.hpp"

#include "align/match_info.hpp"

namespace fsa {

ContigAnalyzer::ContigAnalyzer(Seq::Id tid, const PolDataset& ds)
 : tid_(tid), dataset_(ds), win_slider_(cov_info_, WIN_SIZE, STRIDE), cov_info_(ds.seq_store_.GetSeq(tid))
 , multi_cov_(ds.seq_store_.GetSeq(tid)) {
}

void ContigAnalyzer::Detect() {

    std::ofstream of(std::string("multi_cov_")+Name());
    ComputeCoverage();
    multi_cov_.Dump(of);
    // cov_info_.Scan();
    // cov_info_.Stat();
 
    // LOG(INFO)("ComputeCoverage");
    // win_slider_.Flush();
    // DetectErrors();

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
            assert(m.GetOverlap()->attached > 0);
            //cov_info_.Merge(m, MIN_CLIP*4, max_local_distance_threshold_, MIN_CLIP, 1.0 / m.GetOverlap()->attached);
            multi_cov_.Merge(m, 1.0 / m.GetOverlap()->attached);
            // TODO 参数化
        }

    }
    std::sort(match_.begin(), match_.end(), [](const MatchInfo& a, const MatchInfo& b) {
        return a.GetOverlap()->b_.start < b.GetOverlap()->b_.start || 
            (a.GetOverlap()->b_.start == b.GetOverlap()->b_.start && a.GetOverlap()->b_.end < b.GetOverlap()->b_.end);
    });
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
    return count == 0 ;//|| count < this->win_slider_.SurroundingCoverage(reg) * 0.2;
}


void ContigAnalyzer::DetectErrors() {
    LOG(INFO)("DetectErrors");
    errors_ = win_slider_.DetectErrorRegions2(1000, cov_info_.AvarageCoverage());
    LOG(INFO)("DetectErrors: %zd", errors_.size());
    errors_ = MergeRegions(errors_, 1000);
    LOG(INFO)("DetectErrors: merged %zd", errors_.size());
    
    errors_.erase(std::remove_if(errors_.begin(), errors_.end(), [this](const ErrorRegion& e) {
        // bool torf = CheckRegion(e);
        auto win = win_slider_.Region2Window(e);
        bool bpt = win_slider_.HasBreakpoint(win);
        bool alt = win_slider_.HasAlternate(win);
        LOG(INFO)("check_reg: %s:%zd-%zd, break=%d, alter=%d", Name().c_str(), e.start, e.end, bpt, alt);
        return !win_slider_.HasBreakpoint(win) && !win_slider_.HasAlternate(win);
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

size_t ContigAnalyzer::FirstMatch(size_t pos) {
    auto& match = match_;
    size_t left = 0, right = match.size();
    while (left < right) {
        size_t mid = (left + right) / 2;
        LOG(INFO)("FirstMatch: %zd %zd %zd | %zd %zd", left, right, mid, pos, match[mid].Start());
        if (match[mid].Start() <= pos && match[mid].End() >= pos) {
            right = mid;
        } else if (match[mid].Start() > pos) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    return left;
}


size_t ContigAnalyzer::LastMatch(size_t pos) {
    auto& match = match_;
    size_t left = 0, right = match.size();
    while (left < right) {
        size_t mid = (left + right) / 2;
        LOG(INFO)("LastMatch: %zd %zd %zd | %zd %zd", left, right, mid, pos, match[mid].Start());
        if (match[mid].Start() <= pos) {
            left = mid + 1;
        } else {
            //if (match[mid].Start() > pos) 
            right = mid;
        } 
    }
    return left;
}

std::string ContigAnalyzer::Polish(size_t s, size_t e) {
    assert(s >= 0 && e >= s && cov_info_.Size() >= e);

    std::vector<uint8_t> table(e-s);
    for (size_t i = s; i < e; ++i) {
        table[i-s] = cov_info_.Status(i);
    }

    std::vector<std::array<size_t, 2>> regions;
    size_t start = 0;
    for (size_t i = 1; i < table.size(); ++i) {
        if (table[i] != table[start]) {
            regions.push_back({start, i});
            start = i;
        }
    }
    if (start < table.size()) {
        regions.push_back({start, table.size()});
    } 
    
    // LOG(INFO)("Polish: %zd regions", regions.size());
    std::vector<std::array<size_t, 2>> merged;
    merged.reserve(regions.size());
    const int FLANKING = 10;
    merged.push_back(regions[0]);
    for (size_t i = 1; i < regions.size(); ++i) {
        if (table[regions[i][0]] == table[merged.back()[0]]) {
            merged.back()[1] = regions[i][1];
        } else if (table[regions[i][0]] == 0 && regions[i][1] - regions[i][0] < 2*FLANKING) {
            merged.back()[1] = regions[i][1];
        } else {
            merged.push_back(regions[i]);
        }
    }
    std::string seq;
    for (size_t i = 0; i < regions.size(); ++i) {
        const auto& r = regions[i];
        if (true || table[regions[i][0]] == 0) {
            for (size_t j = regions[i][0]; j < regions[i][1]; ++j) {
  
                auto c = cov_info_.GetBestChoice(j + s);
                if (c < 4) {
                    seq += "ACGT"[c];
                } else if (c == 4) {
                    // deletion, do nothing
                } else if (c == 5) {
                    // insertion, TODO
                } else {
                    LOG(ERROR)("Polish: unexpected base %d at %zd", c, i);
                }
            }

        } else {
            assert(table[r[0]] == 1);
            

            size_t ss = r[0] + s > FLANKING ? r[0] + s - FLANKING : 0;
            size_t ee = r[1] + s + FLANKING < cov_info_.Size() ? r[1] + s + FLANKING : cov_info_.Size();

            std::vector<DnaSeq> segs;
            LOG(INFO)("Polish: error region %zd-%zd, find_position %zd %zd", ss, ee, FirstMatch(ss), match_.size());
            for (size_t im = FirstMatch(ss); im < match_.size() && match_[im].Start() <= ss; ++im) {
                const auto& m = match_[im];
                if (m.MatchedIdentity() >= dataset_.GetLocalQualityThreshold() &&
                    m.MaxLocalDistance() < max_local_distance_threshold_ && m.LClip() < MIN_CLIP && m.RClip() < MIN_CLIP &&
                    m.GetOverlap()->attached > 0 && m.Start() <= ss && m.End() >= ee) {

                    auto qr = m.GetQueryRegion(ss, ee);
                    LOG(INFO)("Query region: %zd-%zd %zd-%zd %s", ss, ee, qr[0], qr[1], dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
                    
                    auto seg = (qr[0] < qr[1]) ?  
                        DnaSeq(dataset_.seq_store_.GetSeq(m.GetOverlap()->a_.id), qr[0], qr[1] - qr[0]) :
                        DnaSeq::ReverseComplement(DnaSeq(dataset_.seq_store_.GetSeq(m.GetOverlap()->a_.id), qr[1], qr[0] - qr[1]));
                    segs.push_back(seg);
                    
                    
                    LOG(INFO)("Polish: seg %d %s", qr[0] < qr[1], seg.ToString()->c_str());
                }
            }
            // for (auto& seg : segs) {
            //     LOG(INFO)("Polish: seg %s", seg.ToString()->c_str());
            // }
            if (segs.size() > 0) {
                seq += segs[0].ToString()->c_str();
            }
        }
    }

    return seq;
}

std::vector<const MatchInfo*> ContigAnalyzer::GetCoverage(size_t pos, int flank) {
    std::vector<const MatchInfo*> cov_matches;

    if (flank < 0) {
        size_t s_pos = pos < -flank ? 0 : pos + flank;
        size_t e_pos = pos ;
        size_t ifirst = FirstMatch(s_pos + dataset_.MaxReadLength());
        size_t ilast = LastMatch(s_pos);
        for (size_t im = ifirst; im < ilast; ++im) {
            const auto& m = match_[im];
            LOG(INFO)("GetCoverage0(right): %zd %zd %zd %s", pos, m.Start(), m.End(), dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
            if ( m.LClip() < MIN_CLIP && m.GetOverlap()->attached > 0 && m.Start() <= s_pos && m.End() >= e_pos) {
                
                LOG(INFO)("GetCoverage1(right): %zd %zd %zd %s", pos, m.Start(), m.End(), dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
                cov_matches.push_back(&m);
                
            }
        }
    } else {
        size_t s_pos = pos ;
        size_t e_pos = pos + flank >= cov_info_.Size() ? cov_info_.Size() : pos + flank;
        size_t ifirst = FirstMatch(s_pos + dataset_.MaxReadLength());
        size_t ilast = LastMatch(s_pos);
        for (size_t im = ifirst; im < ilast; ++im) {
            const auto& m = match_[im];
            LOG(INFO)("GetCoverage0(left): %zd %zd %zd %s", pos, m.Start(), m.End(), dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
            if ( m.RClip() < MIN_CLIP && m.GetOverlap()->attached > 0 && m.Start() <= s_pos && m.End() >= e_pos) {
                
                LOG(INFO)("GetCoverage1(left): %zd %zd %zd %s", pos, m.Start(), m.End(), dataset_.QueryStringById(m.GetOverlap()->a_.id).c_str());
                cov_matches.push_back(&m);
                
            }
        }
    }
    return cov_matches;
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

void ContigAnalyzer::DumpMultiCoverage(std::ofstream& of) {
    const auto &ctg_name = 
    of << ">" << dataset_.QueryStringById(tid_) << "\n";
    multi_cov_.Dump(of);
}

}