#include "window_slider.hpp"

namespace fsa {


void WindowSlider::Flush() {
    std::vector<ErrorRegion> cands;

    WinDivider witr = WinDivider(cov_info_.Size(), win_size_, stride_);
    win_cov_.assign(witr.Size(), WinInfo());
    
    for (size_t i = 0; i < witr.Size(); ++i) {
        auto win = witr.Get(i);
        auto &wc = win_cov_[i];

        size_t clips = 0;
        std::vector<size_t> cov0(win[1] - win[0], 0);
        std::vector<size_t> cov1(win[1] - win[0], 0);
        std::vector<size_t> cov2(win[1] - win[0], 0);

        std::vector<int> ss(win[1] - win[0], 0);

        for (size_t ii = win[0]; ii < win[1]; ++ii) {
            const auto& c = cov_info_.Get(ii);

            cov0[ii-win[0]] = std::accumulate(c.bases0, c.bases0+5, 0);
            cov1[ii-win[0]] = std::accumulate(c.bases1, c.bases1+5, 0);
            cov2[ii-win[0]] = std::accumulate(c.bases2, c.bases2+5, 0);
            
            wc.clips += c.clips;
            wc.ins += c.bases0[5];
            win_cov_[i].count[0] += c.c[0];
            win_cov_[i].count[1] += c.c[1];
            win_cov_[i].count[2] += c.c[2];
        }
        win_cov_[i].c0 = std::accumulate(cov0.begin(), cov0.end(), 0) * 1.0 / cov0.size();
        win_cov_[i].c1 = std::accumulate(cov1.begin(), cov1.end(), 0) * 1.0 / cov0.size();
        win_cov_[i].c2 = std::accumulate(cov2.begin(), cov2.end(), 0) * 1.0 / cov0.size();
        win_cov_[i].min_c[0] = *std::min_element(cov0.begin(), cov0.end());
        win_cov_[i].max_c[0] = *std::max_element(cov0.begin(), cov0.end());
        win_cov_[i].min_c[1] = *std::min_element(cov1.begin(), cov1.end());
        win_cov_[i].max_c[1] = *std::max_element(cov1.begin(), cov1.end());
        win_cov_[i].min_c[2] = *std::min_element(cov2.begin(), cov2.end());
        win_cov_[i].max_c[2] = *std::max_element(cov2.begin(), cov2.end());


        win_cov_[i].count[0] /= win[1] - win[0];
        win_cov_[i].count[1] /= win[1] - win[0];
        win_cov_[i].count[2] /= win[1] - win[0];
        LOG(INFO)("win %zd-%zd %s", win[0], win[1], win_cov_[i].ToString().c_str());
    }
}

auto WindowSlider::DetectErrorRegions(size_t max_gap, const std::array<double,3>& ave_covs) -> std::vector<ErrorRegion>  {

    ComputeCoverageThresholds(0);
    ComputeCoverageThresholds(1);
    ComputeCoverageThresholds(2);
    // Detect error regions based on coverage information
    // - Detect regions with low coverage, high coverage, or significant differences between coverage types
    // 
    const size_t MIN_COV_HQ = 3;
    const double MAX_COV_DIFF = 0.2; 
    const double MIN_COV_DIFF = 0.1; 

    const size_t CAND_WIN_SIZE = 10000;
    const size_t CAND_WIN_COUNT = CAND_WIN_SIZE / stride_;

    auto cand_ave_covs = AveCoverages(CAND_WIN_COUNT);
    auto cand_max_covs = MaxCoverages(CAND_WIN_COUNT);
    auto cand_min_covs = MinCoverages(CAND_WIN_COUNT);



    std::vector<ErrorRegion> cands;
    
    std::vector<double> l_slope;
    std::vector<double> l_slope_rate;
    std::vector<int> kernel = {1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1};
    const int KLEN = kernel.size();
    assert(KLEN % 2 == 1 && KLEN == kernel.size());

    for (size_t i = KLEN/2; i + KLEN/2 < win_cov_.size(); ++i) {
        int s = 0;
        int ss = 0;
        for (size_t ii = 0; ii < KLEN; ++ii) {
            s += int(win_cov_[i + ii - KLEN/2].c1)*kernel[ii];
            ss += int(win_cov_[i + ii - KLEN/2].c1);
        }
        l_slope.push_back(s);
        l_slope_rate.push_back(ss == 0 ? 0.0 : s * 1.0 / ss);
        LOG(DEBUG)("detect_slope_item %.02f %.02f", l_slope.back(), l_slope_rate.back());
    }

    auto mm = ComputeMedianAbsoluteDeviation(l_slope);
    auto low = mm[0] - 3*1.4826*mm[1];
    auto high = mm[0] + 3*1.4826*mm[1];
    LOG(INFO)("PEAK threashold: %.02f %.02f -> %.02f %.02f", mm[0], mm[1], low, high);

    mm = ComputeMedianAbsoluteDeviation(l_slope_rate);
    auto low_rate = mm[0] - 3*1.4826*mm[1];
    auto high_rate = mm[0] + 3*1.4826*mm[1];
    LOG(INFO)("PEAK rate threashold %.02f %.02f -> %.02f %.02f", mm[0], mm[1], low_rate, high_rate);

    for (size_t i = 0 ; i < win_cov_.size(); ++i) {
        auto &winfo = win_cov_[i];
        std::array<size_t,2> win = Window2Region(i);
        assert(winfo.type == 0); // Not yet set type

        if (winfo.min_c[1] * 1.0 < ave_covs[1] * 0.5 || winfo.max_c[1] * 1.0  > ave_covs[1] * 1.5) {
            winfo.type = 1;
        }

        if (win_cov_[i].c1 <= MIN_COV_HQ) {
            winfo.type = 2;
        }

        if ((win_cov_[i].c0- win_cov_[i].c1)*1.0 > std::max<double>(MIN_COV_HQ, win_cov_[i].c1*MAX_COV_DIFF)) {
            winfo.type = 2; 
        } else if ((win_cov_[i].c0- win_cov_[i].c1)*1.0 > std::max<double>(MIN_COV_HQ, win_cov_[i].c1*MIN_COV_DIFF)) {
            winfo.type = 3;
        }

        if (i >= KLEN / 2 && i + KLEN / 2 < win_cov_.size()) {
            // Check slope
            auto sp = l_slope[i-KLEN/2];
            auto spr = l_slope_rate[i-KLEN/2];  
            LOG(INFO)("detect_slope %zd-%zd %.02f < %.02f < %.02f |  %.02f < %.02f < %.02f ", win[0], win[1], low, sp,  high, low_rate, spr, high_rate);
            if (sp < low || sp > high || spr < low_rate || spr > high_rate) {
                winfo.type = 2;
            }
        }
        
        LOG(INFO)("detect %zd-%zd %d %.02f", win[0], win[1], winfo.type);

    }
    
    for (size_t i = 0; i < win_cov_.size(); ++i) {
        std::array<size_t,2> win = Window2Region(i);
        auto &winfo = win_cov_[i];
        LOG(INFO)("type:%zd-%zd %d", win[0], win[1], winfo.type);
        if (winfo.type == 2 || winfo.type == 4) {
            cands.push_back({win[0], win[1], 0});

            for (size_t ii = i; ii > 0; --ii) {
                auto &winfo = win_cov_[ii-1]; 
                std::array<size_t,2> win = Window2Region(ii-1);
                if (winfo.type != 0) {
                    cands.push_back({win[0], win[1], 0});
                } else {
                    break;
                }
            }
            for (size_t ii = i+1; ii < win_cov_.size(); ++ii) {
                auto &winfo = win_cov_[ii];
                std::array<size_t,2> win = Window2Region(ii);
                if (winfo.type != 0) {
                    cands.push_back({win[0], win[1], 0});
                } else {
                    break;
                }
            }
        }
    }


    std::sort(cands.begin(), cands.end(), [](ErrorRegion& a, ErrorRegion &b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    return MergeRegions2(MergeRegions(cands, max_gap));
}


auto WindowSlider::DetectErrorRegions2(size_t max_gap, const std::array<double,3>& ave_covs) -> std::vector<ErrorRegion>  {

    const size_t MIN_COV_HQ = 3;
    const double MAX_COV_DIFF = 0.2; 
    const double MIN_COV_DIFF = 0.1; 
    
    const size_t CAND_WIN_SIZE = 10000;
    const size_t CAND_WIN_COUNT = CAND_WIN_SIZE / stride_;

    auto cand_ave_covs = AveCoverages(CAND_WIN_COUNT);
    auto cand_std_covs = StdCoverages(CAND_WIN_COUNT, cand_ave_covs);
    auto cand_max_covs = MaxCoverages(CAND_WIN_COUNT);
    auto cand_min_covs = MinCoverages(CAND_WIN_COUNT);

    std::vector<double> cand_diff_covs(cand_max_covs.size());
    std::vector<double> cand_slope_covs(cand_max_covs.size());
    std::vector<double> cand_rdiff_covs(cand_max_covs.size());
    std::vector<double> cand_rslope_covs(cand_max_covs.size());
    for (size_t i = 0; i < cand_diff_covs.size(); i++) {

        cand_diff_covs[i] = cand_max_covs[i].first - cand_min_covs[i].first;
        cand_slope_covs[i] = cand_diff_covs[i] / ((int)cand_max_covs[i].second - (int)cand_min_covs[i].second);
        
        cand_rdiff_covs[i] = (cand_max_covs[i].first - cand_min_covs[i].first) / cand_ave_covs[i];
        cand_rslope_covs[i] = cand_rdiff_covs[i] / ((int)cand_max_covs[i].second - (int)cand_min_covs[i].second);
        
        LOG(INFO)("dt_stat(%zd): %.02f-%zd %.02f-%zd, %.02f %.02f, %.02f %.02f, %.02f, %.02f", 
            i, cand_max_covs[i].first, cand_max_covs[i].second, cand_min_covs[i].first, cand_min_covs[i].second,
            cand_diff_covs[i], cand_slope_covs[i],
            cand_rdiff_covs[i], cand_rslope_covs[i],
            cand_ave_covs[i], cand_std_covs[i]
        );
    }

    auto mm_ave = ComputeMedianAbsoluteDeviation(cand_ave_covs);
    auto low_ave = mm_ave[0] - 3*1.4826*mm_ave[1];
    auto high_ave = mm_ave[0] + 3*1.4826*mm_ave[1];
    LOG(INFO)("dt_th_ave: %.02f %.02f -> %.02f %.02f", mm_ave[0], mm_ave[1], low_ave, high_ave);

    auto mm_std = ComputeMedianAbsoluteDeviation(cand_std_covs);
    auto low_std = mm_std[0] - 3*1.4826*mm_std[1];
    auto high_std = mm_std[0] + 3*1.4826*mm_std[1];
    LOG(INFO)("dt_th_std: %.02f %.02f -> %.02f %.02f", mm_std[0], mm_std[1], low_std, high_std);

    auto mm_diff = ComputeMedianAbsoluteDeviation(cand_diff_covs);
    auto low_diff = mm_diff[0] - 3*1.4826*mm_diff[1];
    auto high_diff = mm_diff[0] + 3*1.4826*mm_diff[1];
    LOG(INFO)("dt_th_diff: %.02f %.02f -> %.02f %.02f", mm_diff[0], mm_diff[1], low_diff, high_diff);

    auto mm_slope = ComputeMedianAbsoluteDeviation(cand_slope_covs);
    auto low_slope = mm_slope[0] - 3*1.4826*mm_slope[1];
    auto high_slope = mm_slope[0] + 3*1.4826*mm_slope[1];
    LOG(INFO)("dt_th_slope: %.02f %.02f -> %.02f %.02f", mm_slope[0], mm_slope[1], low_slope, high_slope);

    auto mm_rdiff = ComputeMedianAbsoluteDeviation(cand_rdiff_covs);
    auto low_rdiff = mm_rdiff[0] - 3*1.4826*mm_rdiff[1];
    auto high_rdiff = mm_rdiff[0] + 3*1.4826*mm_rdiff[1];
    LOG(INFO)("dt_th_rdiff: %.02f %.02f -> %.02f %.02f", mm_rdiff[0], mm_rdiff[1], low_rdiff, high_rdiff);

    auto mm_rslope = ComputeMedianAbsoluteDeviation(cand_rslope_covs);
    auto low_rslope = mm_rslope[0] - 3*1.4826*mm_rslope[1];
    auto high_rslope = mm_rslope[0] + 3*1.4826*mm_rslope[1];
    LOG(INFO)("dt_th_rslope: %.02f %.02f -> %.02f %.02f", mm_rslope[0], mm_rslope[1], low_rslope, high_rslope);

    std::vector<ErrorRegion> cands;
    for (size_t i = 0 ; i < win_cov_.size(); ++i) {
        auto &winfo = win_cov_[i];
        std::array<size_t,2> win = Window2Region(i);
        assert(winfo.type == 0); // Not yet set type

        if (win_cov_[i].c1 == 0) {
            winfo.type = 2;
        }
        
        LOG(INFO)("detect %zd-%zd %d %.02f", win[0], win[1], winfo.type, win_cov_[i].c1);
    }

    // for (size_t i = 0; i < cand_std_covs.size(); ++i) {
    //     if (cand_std_covs[i] >= high_std) {
    //         LOG(INFO)("detect_std: %zd, %0.2f > %.02f", i, cand_std_covs[i], high_std);
    //         if (cand_slope_covs[i] <low_slope || cand_slope_covs[i] > high_slope) {
    //             auto s = std::min(cand_max_covs[i].second, cand_min_covs[i].second);
    //             auto e = std::max(cand_max_covs[i].second, cand_min_covs[i].second);
    //             LOG(INFO)("DIFF: %zd %zd %zd", i, s, e);
    //             for (size_t ii = s; ii <= e; ++ii) {
    //                 auto &winfo = win_cov_[ii];
    //                 winfo.type = 4;
    //             }
    //         }
    //     }
    // }
    
    for (size_t i = 0; i < cand_rdiff_covs.size(); ++i) {
        if (cand_rdiff_covs[i] > mm_rdiff[0] + mm_rdiff[1] || cand_diff_covs[i] > mm_diff[0] + mm_diff[1] ) {
            auto s = std::min(cand_max_covs[i].second, cand_min_covs[i].second);
            auto e = std::max(cand_max_covs[i].second, cand_min_covs[i].second);
            LOG(INFO)("DIFF: %zd %zd %zd", i, s, e);
            for (size_t ii = s; ii <= e; ++ii) {
                auto &winfo = win_cov_[ii];
                winfo.type = 4;
            }
        }
    }
    
    for (size_t i = 0; i < win_cov_.size(); ++i) {
        std::array<size_t,2> win = Window2Region(i);
        auto &winfo = win_cov_[i];
        LOG(INFO)("type:%zd-%zd %d", win[0], win[1], winfo.type);
        if (winfo.type == 2 || winfo.type == 4) {
            cands.push_back({win[0], win[1], 0});
        }
    }

    std::sort(cands.begin(), cands.end(), [](ErrorRegion& a, ErrorRegion &b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    return MergeRegions(cands, max_gap);
}

std::vector<ErrorRegion> WindowSlider::MergeRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) const {


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

std::vector<ErrorRegion> WindowSlider::MergeRegions2(const std::vector<ErrorRegion> &regs) const {

    // merge
    std::vector<int8_t> slopes(regs.size(), 0);
    auto calc_slope = [this](const ErrorRegion& reg) {
        auto win = Region2Window(reg);
        win[0] -= 1;
        win[1] += 1;
        LOG(INFO)("www :%d %d %d %d", win[0], win[1], reg.start, reg.end);
        if (win_cov_[win[0]].c0 + win_cov_[win[0]].c0*0.2 < win_cov_[win[1]].c0) {
            return -1;
        } else if (win_cov_[win[0]].c0 > win_cov_[win[1]].c0 + win_cov_[win[1]].c0*0.2) {
            return 1;
        } else {
            return 0;
        }
    };

    std::transform(regs.begin(), regs.end(), slopes.begin(), calc_slope);
    std::vector<ErrorRegion> merged;
    size_t idx = 0;
    while (idx < regs.size()) {
        merged.push_back(regs[idx]);
        LOG(INFO)("mm0 %zd-%zd %d", regs[idx].start, regs[idx].end, slopes[idx]);

        if (slopes[idx] != 0) {
            for (size_t i = idx + 1; i < regs.size(); ++i) {
                if (regs[i].end - merged.back().start >= 1000000) {
                    break;
                }

                if (slopes[i] == 0) continue;
                if (slopes[i] * slopes[idx] == -1) {
                    LOG(INFO)("mmm %zd-%zd %d, %zd-%zd %d", regs[idx].start, regs[idx].end, slopes[idx], regs[i].start, regs[i].end, slopes[i]);
                    bool is_whole = true;
                    auto w0 = Region2Window(regs[idx]);
                    auto w1 = Region2Window(regs[i]);

                    for (size_t ii = w0[1] + 1; ii < w1[0]; ++ii) {
                        LOG(INFO)("mmm-  %.02f %02f %02f", win_cov_[ii].c0, win_cov_[w0[0]].c0, win_cov_[w1[1]].c0);
                        if (slopes[idx] == -1 && win_cov_[ii].c0 < (win_cov_[w0[0]].c0 + win_cov_[w1[1]].c0 )/ 2) {
                            is_whole = false;
                            break;
                        } else if (slopes[idx] == 1 && win_cov_[ii].c0 > (win_cov_[w0[0]].c0 + win_cov_[w1[1]].c0 )/ 2) {
                            is_whole = false;
                            break;
                        }
                    }
                    if (is_whole) {
                        merged.back().end = regs[i].end;
                        idx = i ;
                    } 
                }
                break;
            }
        }
        idx += 1;
    }
    return merged;
}


void WindowSlider::Dump(std::ofstream &of, const std::string& ctg_name) {
    for (size_t i = 0; i < win_cov_.size(); ++i) {
        const auto &w = win_cov_[i];
        size_t s = stride_ * i;
        size_t e = std::min<size_t>(s + win_size_, cov_info_.Size());
        of << ctg_name << ":" << s << '-' << e << " " << w.ToString() << "\n";
    }
}

double WindowSlider::SurroundingCoverage(const ErrorRegion& reg, size_t inv) {
 
    auto win = Region2Window(reg);
 
    double laccu = 0;
    size_t lsz = 0;
    for (size_t i = win[0]; i > 0 && win[0] - i + 1 < inv; --i) {
        laccu = std::max(win_cov_[i-1].c1, laccu);
        lsz += 1;
    }

    double raccu = 0;
    size_t rsz = 0;
    for (size_t i = win[1]+1; i < win_cov_.size() && i - win[1] < inv; ++i) {
        raccu = std::max(win_cov_[i].c1, raccu);
        rsz += 1;
    }
    LOG(INFO)("Sur:%zd-%zd %zd-%zd %.02f %zd %.02f %zd",reg.start, reg.end, win[0], win[1], laccu, lsz, raccu, rsz);
    return std::min((lsz == 0 ? 0.0 : laccu), (rsz == 0 ? 0.0 : raccu));
}

std::array<double,2> WindowSlider::ComputeCoverageThresholds(int type) {
    std::vector<double> covs;
    for (size_t i = 0; i < win_cov_.size(); ++i) {
        covs.push_back(type == 0 ? win_cov_[i].c0 : (type == 1 ? win_cov_[i].c1 : win_cov_[i].c2));
    }
    auto mm = ComputeMedianAbsoluteDeviation(covs);
    auto low  = mm[0] - 3*1.4826*mm[1];
    auto high = mm[0] + 3*1.4826*mm[1];
    LOG(INFO)("Coverage(%d) median mad low high %.02f %.02f %.02f %.02f", type, mm[0], mm[1], low, high);
    return {low, high};
}


std::vector<double> WindowSlider::AveCoverages(size_t winnum) {
    std::vector<double> aves;

    if (win_cov_.size() >= winnum) {
        double ave = std::accumulate(win_cov_.begin(), win_cov_.begin() + winnum, 0.0, [](double sum, const WinInfo& w) {
            return sum + w.c1;
        }) / winnum;
        aves.push_back(ave);

        for (size_t i = winnum; i < win_cov_.size(); ++i) {
            ave += win_cov_[i].c1 / winnum;
            ave -= win_cov_[i - winnum].c1 / winnum;
            aves.push_back(ave);
        }
    }
    return aves;

}

std::vector<double> WindowSlider::StdCoverages(size_t winnum, const std::vector<double>& ave) {
    std::vector<double> std(ave.size());

    for (size_t i = 0; i < ave.size(); ++i) {
        for (size_t iw = 0; iw < winnum; ++iw) {
            std[i] += std::abs(win_cov_[iw+i].c1 - ave[i]) / winnum;
        }
    }
    return std;
}

std::vector<std::pair<double,size_t>> WindowSlider::MaxCoverages(size_t winnum) {
    std::vector<std::pair<double,size_t>> max_covs;

    auto find_max_posion = [](const std::vector<WinInfo>& win_cov, size_t start, size_t end) {
        double max_val = win_cov[start].max_c[1];
        size_t max_pos = start;
        // LOG(INFO)("fff %zd, %zd", start, end);

        for (size_t i = start + 1; i < end; ++i) {
            if (win_cov[i].max_c[1] >= max_val) {
                max_val = win_cov[i].max_c[1];
                max_pos = i;
            }
            // LOG(INFO)("fff %zd, %.02f -- %.02f %zd", i, win_cov[i].max_c[1], max_val, max_pos);
        }
        return std::make_pair(max_val, max_pos);
    };

    if (win_cov_.size() >= winnum) {
        std::pair<double,size_t> mx = {win_cov_[0].max_c[1], 0};
        
        for (size_t i = 1; i < winnum; ++i) {
            if (win_cov_[i].max_c[1] > mx.first) {
                mx.first = win_cov_[i].max_c[1];
                mx.second = i;
            }
        }
        max_covs.push_back(mx);

        for (size_t i = winnum; i < win_cov_.size(); ++i) {
            assert(max_covs.size() == i - winnum + 1);
            if (max_covs.size() == 225) {
                LOG(INFO)("bb %zd, %u, %zd %.02f", i, win_cov_[i].max_c[1], mx.second, mx.first) ;
            }
            if (max_covs.size() == 225) {
                LOG(INFO)("%zd < %zd", i - winnum, mx.second);
            }
            if (i - winnum < mx.second) {
                if (max_covs.size() == 225) {
                    LOG(INFO)("<< %d ", mx.first <= win_cov_[i].max_c[1]);
                    //assert(0);
                }
                if (mx.first <= win_cov_[i].max_c[1]) {
                    mx = {win_cov_[i].max_c[1], i};
                }
            } else {
                mx = find_max_posion(win_cov_, i - winnum + 1, i + 1);
            }
            if (max_covs.size() == 225) {
                LOG(INFO)("cc %.02f , %zd", mx.first, mx.second);
                //assert(0);
            }
            max_covs.push_back(mx);
        }
    }
    for (size_t i = 0; i < max_covs.size(); ++i) {
        double mx = win_cov_[i].max_c[1];
        size_t pos = i;
        for (size_t iw = i+1; iw < i + winnum; ++iw) {
            if (mx <= win_cov_[iw].max_c[1]) {
                mx = win_cov_[iw].max_c[1];
                pos = iw;
            }
        }
        if (!(max_covs[i].first == mx )) {
            LOG(INFO)("ccmax %zd %.02f %zd %.02f %zd", i, max_covs[i].first, max_covs[i].second, mx, pos);
            assert(max_covs[i].first == mx);
        }
    }
    return max_covs;

}


std::vector<std::pair<double,size_t>> WindowSlider::MinCoverages(size_t winnum) {
    std::vector<std::pair<double,size_t>> min_covs;

    auto find_min_posion = [](const std::vector<WinInfo>& win_cov, size_t start, size_t end) {
        double min_val = win_cov[start].min_c[1];
        size_t min_pos = start;
        for (size_t i = start + 1; i < end; ++i) {
            if (min_val >= win_cov[i].min_c[1]) {
                min_val = win_cov[i].min_c[1];
                min_pos = i;
            }
        }
        return std::make_pair(min_val, min_pos);
    };

    if (win_cov_.size() >= winnum) {
        std::pair<double,size_t> mn = {win_cov_[0].min_c[1], 0};
        
        for (size_t i = 1; i < winnum; ++i) {
            if (mn.first >= win_cov_[i].min_c[1]) {
                mn.first = win_cov_[i].min_c[1];
                mn.second = i;
            }
        }
        min_covs.push_back(mn);

        for (size_t i = winnum; i < win_cov_.size(); ++i) {

            if (i - winnum < mn.second) {
                if (mn.first >= win_cov_[i].min_c[1]) {
                    mn = {win_cov_[i].min_c[1], i};
                }
            } else {
                mn = find_min_posion(win_cov_, i - winnum + 1, i + 1);
            }
            min_covs.push_back(mn);
        }
    }
    
    for (size_t i = 0; i < min_covs.size(); ++i) {
        double mx = win_cov_[i].min_c[1];
        size_t pos = i;
        for (size_t iw = i+1; iw < i + winnum; ++iw) {
            if (mx > win_cov_[iw].min_c[1]) {
                mx = win_cov_[iw].min_c[1];
                pos = iw;
            }
        }
        if (!(min_covs[i].first == mx )) {
            LOG(INFO)("ccmin %zd %.02f %zd %.02f %zd", i, min_covs[i].first, min_covs[i].second, mx, pos);
            assert(min_covs[i].first == mx);
        }
    }

    return min_covs;

}

bool WindowSlider::HasBreakpoint(size_t s, size_t e) const {
    for (size_t i = s; i < e; ++i) {
        if (win_cov_[i].min_c[1] < 1) {
            return true;
        }
    }
    return false;
}
bool WindowSlider::HasAlternate(size_t s, size_t e) const {
    const int COV = 10;
    double clips = 0.0;
    double diff = 0.0;
    double cov = 1000;
    for (size_t i = s; i < e; ++i) {
        clips += win_cov_[i].clips;

        double r = win_cov_[i].count[0] == 0 ? 0.0 : win_cov_[i].c0 / win_cov_[i].count[0];
        diff = std::max<double>(diff, win_cov_[i].count[0] - win_cov_[i].c0);
        cov = std::min<double>(cov, win_cov_[i].c0 );
    }

    LOG(INFO)("reg_alter(%zd-%zd): %.02f, %.02f, %.02f", s, e, clips, cov, diff);

    return clips > std::min<double>(COV, cov/2) || diff > std::min<double>(COV, cov/2);


}
}