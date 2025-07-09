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
    }
}

auto WindowSlider::DetectErrorRegions(size_t max_gap, const std::array<double,3>& ave_covs) -> std::vector<ErrorRegion>  {

    // Detect error regions based on coverage information
    // - Detect regions with low coverage, high coverage, or significant differences between coverage types
    // 
    const size_t MIN_COV_HQ = 3;
    const double MAX_COV_DIFF = 0.2; 
    const double MIN_COV_DIFF = 0.1; 


    std::vector<ErrorRegion> cands;
    
    std::vector<double> l_slope;
    std::vector<double> l_slope_rate;
    std::vector<int> kernel = {5, -1, -1, -1, -1, -1};
    const int SLEN = kernel.size();
    assert(SLEN == kernel.size() && SLEN - 1 == kernel[0]);

    for (size_t i = 0; i + SLEN - 1 < win_cov_.size(); ++i) {
        int s = 0;
        int ss = 0;
        for (size_t ii = 0; ii < SLEN; ++ii) {
            s += int(win_cov_[i+ii].c1)*kernel[ii];
            ss += int(win_cov_[i+ii].c1);
        }
        l_slope.push_back(s);
        l_slope_rate.push_back(ss == 0 ? 0.0 : s * 1.0 / ss);
    }

    double mean = 0;
    double mad = 0;
    ComputeMedianAbsoluteDeviation(l_slope, mean, mad);
    LOG(INFO)("PEAK mean %.02f %.02f", mean,mad);
    auto low = mean - 6*1.253*mad;
    auto high = mean + 6*1.253*mad;

    
    double mean_rate = 0;
    double mad_rate = 0;
    ComputeMedianAbsoluteDeviation(l_slope_rate, mean_rate, mad_rate);
    LOG(INFO)("PEAK rate mean %.02f %.02f", mean_rate, mean_rate);
    auto low_rate = mean_rate - 6*1.253*mad_rate;
    auto high_rate = mean_rate + 6*1.253*mad_rate;

    for (size_t i = 0; i + SLEN - 1 < win_cov_.size(); ++i) {
        auto &winfo = win_cov_[i];
        auto p = l_slope[i];
        auto r = l_slope_rate[i];
        std::array<size_t,2> win = Window2Region(i);
        LOG(INFO)("PEAK %zd-%zd %.02f %.02f %.02f |  %.02f %.02f %.02f ", win[0], win[1], p, low, high, r, low_rate, high_rate);
        if (p < low && win_cov_[i].c1 <= ave_covs[1] || p > high && win_cov_[i].c1 >= ave_covs[1] ||
            (r < low_rate && win_cov_[i].c1 <= ave_covs[1] || r > high_rate && win_cov_[i].c1 >= ave_covs[1])) {
                
            LOG(INFO)("set_type_t: %zd-%zd %d", win[0], win[1], win_cov_[i].c1);
            winfo.type = 4;
        }
    }

        
    for (size_t i = 0; i < win_cov_.size(); ++i) {
        auto &winfo = win_cov_[i];
        std::array<size_t,2> win = Window2Region(i);
        if (winfo.min_c[1] * 1.0 < ave_covs[1] * 0.5 || winfo.max_c[1] * 1.0  > ave_covs[1] * 2) {
            if (winfo.type == 0) {
                if (winfo.c1 < winfo.count[1] *0.5) {
                    winfo.type = 2;
                } else {
                    winfo.type = 1;
                }
            } 
        } 

        if (win_cov_[i].c1 <= MIN_COV_HQ) {
            LOG(INFO)("add_cand_win0 %zd-%zd %.02f %.02f", win[0], win[1], win_cov_[i].c1 ,win_cov_[i].c0);
            //cands.push_back({win[0], win[1], 0});
            winfo.type = 2;
        }
        if ((win_cov_[i].c0- win_cov_[i].c1)*1.0 > std::max<double>(MIN_COV_HQ, win_cov_[i].c1*MAX_COV_DIFF)) {
            LOG(INFO)("add_cand_win1 %zd-%zd %zd %zd", win[0], win[1], win_cov_[i].c1 ,win_cov_[i].c0);    
            winfo.type = 2; 
        } else if ((win_cov_[i].c0- win_cov_[i].c1)*1.0 > std::max<double>(MIN_COV_HQ, win_cov_[i].c1*MIN_COV_DIFF)) {
            winfo.type = 3;
        }
        
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
    //return MergeRegions(cands, max_gap);
}

std::vector<ErrorRegion> WindowSlider::DetectSimpleRegions() const {
    
    std::vector<ErrorRegion> cands;

    for (size_t i = 0; i < win_cov_.size(); ++i) {
        std::array<size_t,2> win = {i*stride_, std::min(i*stride_+win_size_, cov_info_.Size()) };
        const auto& wc = win_cov_[i];

        if (wc.clips == 0 && wc.ins == 0) {

            cands.push_back({win[0], win[1], 0});
        }
    
    }
    return MergeRegions(cands, 0);
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


std::vector<ErrorRegion> WindowSlider::ExtendRegions(const std::vector<ErrorRegion> &regs, size_t max_gap) const  {
    
    const size_t EXT_SIZE = 5;
    std::vector<ErrorRegion> extended;
    for (const auto &r : regs) {
        auto win = Region2Window(r);
    }

    return extended;
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

//     auto check_adjacent_reg = [this](const ErrorRegion& reg0, const ErrorRegion& reg1) {
//         if (reg0.end + 1000000 < reg1.start) return false;

//         size_t s0 = reg0.start / stride_ - 1;
//         size_t e0 = 1+ (reg0.end > win_size_ ? reg0.end - win_size_ : 0) / stride_;
//         size_t s1 = reg1.start / stride_ - 1;
//         size_t e1 = 1 + (reg1.end > win_size_ ? reg1.end - win_size_ : 0) / stride_;

//         int t0 = 0;
//         if (win_cov_[s0].c1 + win_cov_[s0].c1*0.2 < win_cov_[e0].c1) {
//             t0 = -1;
//         } else if (win_cov_[s0].c1 > win_cov_[e0].c1 + win_cov_[e0].c1*0.2) {
//             t0 = 1;
//         }

//         int t1 = 0;
//         if (win_cov_[s1].c1 + win_cov_[s1].c1*0.2 < win_cov_[e1].c1) {
//             t1 = -1;
//         } else if (win_cov_[s1].c1 > win_cov_[e1].c1 + win_cov_[e1].c1*0.2) {
//             t1 = 1;
//         }
//         LOG(INFO)("mmm %zd-%zd %d, %zd-%zd %d", reg0.start, reg0.end, t0, reg1.start, reg1.end, t1);

//         if (t0 == -1 && t1 == 1) {
//             for (size_t i = e0; i <= s1; ++i) {
//                 LOG(INFO)("mmm-  %d %d %d", win_cov_[i].c1, win_cov_[s0].c1, win_cov_[e1].c1);
//                 if (win_cov_[i].c1 < (win_cov_[s0].c1 + win_cov_[e1].c1 )/ 2) {
//                     return false;
//                 }
//             }
//             return true;

//         }

//         if (t0 == 1 && t1 == -1) {
//             for (size_t i = e0; i <= s1; ++i) {
//                 if (win_cov_[i].c1 > (win_cov_[s0].c1 + win_cov_[e1].c1 )/ 2) {
//                     return false;
//                 }
//             }
//             return true;

//         }
//         return false;
//     };

//     std::vector<ErrorRegion> merged;
//     if (regs.size() > 0) {
//         merged.push_back(regs[0]);

//         for (size_t i = 1; i < regs.size(); ++i) {
//             if (check_adjacent_reg(merged.back(), regs[i])) {
//                 assert(merged.back().end <= regs[i].end);
//                 merged.back().end = regs[i].end;
//             } else {
//                 merged.push_back(regs[i]);
//             }
//         }
//     }
//     return merged;
// }



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
}