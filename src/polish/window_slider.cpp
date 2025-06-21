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

        std::vector<int> ss(win[1] - win[0], 0);

        for (size_t i = win[0]; i < win[1]; ++i) {
            const auto& c = cov_info_.Get(i);

            cov0[i-win[0]] = std::accumulate(c.cbases, c.cbases+5, 0);
            cov1[i-win[0]] = std::accumulate(c.abases, c.abases+5, 0);
            
            wc.clips += c.clips;
            wc.ins += c.abases[5];
        }
        win_cov_[i].c0 = std::accumulate(cov0.begin(), cov0.end(), 0);
        win_cov_[i].c1 = std::accumulate(cov1.begin(), cov1.end(), 0);
        win_cov_[i].min_c = *std::min_element(cov0.begin(), cov0.end());
    }
}



auto WindowSlider::DetectErrorRegions(size_t max_gap) -> std::vector<ErrorRegion>  {

    std::vector<ErrorRegion> cands;

    for (size_t i = 0; i < win_cov_.size(); ++i) {
        std::array<size_t,2> win = {i*stride_, std::min(i*stride_+win_size_, cov_info_.Size()) };

        // 检测是否有为0的位置，即断开位置
        if (win_cov_[i].min_c == 0) {
            cands.push_back({win[0], win[1], 0});
        }

        // 检查c0是否急剧变化，

        // c0何c1有较大差别
        
        if ((win_cov_[i].c1- win_cov_[i].c0)*1.0 > win_cov_[i].c0*0.2) {
            cands.push_back({win[0], win[1], 1}); 
        }
    
    }
    return MergeRegions(cands, max_gap);
}

auto WindowSlider::DetectErrorRegions1(size_t max_gap) -> std::vector<ErrorRegion>  {

    const size_t MIN_COV_HQ = 3;
    const double MAX_COV_DIFF = 0.2; 
    const double MIN_COV_DIFF = 0.1; 

    std::vector<ErrorRegion> cands;

    for (size_t i = 0; i < win_cov_.size(); ++i) {
        std::array<size_t,2> win = Window2Region(i);

        if (win_cov_[i].c0 / win_size_ <= MIN_COV_HQ) {
            LOG(INFO)("add_cand_win0 %zd-%zd %zd %zd", win[0], win[1], win_cov_[i].c0 ,win_cov_[i].c1);
            cands.push_back({win[0], win[1], 0});
        }
        if ((win_cov_[i].c1- win_cov_[i].c0)*1.0 > std::max<double>(MIN_COV_HQ*win_size_, win_cov_[i].c0*MAX_COV_DIFF)) {
            LOG(INFO)("add_cand_win1 %zd-%zd %zd %zd", win[0], win[1], win_cov_[i].c0 ,win_cov_[i].c1);
            for (size_t ii = i >= 3 ? i - 3 : 0; ii < (i + 4 < win_cov_.size() ? i+3 : win_cov_.size()); ++ii) {
                if ((win_cov_[ii].c1- win_cov_[ii].c0)*1.0 > std::max<double>(MIN_COV_HQ*win_size_, win_cov_[ii].c0*MIN_COV_DIFF)) {
                    std::array<size_t,2> win = Window2Region(ii);
                    cands.push_back({win[0], win[1], 1}); 
                }
            }
        }
    
    }
    std::sort(cands.begin(), cands.end(), [](ErrorRegion& a, ErrorRegion &b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    return MergeRegions2(MergeRegions(cands, max_gap));
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


std::vector<ErrorRegion> WindowSlider::MergeRegions2(const std::vector<ErrorRegion> &regs) const {

    // merge
    auto check_adjacent_reg = [this](const ErrorRegion& reg0, const ErrorRegion& reg1) {
        if (reg0.end + 1000000 < reg1.start) return false;

        size_t s0 = reg0.start / stride_ - 1;
        size_t e0 = 1+ (reg0.end > win_size_ ? reg0.end - win_size_ : 0) / stride_;
        size_t s1 = reg1.start / stride_ - 1;
        size_t e1 = 1 + (reg1.end > win_size_ ? reg1.end - win_size_ : 0) / stride_;

        int t0 = 0;
        if (win_cov_[s0].c0 + win_cov_[s0].c0*0.2 < win_cov_[e0].c0) {
            t0 = -1;
        } else if (win_cov_[s0].c0 > win_cov_[e0].c0 + win_cov_[e0].c0*0.2) {
            t0 = 1;
        }

        int t1 = 0;
        if (win_cov_[s1].c0 + win_cov_[s1].c0*0.2 < win_cov_[e1].c0) {
            t1 = -1;
        } else if (win_cov_[s1].c0 > win_cov_[e1].c0 + win_cov_[e1].c0*0.2) {
            t1 = 1;
        }
        LOG(INFO)("mmm %zd-%zd %d, %zd-%zd %d", reg0.start, reg0.end, t0, reg1.start, reg1.end, t1);

        if (t0 == -1 && t1 == 1) {
            for (size_t i = e0; i <= s1; ++i) {
                LOG(INFO)("mmm-  %d %d %d", win_cov_[i].c0, win_cov_[s0].c0, win_cov_[e1].c0);
                if (win_cov_[i].c0 < (win_cov_[s0].c0 + win_cov_[e1].c0 )/ 2) {
                    return false;
                }
            }
            return true;

        }

        if (t0 == 1 && t1 == -1) {
            for (size_t i = e0; i <= s1; ++i) {
                if (win_cov_[i].c0 > (win_cov_[s0].c0 + win_cov_[e1].c0 )/ 2) {
                    return false;
                }
            }
            return true;

        }
        return false;
    };

    std::vector<ErrorRegion> merged;
    if (regs.size() > 0) {
        merged.push_back(regs[0]);

        for (size_t i = 1; i < regs.size(); ++i) {
            if (check_adjacent_reg(merged.back(), regs[i])) {
                assert(merged.back().end <= regs[i].end);
                merged.back().end = regs[i].end;
            } else {
                merged.push_back(regs[i]);
            }
        }
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
}