#include "stat_readinfo.hpp"

#include "../utility.hpp"
#include "../utils/logger.hpp"

namespace fsa {


double CalcLocalOverhangThreshold(std::vector<std::array<double,2>> &overhangs, AsmOptions& opts_) {
        
    std::sort(overhangs.begin(), overhangs.end(), [](const std::array<double,2>& a, const std::array<double,2> &b){
        return  a[0]*1.0/ a[1] < b[0] *1.0/ b[1];
    });
    double overhang_threshold = opts_.filter0.max_overhang;

    
    double median = 0;
    double mad = 0;
    size_t size = std::min((int)overhangs.size(), 40);  // TODO
    std::vector<std::array<double,2>> ohs(size);
    std::transform(overhangs.begin(), overhangs.begin()+size, ohs.begin(),  [](const std::array<double,2> &a) {
        return std::array<double,2>({a[0], a[1]});
    });

    if (ohs.size() > 0) {
        ComputeMedianAbsoluteDeviation(ohs, median, mad);
        if (mad >= 1.0e-10) {
            overhang_threshold = median + 6*1.4826*mad;
            overhang_threshold = std::min<double>(overhang_threshold, opts_.filter0.max_overhang);

        } else {
            ComputeMeanAbsoluteDeviation(ohs, median, mad);
            overhang_threshold = median + 6*1.253*mad;
            overhang_threshold = std::min<double>(overhang_threshold, opts_.filter0.max_overhang);
        }
    }
    return overhang_threshold;
}


void ReadStatInfo::Stat(int i, AsmOptions& opts_) {
    auto& ori = *this;
    auto& nri = *this;
    nri.id = i;
    nri.count = ori.identities.size();
    nri.len = ori.len;

    const int winsize = 1000;
    std::vector<std::vector<double>> win_idents((ori.len + winsize  / 2) / winsize);
    std::sort(ori.identities.begin(), ori.identities.end(), [](const IdentityInfo& a, const IdentityInfo& b) {
        return (a.end - a.start) * a.identity > (b.end - b.start) * b.identity;
    });
    for (auto &i : ori.identities) {
        size_t s = (i.start + winsize / 2) / winsize;
        size_t e = (i.end + winsize / 2) / winsize;
        std::for_each(win_idents.begin()+s, win_idents.begin()+e, [&i](std::vector<double>& v) {
            v.push_back(i.identity);
        });
    }


    nri.identity_threshold.resize(win_idents.size());
    std::transform(win_idents.begin(), win_idents.end(), nri.identity_threshold.begin(), [this, i, opts_](std::vector<double>& a){
        if (a.size() > 0) {
            double median, mad;
            //std::sort(a.begin(), a.end(), [](double i0, double i1) { return i0 > i1; });

            std::vector<double> tmp(a.begin(), a.begin() + std::min<size_t>(opts_.coverage , a.size()));
            ComputeMedianAbsoluteDeviation(tmp, median, mad);
            return std::max(opts_.filter0.min_identity, median-6*1.4826*mad);
        } else {
            return opts_.filter0.min_identity;
        }
    });

    std::vector<std::array<double,2>> lohs;   
    std::vector<std::array<double,2>> rohs;
    for (auto &i : ori.overhangs) {
        if (i.overhang[0] >= 0) {
            lohs.push_back({(double)i.overhang[0], (double)i.alen});
        }
        if (i.overhang[1] >= 0) {
            rohs.push_back({(double)i.overhang[1], (double)i.alen});
        }
    }

    nri.overhang_l_threshold = CalcLocalOverhangThreshold(lohs, opts_);
    nri.overhang_r_threshold = CalcLocalOverhangThreshold(rohs, opts_);
}


void StatReadInfo::Add(Overlap& o) {

    read_infos_[o.a_.id].len = o.a_.len;
    read_infos_[o.b_.id].len = o.b_.len;

    read_infos_[o.a_.id].identities.push_back({(float)o.identity_, o.b_.len, o.a_.start, o.a_.end});
    read_infos_[o.b_.id].identities.push_back({(float)o.identity_, o.a_.len, o.b_.start, o.b_.end});
    
    auto oh = o.Overhang2();
    if (oh[0] != 0) {
        int16_t loh = (oh[0] & 0x1) == 0 ? -1 : o.a_.start - 0;
        int16_t roh = (oh[0] & 0x2) == 0 ? -1 : o.a_.len - o.a_.end;
        read_infos_[o.a_.id].overhangs.push_back({{loh, roh}, o.b_.len, o.a_.end - o.a_.start});
    }
    if (oh[1] != 0) {
        int16_t loh = (oh[1] & 0x1) == 0 ? -1 : o.b_.start - 0;
        int16_t roh = (oh[1] & 0x2) == 0 ? -1 : o.b_.len - o.b_.end;
        read_infos_[o.b_.id].overhangs.push_back({{loh, roh}, o.a_.len, o.b_.end - o.b_.start});
    }

}

void StatReadInfo::Merge(StatReadInfo &sri) {
    for (auto &ri : sri.read_infos_) {
        auto iter = read_infos_.find(ri.first);
        if (iter != read_infos_.end()) {
            iter->second.identities.insert(iter->second.identities.end(), std::make_move_iterator(ri.second.identities.begin()), std::make_move_iterator(ri.second.identities.end()));
            iter->second.overhangs.insert(iter->second.overhangs.end(), std::make_move_iterator(ri.second.overhangs.begin()), std::make_move_iterator(ri.second.overhangs.end()));
        } else {
            read_infos_.insert(ri);
        }
    }
    sri.Clear();
}

void StatReadInfo::Clear() {
    read_infos_.clear();
}



double ReadStatInfo::IdentityThreshold(int start, int end) const {
    const int winsize = 1000;

    size_t s = (start + winsize / 2) / winsize;
    size_t e = (end + winsize / 2) / winsize;
    assert (s >= 0 && s <= e && e <= identity_threshold.size());

    size_t count = 0;
    double sum = 0;
    for (auto i = s; i< e; ++i) {
        if (identity_threshold[i] > 0) {
            count ++;
            sum += identity_threshold[i];
        }
    }

    return sum / count;
}

bool ReadStatInfo::CheckIdentity(const Overlap& o) const  {
    if (identity_threshold.size() > 0) {
        const auto& r = o.GetRead(id);
        return o.identity_ >= IdentityThreshold(r.start, r.end);
    }
    return true;
     
}


bool ReadStatInfo::CheckOverhang(const Overlap& o) const  {
    
    if (identity_threshold.size() > 0) {
        const auto& r = o.GetRead(id);
        return o.identity_ >= IdentityThreshold(r.start, r.end);
    }
    return true;
     
}

} // namespace