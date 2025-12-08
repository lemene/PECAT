#pragma once
#include "match_info.hpp"

#include "overlap.hpp"
#include "sequence.hpp"

namespace fsa {


class MatchInfo {
public:
    MatchInfo(const Overlap* ol, const DnaSeq& q, const DnaSeq& t);
    struct InfoItem {
        InfoItem() : ref(0), base(0), ins(0) {}
        uint32_t ref : 2;       // 0,1,2,3 ->A,C,G,T    base in reference
        uint32_t base : 3;   // 0,1,2,3,4 ->A,C,G,T,-   base in query
        uint32_t ins : 19 ;  // index of insert 
        uint32_t pos;       // position in query
    };


    size_t Start() const { return ol_->b_.start; }
    size_t End() const { return ol_->b_.end; }
    size_t Len() const { return ol_->b_.len; }
    size_t LClip() const { 
        auto qclip = ol_->SameDirect() ? ol_->a_.start : ol_->a_.len - ol_->a_.end; 
        return std::min(qclip, ol_->b_.start);
    }

    size_t RClip() const { 
        auto qclip = !ol_->SameDirect() ? ol_->a_.start : ol_->a_.len - ol_->a_.end; 
        return std::min(qclip, ol_->b_.len - ol_->b_.end);
    }


    size_t Size() const { return match_.size(); }
    const InfoItem& Get(size_t i) const  { return match_[i]; }
    const size_t GetInssize(size_t i) const { return i == 0 ? 0 : insert_[i-1][1]-insert_[i-1][0]; }
    double MaxLocalDistance() const { return max_local_distance_; }
    void CalculateMaxLocalDistance(size_t win_size) ;
    std::vector<std::array<size_t,2>> LocalDistance(size_t win_size) const;
    std::vector<double> LocalIdentity(size_t win_size) const;
    double Identity() const { return ol_->Identity(); }
    double MatchedIdentity() const { return matched_identity_; }
    std::vector<std::array<size_t,2>> GetHighQualityRegions(size_t win_size, double max_dist, size_t min_clip, size_t min_intv) const;
    std::vector<std::array<size_t,2>> GetHighQualityRegions(size_t half_win_size, double max_dist) const;
    const Overlap* GetOverlap() const { return ol_; }
    std::array<size_t,2> GetQueryRegion(size_t s, size_t e) const;

protected:
    void ComputeMatchedIdentity(size_t large_size = 500);
protected:
    std::vector<InfoItem> match_;
    std::vector<std::array<size_t,2>> insert_;
    const Overlap* ol_ { nullptr };
    double matched_identity_;   // without larget indels
    double max_local_distance_ { 0.0 };
};
}
