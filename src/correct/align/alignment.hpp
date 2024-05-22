#pragma once

#include "../../sequence.hpp"

namespace fsa {


class Alignment {
public:
    Alignment(Seq::Id t, Seq::Id q) : tid(t), qid(q) { }
    Alignment(const DnaSeq* t = nullptr, const DnaSeq* q = nullptr)
        : target(t), query(q) {
    }

    void Reset(const DnaSeq* t = nullptr, const DnaSeq* q = nullptr) {
        target = t; query = q;
        aligned_target = "";
        aligned_query = "";
        target_end = target_start = 0;
    }

    size_t AlignSize()  const { return target_end - target_start; }
    size_t TargetSize() const { return target == nullptr ? 0 : target->Size(); }
    size_t QuerySize() const { return query == nullptr ? 0 : query->Size(); }
    double Identity() const { return 100*(1- distance * 1.0 / aligned_target.size()); }
    double IdentityIgnoreHomo(size_t len) const;
    std::array<char,2> GetAlign(size_t i) const { return {aligned_query[i], aligned_target[i]}; }

    bool Valid() const { return target_end - target_start > 0; }
    void Swap(bool sameDirect=true); 
    void Rearrange() { Rearrange(aligned_query, aligned_target); }
    static void Rearrange(std::string &alq, std::string &alt);
    static void Rearrange1(std::string &alq, std::string &alt);
    bool TrimEnds(size_t checklen=2000, int stub=8);
    
    void ComputeDistance(size_t win_size);
    uint16_t MaxLocalDistance() const { return local_distances[max_local_distance_position]; }
    std::pair<bool, uint16_t> MaxLocalDistance (size_t s, size_t e) const;
    size_t MaxLocalDistancePosition() const { return max_local_distance_position + target_start; }
    double MaxLocalIdentity_100(size_t win_size) const { return 100.0 - (MaxLocalDistance() * 100.0 / win_size); }

    Seq::Id tid { Seq::NID };
    Seq::Id qid { Seq::NID };
    size_t target_start {0};
    size_t target_end {0};
    size_t query_start {0};
    size_t query_end {0};
    size_t distance {0};
    std::string aligned_target;
    std::string aligned_query;

    std::vector<int16_t> local_distances;
    size_t max_local_distance_position { 0 };

    const DnaSeq* target { nullptr };
    const DnaSeq* query { nullptr };
};

} // namespace fsa {
