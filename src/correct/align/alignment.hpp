#ifndef FSA_CORRECT_ALIGNMENT_HPP
#define FSA_CORRECT_ALIGNMENT_HPP

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
    std::array<char,2> GetAlign(size_t i) const { return {aligned_query[i], aligned_target[i]}; }

    bool Valid() const { return target_end - target_start > 0; }
    void Swap(bool sameDirect=true); 
    void Rearrange() { Rearrange(aligned_query, aligned_target); }
    static void Rearrange(std::string &alq, std::string &alt);
    static void Rearrange1(std::string &alq, std::string &alt);
    bool TrimEnds(size_t checklen=2000, int stub=8);

    Seq::Id tid { Seq::NID };
    Seq::Id qid { Seq::NID };
    size_t target_start {0};
    size_t target_end {0};
    size_t query_start {0};
    size_t query_end {0};
    size_t distance {0};
    std::string aligned_target;
    std::string aligned_query;

    const DnaSeq* target { nullptr };
    const DnaSeq* query { nullptr };
};

} // namespace fsa {

#endif // FSA_CORRECT_ALIGNMENT_HPP