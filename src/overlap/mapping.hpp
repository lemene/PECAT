#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "overlap_store.hpp"

namespace fsa {

class Mapping {
public:
    Mapping(const OverlapStore &ols) : ol_store_(ols) {}

    std::unordered_set<Seq::Id> GetMappedReads() const;
    const std::string& QueryNameById(Seq::Id id) { return ol_store_.GetStringPool().QueryStringById(id); }
    void BuildIndex();

    struct Pair : public Overlap {
        Pair(const Overlap* q, const Overlap* t) : query(q), target(t) {
            ToOverlap();
        }

        // 0: match, 1: insert, 2: delect, 3: mismatch; consistent with edlib
        std::vector<uint8_t> AlignBases(Seq::Id tid, const DnaSeq &qseq, const DnaSeq &tseq);
        std::vector<uint8_t> AlignBases00(const DnaSeq &qseq, const std::vector<int> &al_q_2_ref, 
                                          const DnaSeq &tseq, const std::vector<int> &al_t_2_ref);
protected:
        void ToOverlap();

        std::vector<uint8_t> RealignBlock(
            const DnaSeq &qseq, const std::vector<int> &al_q_2_ref, 
            const DnaSeq &tseq, const std::vector<int> &al_t_2_ref, 
            size_t start, size_t end);
public:
        const Overlap* query {nullptr};
        const Overlap* target {nullptr};
        std::vector<int> aligned_target;
        std::vector<int> aligned_query; 
        size_t start {0};
        size_t end {0};
    };

    std::vector<Pair> QueryOverlaps(const std::string &name) const;
    std::vector<Pair> QueryOverlaps(Seq::Id id) const;
    
protected:
    void BuildTargetIndex();
    void BuildQueryRange();
    void BuildQueryIndex();
protected:

    const OverlapStore &ol_store_;
    std::vector<const Overlap*> sorted_by_start_;
    std::unordered_map<Seq::Id, std::array<size_t, 2>> targets_;
    std::unordered_map<Seq::Id, std::array<size_t, 2>> queries_;
    
    std::vector<std::array<size_t, 2>> query_ranges_;

};

class QueryGrouper {
public:
    QueryGrouper(const OverlapStore &ols) : ol_store_(ols) {}

    void BuildIndex();
    size_t QuerySize() const { return queries_.size(); }
    std::array<size_t, 2> GetQuery(size_t i) { return queries_[i]; }
    const Overlap* GetOverlap(size_t i) { return sorted_[i]; }

protected:

    const OverlapStore &ol_store_;
    std::vector<const Overlap*> sorted_;
    std::unordered_map<Seq::Id, std::array<size_t,2>> index_;
    std::vector<std::array<size_t,2>> queries_;

};

} // namespace fsa