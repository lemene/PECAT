#pragma once

#include <string>
#include <unordered_set>
#include <unordered_map>

#include "overlap_store.hpp"

namespace fsa {

class Mapping {
public:
    Mapping(const OverlapStore &ols) : ol_store_(ols) {}

    std::unordered_set<Seq::Id> GetMappedReads() const;
    const std::string& QueryNameById(Seq::Id id) { return ol_store_.GetStringPool().QueryStringById(id); }
    void BuildIndex();
    void QueryOverlaps(const std::string &name);

    struct Pair {
        Overlap ToOverlap();

        const Overlap* query;
        const Overlap* target;
    };
    std::vector<Pair> QueryOverlaps(Seq::Id id);
    
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

} // namespace fsa