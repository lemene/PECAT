#pragma once

#include <unordered_map>
#include <unordered_set>

#include "alignment.hpp"
#include "../../sequence.hpp"

namespace fsa {

class AlignmentCache {
public:
    bool GetAlignment(Seq::Id qit, Seq::Id tid, bool isSameDirect, Alignment &al) ;
    void SetAlignment(Seq::Id qit, Seq::Id tid, bool isSameDirect, const Alignment &r);
    bool HasAlignment(Seq::Id qit, Seq::Id tid, bool isSameDirect) const;
    void Clear() { cache_.clear(); ids_.clear();}
    
    void Reset(const std::vector<Seq::Id> &ids, size_t size) { 
        Clear();
        ids_.insert(ids.begin(), ids.begin()+size); 
    }
protected:
    long long int ToKey(Seq::Id qit, Seq::Id tid, bool same_direct) const {
        long long int k = ((long long int)qit << 32) + tid;
        assert(k > 0);  // TODO Save one bit
        return same_direct ? k : -k;
    }

protected:
    //std::unordered_map<std::array<Seq::Id, 2>, Alignment, ArrayHash<Seq::Id, 2>, ArrayEqual<Seq::Id, 2>> cache_;
    std::unordered_map<long long int, Alignment> cache_;
    std::unordered_set<Seq::Id> ids_;
};

} // namespace fsa {
