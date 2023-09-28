#include "alignment_cache.hpp"

namespace fsa {

bool AlignmentCache::GetAlignment(Seq::Id qid, Seq::Id tid, bool direct, Alignment &al)  {
    auto iter = cache_.find(ToKey(tid, qid, direct));
    if (iter != cache_.end()) {
        al = iter->second;
        cache_.erase(iter);
        if (al.Valid()) al.Swap(direct);
        return true;
    } else {
        return false;
    }
}


void AlignmentCache::SetAlignment(Seq::Id qid, Seq::Id tid, bool direct, const Alignment &r) {

    if (ids_.find(qid) != ids_.end() && ids_.find(tid) != ids_.end()) {
        cache_[ToKey(qid, tid, direct)] = r;
    }
}

bool AlignmentCache::HasAlignment(Seq::Id qid, Seq::Id tid, bool direct) const {
    return cache_.find(ToKey(qid, tid, direct)) != cache_.end() || cache_.find(ToKey(tid, qid, direct)) != cache_.end();
}

} // namespace fsa {