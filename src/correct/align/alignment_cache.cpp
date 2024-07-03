#include "alignment_cache.hpp"
#include "utils/logger.hpp"

namespace fsa {

bool AlignmentCache::GetAlignment(Seq::Id qid, Seq::Id tid, Alignment &al)  {
    auto iter = cache_.find(ToKey(tid, qid));
    if (iter != cache_.end()) {
        al = iter->second;
        cache_.erase(iter);
        if (al.Valid()) al.Swap();
        return true;
    } else {
        return false;
    }
}


void AlignmentCache::SetAlignment(Seq::Id qid, Seq::Id tid, const Alignment &r) {

    if (ids_.find(qid) != ids_.end() && ids_.find(tid) != ids_.end()) {
        cache_[ToKey(qid, tid)] = r;
    }
}

bool AlignmentCache::HasAlignment(Seq::Id qid, Seq::Id tid) const {
    return cache_.find(ToKey(qid, tid)) != cache_.end() || cache_.find(ToKey(tid, qid)) != cache_.end();
}

} // namespace fsa {