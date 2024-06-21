#include "mapping.hpp"

namespace fsa {


void Mapping::BuildIndex() {
    assert(sorted_by_start_.size() == 0 && "Rebuild");
    LOG(INFO)("sort overlaps by start position");
    sorted_by_start_.reserve(ol_store_.Size());
    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        sorted_by_start_.push_back(&ol_store_.Get(i));
    }

    std::sort(sorted_by_start_.begin(), sorted_by_start_.end(), [](const Overlap* a, const Overlap *b) {
        return a->b_.id < b->b_.id || (a->b_.id == b->b_.id && a->b_.start < b->b_.start);
    });

    if (sorted_by_start_.size() > 0) {
        BuildTargetIndex();
        BuildQueryRange();
        BuildQueryIndex();
    }
}

void Mapping::BuildTargetIndex() {
    LOG(INFO)("BuildTargetIndex");
    assert(sorted_by_start_.size() > 0);

    std::pair<Seq::Id, std::array<size_t, 2>> curr = {sorted_by_start_[0]->b_.id, {0, 1}};
    for (size_t i = 0; i < sorted_by_start_.size(); ++i) {
        if (sorted_by_start_[i]->b_.id != curr.first) {
            targets_.insert(curr);
            curr.first = sorted_by_start_[i]->b_.id;
            curr.second[0] = i;
            curr.second[1] = i+1;
        } else {
            curr.second[1] = i+1;
        }
    }
    targets_.insert(curr);
}

void Mapping::BuildQueryRange() {
    LOG(INFO)("BuildQueryRange");
    assert(sorted_by_start_.size() > 0);

    query_ranges_.reserve(sorted_by_start_.size());

    size_t curr = 0;
    for (size_t i = 0; i < sorted_by_start_.size(); ++i) {
        if (sorted_by_start_[curr]->b_.id == sorted_by_start_[i]->b_.id) {
            curr = i;
        }

        if (sorted_by_start_[curr]->b_.end <= sorted_by_start_[i]->b_.start) {
            for (; curr <= i; ++curr) {
                if (sorted_by_start_[curr]->b_.end > sorted_by_start_[i]->b_.start) {
                    break;
                }
            }
        }
        query_ranges_.push_back({i, curr});
    }
    LOG(INFO)("BuildQueryRange1");

    std::sort(query_ranges_.begin(), query_ranges_.end(), [this](const std::array<size_t, 2> &a, const std::array<size_t, 2> &b) {
        const Overlap* aol = sorted_by_start_[a[0]];
        const Overlap* bol = sorted_by_start_[b[0]];
        return (aol->a_.id < bol->a_.id) || 
               (aol->a_.id == bol->a_.id && aol->b_.id < bol->b_.id) ||
               (aol->a_.id == bol->a_.id && aol->b_.id == bol->b_.id  && aol->b_.start < bol->b_.start);
    });
}

void Mapping::BuildQueryIndex() {
    LOG(INFO)("BuildQueryIndex");
    
    assert(query_ranges_.size() > 0);
    std::pair<Seq::Id, std::array<size_t, 2>> curr = {sorted_by_start_[query_ranges_[0][0]]->a_.id, {0, 1}};
    for (size_t i = 0; i < sorted_by_start_.size(); ++i) {
        if (sorted_by_start_[query_ranges_[i][0]]->a_.id != curr.first) {
            queries_.insert(curr);
            curr.first = sorted_by_start_[query_ranges_[i][0]]->a_.id;
            curr.second[0] = i;
            curr.second[1] = i+1;
        } else {
            curr.second[1] = i+1;
        }
    }
    assert(queries_.find(curr.first) == queries_.end());
    queries_.insert(curr);
}

std::unordered_set<Seq::Id> Mapping::GetMappedReads() const {
    std::unordered_set<Seq::Id> mapped;

    for (size_t i = 0; i < ol_store_.Size(); ++i) {
        const auto& ol = ol_store_.Get(i);

        if (ol.AlignedSize() >= 0.9 * ol.a_.len) {
            mapped.insert(ol.a_.id);
        }
    }
    return mapped;
}

void Mapping::QueryOverlaps(const std::string &name) {
    auto id = ol_store_.GetStringPool().QueryIdByString(name);
    if (id != StringPool::NID) {
        QueryOverlaps(id);
    }
}

std::vector<Mapping::Pair> Mapping::QueryOverlaps(Seq::Id id) {
    assert (id != StringPool::NID) ;

    std::vector<std::array<const Overlap*, 2>> ols;

    auto tgt = queries_.find(id);
    if (tgt != queries_.end()) {
        
        for (size_t i = tgt->second[0]; i < tgt->second[1]; ++i) {
            auto range = query_ranges_[i];
            const Overlap* tol = sorted_by_start_[i];

            for (size_t ii = range[1]; ii < sorted_by_start_.size(); ++ii) {
                auto qol = sorted_by_start_[ii];
                const int offset = 3000;
                if (qol->b_.end >= tol->b_.start + offset && qol->b_.start + offset <= tol->b_.end ) {
                    ols.push_back({qol, tol});
                }

                if (qol->b_.start > tol->b_.end) {
                    break;
                }
            }
        }

    }

    
}

Overlap Mapping::Pair::ToOverlap() const {
    Overlap ol;
    ol.a_.id = query->a_.id;
    ol.a_.len = query->a_.len;
    o.a_.strand = query->SameDirect() == target->SameDirect() ? 0 : 1;
    
    ol.b_.id = target->a_.id;
    ol.b_.len = target->a_.len;
    ol.b_.strand = 0;
    
}
} // namespace fsa