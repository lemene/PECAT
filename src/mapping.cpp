#include "mapping.hpp"

namespace fsa {


void Mapping::Load(const std::string &fname) {
    ol_store_.Load(fname);

}

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
    auto qid = ol_store_.GetStringPool().QueryIdByString(name);
    LOG(INFO)("QID = %d %zd", qid, queries_.size());
    if (qid != StringPool::NID) {
        auto range = queries_.find(qid);
        if (range != queries_.end()) {
            auto qt = sorted_by_start_[range->second[0]];

            LOG(INFO)("QT = %d %d %d", qt->b_.start, qt->b_.end, qt->a_.len);
            for (size_t i = range->second[1]; i < sorted_by_start_.size(); ++i) {
                auto qq = sorted_by_start_[i];

                LOG(INFO)("QQ = %d %d %d", qq->b_.start, qq->b_.end, qq->a_.len);
                if (qq->b_.end > qt->b_.start) {
                    printf("%s: %s\n", ol_store_.GetStringPool().QueryStringById(qq->a_.id).c_str(), qq->ToM4Line().c_str());
                }

                if (qq->b_.start > qt->b_.end) {
                    break;
                }
            }
        }
    }
}

} // namespace fsa