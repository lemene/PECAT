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
        if (sorted_by_start_[curr]->b_.id != sorted_by_start_[i]->b_.id) {
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

std::vector<Mapping::Pair> Mapping::QueryOverlaps(const std::string &name) const {
    auto id = ol_store_.GetStringPool().QueryIdByString(name);
    return id != StringPool::NID ? QueryOverlaps(id) : std::vector<Pair>();
}

std::vector<Mapping::Pair> Mapping::QueryOverlaps(Seq::Id id) const {
    assert (id != StringPool::NID) ;

    std::vector<Mapping::Pair> ols;
    auto tgt = queries_.find(id);
    if (tgt != queries_.end()) {
        
        for (size_t i = tgt->second[0]; i < tgt->second[1]; ++i) {
            auto range = query_ranges_[i];
            const Overlap* tol = sorted_by_start_[range[0]];

            for (size_t ii = range[1]; ii < sorted_by_start_.size(); ++ii) {
                if (ii == range[0]) continue;

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

    return ols;
}

Overlap Mapping::Pair::ToOverlap() const {
    Overlap ol;
    ol.a_.id = query->a_.id;
    ol.a_.len = query->a_.len;
    ol.a_.strand = query->SameDirect() == target->SameDirect() ? 0 : 1;
    
    ol.b_.id = target->a_.id;
    ol.b_.len = target->a_.len;
    ol.b_.strand = 0;
    ol.identity_ = std::min(query->identity_, target->identity_);

    size_t als_len = std::max<size_t>(query->b_.end, target->b_.end) - 
                     std::min<size_t>(query->b_.start, target->b_.start);
    std::vector<std::array<int,2>> als(als_len, {-1, -1});

    size_t als_start = std::min<size_t>(query->b_.start, target->b_.start);
    auto align_cigar = [](const Overlap* query, decltype(als)& als, size_t als_start, size_t idx) {
        size_t qcurr = 0;
        size_t tcurr = 0;
        for (const auto& d : query->detail_) {
            switch (d.type) {
            case 'M':
            case '=':
            case 'X':
                
                for (size_t i = 0; i < d.len; ++i) {
                    auto p = query->a_.strand == 0 ? (query->a_.start + qcurr + i) : (query->a_.end - qcurr - i);
                    
                    als[query->b_.start + tcurr + i-als_start][idx] = p;
                }
                qcurr += d.len;
                tcurr += d.len;
                break;

            case 'D':
                for (size_t i = 0; i < d.len; ++i) {
                    auto p = query->a_.strand == 0 ? query->a_.start + qcurr : query->a_.end - qcurr ;
                    als[query->b_.start + tcurr + i-als_start][idx] = p;
                }
                tcurr += d.len;
                break;

            case 'I':
                qcurr += d.len;
                break;

            default:
                LOG(ERROR)("Not support cigar type '%c'.", d.type);
            }
        }
    };

    align_cigar(query, als, als_start, 0);
    align_cigar(target, als, als_start, 1);


    // for (size_t i = 0; i < als.size(); ++i) {
    //     LOG(INFO)("-  %d %d", als[i][0], als[i][1]);
    // }

    const size_t N = 1;
    size_t start = std::max<size_t>(query->b_.start, target->b_.start) - als_start;
    size_t end = std::min<size_t>(query->b_.end, target->b_.end) - als_start;
    //LOG(INFO)("S-E0: %d-%d", start,end);
    for (; start + N < end; ++start) {
        if (std::abs(als[start][0] - als[start+N][0]) == N &&
            std::abs(als[start][1] - als[start+N][1]) == N) {
            
            break;
        }
    }

    for (; start + N < end; end--) {
        if (std::abs(als[end-1][0] - als[end-1-N][0]) == N &&
            std::abs(als[end-1][1] - als[end-1-N][1]) == N) {
            break;
        }
    }
        

    //LOG(INFO)("S-E1: %d-%d", start,end);
    
    //assert(end > start + 100);
    if (query->a_.strand == 0) {
        ol.a_.start = als[start][0];
        ol.a_.end = als[end-1][0] + 1;
    } else {
        ol.a_.start = als[end-1][0]-1;
        ol.a_.end = als[start][0];
    }
    
    if (target->a_.strand == 0) {
        ol.b_.start = als[start][1];
        ol.b_.end = als[end-1][1] + 1;
    } else {
//        printf("tss: %d %d %d %d\n", start, end, als[start][0]-1, als[end-1][0]);
        ol.b_.start = als[end-1][1]-1;
        ol.b_.end = als[start][1];
    }
    return ol;
    
}
} // namespace fsa