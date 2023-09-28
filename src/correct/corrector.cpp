#include "corrector.hpp"

namespace fsa {

void OverlapGrouper::BuildIndex(size_t thread_size, const std::unordered_set<int>& read_ids) {

    sorted_.assign(ol_store_.Size()*2, nullptr);

    for (size_t i = 0; i < ol_store_.Size(); ++i)  {
        auto ol = &ol_store_.Get(i);
        sorted_[i] = ol;
        sorted_[i+ol_store_.Size()] = ol;
    }
    std::sort(sorted_.begin(), sorted_.begin() + sorted_.size()/2, [](const Overlap* a, const Overlap *b) { 
        return (a->a_.id < b->a_.id) || 
               (a->a_.id == b->a_.id && a->b_.id < b->b_.id) ||
               (a->a_.id == b->a_.id && a->b_.id == b->b_.id && a->AlignedSize() > b->AlignedSize()) ||
               (a->a_.id == b->a_.id && a->b_.id == b->b_.id && a->AlignedSize() == b->AlignedSize() && a->SameDirect() && !b->SameDirect());
    });

    std::sort(sorted_.begin() + sorted_.size()/2, sorted_.end(), [](const Overlap* a, const Overlap *b) { 
        return (a->b_.id < b->b_.id) || 
               (a->b_.id == b->b_.id && a->a_.id < b->a_.id) ||
               (a->b_.id == b->b_.id && a->a_.id == b->a_.id && a->AlignedSize() > b->AlignedSize()) ||
               (a->b_.id == b->b_.id && a->a_.id == b->a_.id && a->AlignedSize() == b->AlignedSize() && a->SameDirect() && !b->SameDirect()) ;
    });

    assert(ol_store_.Size() > 0 && sorted_.size() == 2*ol_store_.Size());

    auto gp_s_qry = 0;
    auto gp_id_qry = sorted_[gp_s_qry]->a_.id;
    for (size_t i = 0; i < sorted_.size()/2; ++i) {
        if (gp_id_qry != sorted_[i]->a_.id) {
            auto &idx = index_[gp_id_qry];
            idx.by_qurey[0] = gp_s_qry;
            idx.by_qurey[1] = i;    // group_end

            gp_s_qry = i;
            gp_id_qry = sorted_[gp_s_qry]->a_.id;
        }
    }
    index_[gp_id_qry].by_qurey[0] = gp_s_qry;
    index_[gp_id_qry].by_qurey[1] = sorted_.size() / 2;

    auto gp_s_tgt = sorted_.size() / 2;             // group start position in target.
    auto gp_id_tgt = sorted_[gp_s_tgt]->b_.id;      // group id in target

    for (size_t i = sorted_.size()/2; i < sorted_.size(); ++i) {
        if (gp_id_tgt != sorted_[i]->b_.id) {
            auto &idx = index_[gp_id_tgt];
            idx.by_target[0] = gp_s_tgt;
            idx.by_target[1] = i;

            gp_s_tgt = i;
            gp_id_tgt = sorted_[gp_s_tgt]->b_.id;
        }
    }
    index_[gp_id_tgt].by_target[0] = gp_s_tgt;
    index_[gp_id_tgt].by_target[1] = sorted_.size();
}

auto OverlapGrouper::Get(int id) -> Group {
    Group group(id);
    auto &ols = group.ols;
    auto &index = group.index;

    auto info = index_.find(id);
    if (info != index_.end()) {
        const Index& sindex = info->second;

        if (sindex.by_qurey[1] > sindex.by_qurey[0]) {
            size_t s = ols.size();
            ols.insert(ols.end(), sorted_.begin() + sindex.by_qurey[0], sorted_.begin() + sindex.by_qurey[1]);

            index.push_back({s, ols.size()});
            for (size_t i = s; i < ols.size(); ++i) {
                if (ols[i]->b_.id != ols[index.back()[0]]->b_.id) {
                    index.back()[1] = i;
                    index.push_back({i, ols.size()});
                }
            }

        }

        if (sindex.by_target[1] > sindex.by_target[0]) {
            size_t s = ols.size();
            ols.insert(ols.end(), sorted_.begin() + sindex.by_target[0], sorted_.begin() + sindex.by_target[1]);

            index.push_back({s, ols.size()});
            for (size_t i = s; i < ols.size(); ++i) {
                if (ols[i]->a_.id != ols[index.back()[0]]->a_.id) {
                    index.back()[1] = i;
                    index.push_back({i, ols.size()});
                }
            }
        }
    }
    return group;
}

void OverlapGrouper::Group::Sort(double opt_ohwt) {
    assert(!Empty());

    auto weights = GetWeight(opt_ohwt);

    std::sort(index.begin(), index.end(), [&weights](const std::array<size_t, 2> &a, const std::array<size_t,2> &b) {
        return weights[a[0]] > weights[b[0]];
    });
}

std::vector<double> OverlapGrouper::Group::GetWeight(double opt_ohwt) {
    assert(!Empty());
    size_t target_length = ols[0]->GetRead(id).len;
    std::vector<double> cand_cov_wts (target_length+1);

    double wtsum = 0.0;
    for (size_t i = 0; i < Size(); ++i) {
        auto o = Get(i, 0);
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        double ohwt = opt_ohwt * o->identity_ / 100;
        double olwt = o->identity_ / 100;

        auto mr = o->MappingTo<2>(t, {0, q.len});
        auto start = std::max(0, mr[0] < mr[1] ? mr[0] : mr[1]);
        auto end =   std::min(t.len, mr[0] >= mr[1] ? mr[0] : mr[1]);
        // start -- t.start -- t.end -- end
        assert(start <= t.start && t.end <= end);

        cand_cov_wts[start]   += ohwt;
        cand_cov_wts[t.start] += (olwt - ohwt);
        cand_cov_wts[t.end]   -= (olwt - ohwt);
        cand_cov_wts[end]     -= ohwt;

        wtsum += olwt;
    }

    for (size_t i=1; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] += cand_cov_wts[i-1];
    }
    assert(std::abs(cand_cov_wts.back()) < 0.0000001);  // cand_cov_wts.back() == 0

    for (size_t i=0; i<cand_cov_wts.size(); ++i) {
        cand_cov_wts[i] = wtsum - cand_cov_wts[i];
    }

    std::vector<double> weights(ols.size(), 0.0);
    for (size_t i = 0; i < index.size(); ++i) {
        auto o = ols[index[i][0]];
        auto &t = o->GetRead(id);
        auto &q = o->GetOtherRead(id);

        double ohwt = opt_ohwt * o->identity_ / 100;
        double olwt = o->identity_ / 100;

        auto mr = o->MappingTo<2>(t, {0, q.len});
        auto start = std::max(0, mr[0] < mr[1] ? mr[0] : mr[1]);
        auto end =   std::min(t.len, mr[0] >= mr[1] ? mr[0] : mr[1]);
        // start -- t.start -- t.end -- end
        assert(start <= t.start && t.end <= end);

        double wt = std::accumulate(cand_cov_wts.begin()+start, cand_cov_wts.begin()+t.start, 0.0) * ohwt +
                    std::accumulate(cand_cov_wts.begin()+t.start, cand_cov_wts.begin()+t.end, 0.0) * olwt + 
                    std::accumulate(cand_cov_wts.begin()+t.end, cand_cov_wts.begin()+end, 0.0) * ohwt;

        weights[index[i][0]] = wt;
    }
    return weights;
}

} // namespace fsa
