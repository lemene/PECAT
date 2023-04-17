#include "local_phaser.hpp"
#include <algorithm>

#include "read_haplotype.hpp"
#include "../read_store.hpp"

namespace fsa {

static thread_local bool local_print_rubbish = false;
static thread_local const ReadStore* local_read_store = nullptr;

LocalPhaser::LocalPhaser(const PhsOptions::PhaserOptions& opts, ReadOffset tid, const std::unordered_set<ReadOffset> &qids, 
    const std::unordered_map<ReadOffset, ReadInfo>& rdinfos,
    const ReadStore &read_store, bool correct) 
 :  opts_(opts), rd_store_(read_store), correct_(correct) {

    local_print_rubbish = read_store.QueryNameById(tid.id) == opts_.debug;
    local_read_store = &read_store;

    DEBUG_local_printf("---- stat : %s\n", read_store.QueryNameById(tid.id).c_str());
    target_ = &rdinfos.find(tid)->second;
    for (size_t i=0; i<target_->vars.size(); ++i) {
        mapper[target_->vars[i][0]] = i;
    }

    for (const auto &m : mapper) {
        DEBUG_local_printf("mapper: %d %d\n", m.first, m.second);
    }

    coverages_.assign(target_->vars.size(), 0);
    for (auto i : qids) {
        auto q = &rdinfos.find(i)->second;
        queries_.push_back({q, MapVariants(q->vars)});
    }    
    queries_.push_back({target_, MapVariants(target_->vars)});
    
    for (auto &q : queries_) {
        if (q.query == target_) {
            target_in_queires_ = &q;
        }
    }
}

std::vector<int> LocalPhaser::MapVariants(const std::vector<std::array<int, 5>>& vars) {
    std::vector<int> mapped(target_->vars.size(), 0);

    for (auto & v : vars) {
        auto iter = mapper.find(v[0]);
        if (iter != mapper.end()) {
            if (v[3] != -1 && target_->vars[iter->second][3] != -1) {
                if (v[3] == target_->vars[iter->second][3]) {
                    mapped[iter->second] = 1;
                } else {
                    mapped[iter->second] = -1;
                }
            } else {
                mapped[iter->second] = 0;
            }

            if (v[3] != -1) {
                coverages_[iter->second] ++;
            }
        }
    }

    return mapped;
}


void PrintCentroid(const std::vector<int>& cs, const std::string &msg) {
    DEBUG_local_printf("%s:", msg.c_str());
    for (auto c : cs) {
        DEBUG_local_printf("%2d,", c);
    }
    DEBUG_local_printf("\n");
}

auto LocalPhaser::Group::SelectCentroids() const -> std::vector<Centroid>   {
    std::vector<const Query*> centroids { nullptr, nullptr };

    std::vector<std::array<int,3>> vscores ((*queries.begin())->vars.size(), {0, 0, 0});

    for (auto q : queries) {
        for (size_t i = 0; i < vscores.size(); ++i) {
            vscores[i][1 + q->vars[i]]++;
        }
    }

    // Eliminate randomness
    std::vector<const Query*> qs(queries.begin(), queries.end());
    std::sort(qs.begin(), qs.end(), [](const Query* a, const Query* b) {
        return a->query->id < b->query->id;
    });

    double best_score0 = 0;
    for (auto q : qs) {
        double score = 0.0;
        for (size_t i = 0; i < q->vars.size(); ++i) {
            if (q->vars[i] != 0) {
                if (vscores[i][0] != 0 && vscores[i][2] != 0) {
                    score += vscores[i][0] + vscores[i][2] - std::abs(vscores[i][0] - vscores[i][2]);
                }
            }
        }
        if (centroids[0] == nullptr || score > best_score0) {
            centroids[0] = q;
            best_score0 = score;
        }
    }

    size_t best_score1 = 0;
    for (auto q : qs) {
        auto s = centroids[0]->Distance(q);
        if (centroids[1] == nullptr || s > best_score1) {
            centroids[1] = q;
            best_score1 = s;
        }
    }

    PrintCentroid(centroids[0]->vars, "init query0");
    PrintCentroid(centroids[1]->vars, "init query1");

    return { Group(centroids[0]).GetCentroid(), Group(centroids[1]).GetCentroid() };
}



void LocalPhaser::Group::Centroid::SetRange2(const Group& g, const Centroid& alt, const Group& alt_g, const PhsOptions::PhaserOptions& opts) {
    
    auto cs = Values(opts.min_support_rate, opts.min_support_count, alt);
    auto alt_cs = alt.Values(opts.min_support_rate, opts.min_support_count, *this);

    std::vector<std::array<size_t, 3>> ranges = SplitRange(cs, g, alt_cs, alt_g, opts.min_link_support_rate, opts.min_link_support_count);
    for (size_t i = 0; i < ranges.size(); ++i) {
        DEBUG_local_printf("rs0(%zd) %zd %zd\n", i, ranges[i][0], ranges[i][1]);
    }

    std::vector<std::array<size_t, 3>> finals = CombineRange(ranges, cs, g, alt_cs, alt_g, opts);
 
    range[0] = g.vars.size();
    range[1] = 0;
    for (size_t i = 0; i < finals.size(); ++i) {
        DEBUG_local_printf("rs1(%zd) %zd %zd\n", i, finals[i][0], finals[i][1]);
        if (finals[i][1] - finals[i][0] >= (size_t)opts.min_inconsist_count) {
            range[0] = std::min(range[0], finals[i][0]);
            range[1] = std::max(range[1], finals[i][1]);
        }
    }
    if (range[0] > range[1]) {
        range[0] = 0;
        range[1] = 0;
    }

    // std::sort(finals.begin(), finals.end(), [](const std::array<size_t, 3>& a, const std::array<size_t, 3>& b) {
    //     return a[2] > b[2];
    // });

    // if (finals.size() > 0) {
    //     range[0] = finals[0][0];
    //     range[1] = finals[0][1];
    // } else {
    //     range[0] = 0;
    //     range[1] = 0;
    // }
}

std::vector<std::array<size_t, 3>> LocalPhaser::Group::Centroid::SplitRange(const std::vector<int> &cs, const Group &group, const std::vector<int> &alt_cs, const Group &alt_group, double rate, int count)
{
    std::vector<std::array<size_t, 3>> ranges;
    size_t start = -1;
    size_t curr = -1;
    size_t num = 0;
    for (size_t i = 0; i < cs.size(); ++i) {
        if (cs[i] != 0) {
            if (start == (size_t)-1) {
                start = i;
                curr = i;
                num = 1;
            } else {
                std::array<size_t, 2> stat = {0, 0};
                for (auto q : group.queries) {
                    if (q->vars[curr] != 0 && q->vars[i] != 0) {
                        if (q->vars[curr] == cs[curr] && q->vars[i] == cs[i]) {
                            stat[0] ++;
                        } else {
                            stat[1] ++;
                        }
                    }
                }
                DEBUG_local_printf("SetRange: (%zd - %zd) = (%zd %zd)\n", curr, i, stat[0], stat[1]);
                if (stat[0] - stat[1] < std::max<size_t>(count, std::ceil(rate*(stat[0]+stat[1])))) {
                    
                    DEBUG_local_printf("SetRange -- : %d,%d  %d,%d\n", alt_cs[curr] , cs[curr], alt_cs[i], cs[i]);
                    if (alt_cs[curr] * cs[curr] == -1 && alt_cs[i] * cs[i] == -1) {
                        
                        std::array<size_t, 2> alt_stat = {0, 0};
                        for (auto q : alt_group.queries) {
                            if (q->vars[curr] != 0 && q->vars[i] != 0) {
                                if (q->vars[curr] == alt_cs[curr] && q->vars[i] == alt_cs[i]) {
                                    alt_stat[0] ++;
                                } else {
                                    alt_stat[1] ++;
                                }
                            }
                        }
                DEBUG_local_printf("SetRange-alt: (%zd - %zd) = (%zd %zd)\n", curr, i, alt_stat[0], alt_stat[1]);
                        stat[0] += alt_stat[0];
                        stat[1] += alt_stat[1];

                    }
                }
                //DEBUG_local_printf("SetRange2: (%zd - %zd) = (%zd %zd)\n", curr, i, stat[0], stat[1]);

                if (stat[0] - stat[1] >= std::max<size_t>(count, std::ceil(rate*(stat[0]+stat[1])))) {
                    curr = i;
                    num++;
                } else {
                    ranges.push_back({start, curr+1, num});
                    start = i;
                    curr = i;
                    num = 1;
                }
            }
        }
    }
    if (start != (size_t)-1) {
        ranges.push_back({start, curr+1, num});
    }

    return ranges;
}

std::vector<std::array<size_t,3>> LocalPhaser::Group::Centroid::CombineRange(const std::vector<std::array<size_t,3>>& ranges, const std::vector<int> &cs, const Group& group, const std::vector<int>& alt_cs, const Group& alt_group, const PhsOptions::PhaserOptions& opts) {
 
    auto linkmaps = GetRangeLink(ranges, cs, group, opts.min_link_valid_rate, opts.min_link_valid_count);
    auto alt_linkmaps = GetRangeLink(ranges, alt_cs, alt_group, opts.min_link_valid_rate, opts.min_link_valid_count);

    std::vector<std::array<size_t, 3>> finals;

    // for (size_t ir = 0; ir < ranges.size();) {
    //     size_t end = ir;
    //     for (size_t jr = ir + 1; jr < ranges.size(); ++jr) {
    //         int count = opts.min_link_support_count;
    //         int rate = opts.min_link_support_rate;
    //         std::array<size_t,2> stat = linkmaps[jr * ranges.size() + ir];
    //         if (stat[0] - stat[1] >= std::max<size_t>(count, std::ceil(rate*(stat[0]+stat[1])))) {
    //             end = jr;
    //         } else {
    //             const auto &alt_stat = alt_linkmaps[jr * ranges.size() + ir];
    //             stat[0] += alt_stat[0];
    //             stat[1] += alt_stat[1];

    //             if (stat[0] - stat[1] >= std::max<size_t>(count, std::ceil(rate*(stat[0]+stat[1])))) {
    //                 end = jr;
    //             }
    //         }
    //         DEBUG_local_printf("map (%zd %zd) %zd %zd\n", ir, jr, stat[0], stat[1]);
            
    //         DEBUG_local_printf("map %d %d %zd\n", count, (int)std::ceil((rate-(1-rate))*(stat[0]+stat[1])), end);
    //     }

    //     std::array<size_t, 3> nr = ranges[ir]; ir++;
    //     for (; ir <= end; ++ir) {
    //         nr[1] = ranges[ir][1];
    //         nr[2] += ranges[ir][2];
    //     }
    //     finals.push_back(nr);
    // }

    std::vector<std::array<size_t,2>> linked;
    for (size_t ir = 0; ir < ranges.size(); ++ir) {
        size_t end = ir;
        for (size_t jr = ir + 1; jr < ranges.size(); ++jr) {
            int count = opts.min_link_support_count;
            double rate = opts.min_link_support_rate;
            std::array<size_t,2> stat = linkmaps[jr * ranges.size() + ir];
            
            if ((int)stat[0] - (int)stat[1] >= std::max<int>(count, std::ceil(rate*(stat[0]+stat[1])))) {
                end = jr;
            } else {
                const auto &alt_stat = alt_linkmaps[jr * ranges.size() + ir];
                stat[0] += alt_stat[0];
                stat[1] += alt_stat[1];

                if ((int)stat[0] - (int)stat[1] >= std::max<int>(count, std::ceil(rate*(stat[0]+stat[1])))) {
                    end = jr;
                }
            }
            DEBUG_local_printf("map (%zd %zd) %zd %zd\n", ir, jr, stat[0], stat[1]);
            
            DEBUG_local_printf("map %d %d %zd\n", count, (int)std::ceil((rate-(1-rate))*(stat[0]+stat[1])), end);
        }
        linked.push_back({ir, end});
        DEBUG_local_printf("linked: %zd - %zd\n", ir, end);
    }

    size_t irange = 0;
    for (size_t i = 1; i < linked.size(); ++i) {
        assert(i > irange && linked[i][0] > linked[irange][0]);
        if (linked[i][0] <= linked[irange][1]) {
            linked[irange][1] = std::max<size_t>(linked[i][1], linked[irange][1]);
        } else {
            irange ++;
            linked[irange] = linked[i];
        }
    }

    assert(irange <= linked.size());
    if (ranges.size() > 0) {
    for (size_t i = 0; i <= irange; ++i) {

        std::array<size_t, 3> nr = ranges[linked[i][0]]; 
        for (size_t j = linked[i][0]+1; j <= linked[i][1]; ++j) {
            nr[1] = ranges[j][1];
            nr[2] += ranges[j][2];
        }
        DEBUG_local_printf("Final Range: (%zd - %zd) (%zd, %zd, %zd)\n", linked[i][0], linked[i][1], nr[0], nr[1], nr[2]);
        finals.push_back(nr);
    }
    }

    return finals;
}

std::vector<std::array<size_t,2>> LocalPhaser::Group::Centroid::GetRangeLink(const std::vector<std::array<size_t,3>>& ranges, const std::vector<int> &cs, const Group& group, double rate, int count) {
 
    std::vector<std::array<size_t,2>> linkmaps(ranges.size()*ranges.size(), {0, 0});
    for (auto q : group.queries) {
        std::vector<std::array<size_t,3>> bs(ranges.size());
        for (size_t ir = 0; ir < ranges.size(); ++ir) {
            std::array<size_t, 3> stat = {0, 0, 0};
            for (size_t p = ranges[ir][0]; p < ranges[ir][1]; ++p) {
                stat[cs[p] * q->vars[p] + 1] += 1;
            }
            bs[ir] = stat;

            DEBUG_local_printf("ToC: %s | %zd = %d %d\n", local_read_store->QueryNameById(q->query->id).c_str(), ir, stat[0], stat[2]);
        }
        for (size_t ib = 0; ib < bs.size(); ++ib) {
            if (bs[ib][0] + bs[ib][2] < (size_t)count) continue;
            for (size_t jb = 0; jb < ib; ++jb) {
                if (bs[jb][0] + bs[jb][2] < (size_t)count) continue;

                std::array<size_t,2> &stat = linkmaps[ib*ranges.size() + jb];
                size_t th_ib = std::max<size_t>(count, std::ceil((bs[ib][0]+bs[ib][2])*rate));
                size_t th_jb = std::max<size_t>(count, std::ceil((bs[jb][0]+bs[jb][2])*rate));
                if (bs[ib][2] - bs[ib][0] >= th_ib && bs[jb][2] - bs[jb][0] >= th_jb) {
                    stat[0] += 1;
                } else {
                    stat[1] += 1;
                }
                //} else if (bs[ib][0] - bs[ib][2] >= th_ib && bs[jb][0] - bs[jb][2] >= th_jb) {
                //    stat[1] += 1;
                //}
            }
        }
    }
    return linkmaps;

}

auto LocalPhaser::Group::Centroid::Compare(const Query &q) const -> CompareResult {
    assert(locs.size() == q.vars.size());

    CompareResult rs;
    //for (size_t i = 0; i < locs.size(); ++i) {
    for (size_t i = range[0]; i < range[1]; ++i) {
        if (locs[i].cov > 0) {
            rs.asize ++;
        }
        if (q.vars[i] != 0) {
            rs.bsize ++;
        }

        if (locs[i].cov > 0 && q.vars[i] != 0) {
            rs.common_size ++;

            auto v = locs[i].accu <= 0 ? -1 : 1;
            if (v * q.vars[i] > 0) {
                rs.same ++;
            } else {
                rs.diff ++;
            }
        }
    }
    return rs;
}

std::vector<int> LocalPhaser::Group::Centroid::Values(double rate, int count) const {
    assert(locs.size() > 0);

    std::vector<int> centroid(locs.size());
    for (size_t i = 0; i < centroid.size(); ++i) {
        centroid[i] = locs[i].GetValidType(rate, count);
    }

    return centroid;
}

std::vector<int> LocalPhaser::Group::Centroid::Values(double rate, int count, const Centroid& alt) const {
    assert(locs.size() > 0 && locs.size() == alt.locs.size());

    std::vector<int> centroid(locs.size());
    for (size_t i = 0; i < centroid.size(); ++i) {
        centroid[i] = locs[i].GetValidType(rate, count);
        if (centroid[i] == 0) { // 是否互斥
            auto a = locs[i].GetNotableType(rate);
            auto b = alt.locs[i].GetNotableType(rate);
            if (a * b == -1) {
                Loc c(locs[i].accu - alt.locs[i].accu, locs[i].cov + alt.locs[i].cov);
                centroid[i] = c.GetValidType(rate, count);
            }
        }
    }

    return centroid;
}

size_t LocalPhaser::DistributeQuery(const std::vector<Group::Centroid>& cs, const Query* q) const {

    std::vector<std::pair<size_t, double>> result;
    for (size_t i = 0; i < cs.size(); ++i) {
        auto cmp = cs[i].Compare(*q);

        if ((int)cmp.common_size >= opts_.min_similar_count && cmp.Coverage() >= opts_.min_similar_coverage && cmp.Similary() >= opts_.min_similar_rate) { // TODO
            result.push_back({i, cmp.Similary()});
        }
    }
    
    assert(result.size() <= 2);
    if (result.size() == 1) {
        return result[0].first;
    } else if (result.size() == 2) {
        if (result[0].second > result[1].second + opts_.min_similar_diff) {
            return result[0].first;
        } else if (result[0].second + opts_.min_similar_diff < result[1].second) { 
            return result[1].first;
        } else {
            return cs.size();
        }

    } else {
        //
        return cs.size();   
    }
}


auto LocalPhaser::SplitGroups(const Group& g0, const std::vector<Group::Centroid> &centroids)
    -> std::vector<Group> {

    std::vector<Group> groups(centroids.size() + 1);

    for (auto q : g0.queries) {
        auto g = DistributeQuery(centroids, q);
        groups[g].Merge(Group(q));
    }

    return groups;
}

size_t LocalPhaser::AdjustGroups(std::vector<Group> &groups) {
    assert(groups.size() == 3 && groups[0].Size() > 0 && groups[1].Size() > 0); // group0, group1, untag

    std::vector<Group::Centroid> centroids = { groups[0].GetCentroid(), groups[1].GetCentroid() };
    centroids[0].SetRange2(groups[0], centroids[1], groups[1], opts_);
    centroids[1].SetRange2(groups[1], centroids[0], groups[0], opts_);

    std::unordered_map<const Query*, std::array<size_t,2>> changes;
    for (size_t ig = 0; ig < groups.size(); ++ig) {
        for (const auto q : groups[ig].queries) {
            auto ng = DistributeQuery(centroids, q);
            if (ng != ig) {
                changes[q] = {ig, ng};
            }
        }
    }

    for (const auto& i : changes) {
        auto q = i.first;
        auto& ch = i.second;
        assert(ch[0] >= 0 && ch[0] < groups.size() && ch[1] >=0 && ch[1] < groups.size());
        groups[ch[0]].Remove(q);
        groups[ch[1]].Merge(q);
    }
    return changes.size();
}

std::vector<LocalPhaser::Group> LocalPhaser::BiKMeans(const Group& g0) {

    auto groups = SplitGroups(g0, g0.SelectCentroids());

    for (size_t i = 0; i < 10; ++i) {
        
        if (groups[0].Size() == 0 || groups[1].Size() == 0 || AdjustGroups(groups) == 0) {
            break;
        }
    }

    return groups;
} 
    
std::array<double,2> LocalPhaser::Group::CalcGroupSSE(const PhsOptions::PhaserOptions& opts) const {
    auto centroid = GetCentroid().Values(0.5, 1);
    std::array<double, 2> max_dis = { 0.0, 0.0 };
    std::array<size_t,3> all_stat = {0, 0, 0};
    for (auto q : queries) {
        auto stat = q->Stat(centroid);
        all_stat[0] += stat[0];
        all_stat[1] += stat[1];
        all_stat[2] += stat[2];

        max_dis[0] = std::max<double>(max_dis[0], all_stat[0]);
        max_dis[1] = std::max<double>(max_dis[1], all_stat[0]*1.0 / (all_stat[0]+all_stat[2]));
    }
    DEBUG_local_printf("group sse: (%zd,%zd,%zd) / %zd\n", all_stat[0], all_stat[1], all_stat[2], queries.size());

    return {1 - Query::Similary(all_stat), all_stat[0]*1.0 };
}

std::vector<std::array<int,2>> LocalPhaser::Group::StatVars(const std::unordered_map<int, int>& mapper) const {
    assert(vars.size() > 0);
    
    std::vector<std::array<int,2>> accu_vars(vars.size(), {0,0});

    for (const auto &q : queries) {
        for (const auto &v : q->query->vars) {
            if (v[3] >= 0) {
                auto m = mapper.find(v[0]);
                if (m != mapper.end()) {
                    accu_vars[m->second][v[3]] ++;
                }
            }
        }
    }
    return accu_vars;
}

void LocalPhaser::Run() {
    
    for (auto &q : queries_) {
        PrintGroups("init", {Group(&q)});
    }
    auto groups = Combine3(Divide(Group(queries_)));

    std::vector<size_t> inconsist;
    for (size_t i = 1; i < groups.size(); ++i) {
        if (groups[0].IsInconsistentWithGroup0(groups[i], opts_, coverages_)) {
            FindInconsistent(groups[0], groups[i]);
            inconsist.push_back(i);
        }
    }
    if (inconsist.size() > 0) {
        FindConsistent(groups[0]);
    }

    if (correct_) {
        if (inconsist.size() > 0) {
            CorrectTarget(groups, inconsist);
        } else {
            CorrectTarget(groups);
        }
    }
}


void LocalPhaser::Run2() {
    
    std::vector<Group> groups0;
    for (const auto &q : queries_) {
        groups0.push_back(Group(&q));
    }
        PrintGroups("groups0", groups0);
    auto groups1 = MergeGroups(groups0, 0.8, 0.7, 1);
        PrintGroups("groups1", groups1);
    auto groups2 = MergeGroups(groups1, 0.7, 0.6, 2);
    
        PrintGroups("groups2", groups2);
    auto groups3 = MergeGroups(groups2, 0.6, 0.5, 4);
        PrintGroups("groups3", groups3);

    auto groups = Combine(groups3);
    std::vector<size_t> inconsist;
    for (size_t i = 1; i < groups.size(); ++i) {
        if (groups[0].IsInconsistentWithGroup0(groups[i], opts_, coverages_)) {
            FindInconsistent(groups[0], groups[i]);
            inconsist.push_back(i);
        }
    }
    FindConsistent(groups[0]);

    if (correct_) {
        if (inconsist.size() > 0) {
            CorrectTarget(groups, inconsist);
        } else {
            CorrectTarget(groups);
        }
    }
}

auto LocalPhaser::Divide(const Group& g) -> std::vector<Group> {
    std::vector<Group> groups;
    std::list<Group> work({g});
    while (work.size() > 0) {
        DEBUG_local_printf("working %zd %zd\n", work.size(), groups.size());

        const auto& w = work.front();
        auto sse = w.CalcGroupSSE(opts_);

        DEBUG_local_printf("w(%zd), SSE: %f %f <-> %f %f\n", w.Size(), sse[0], sse[1], opts_.min_inconsist_rate / 3, opts_.min_inconsist_count*1.0 * opts_.min_support_count/2 / 3);
        if ((int)w.Size() <= opts_.min_support_count /2 || (
            sse[1] < opts_.min_inconsist_count*1.0 * opts_.min_support_count/2 / 3 &&
            sse[0] < opts_.min_inconsist_rate / 3)) {
            groups.push_back(w);
        } else {
            auto subs = BiKMeans(w);
            PrintGroups("sub: ", std::vector<Group>(subs.begin(), subs.end()));
            
            if (subs[0].Size() == 0) {
                if (subs[1].Size() != 0)  groups.push_back(subs[1]);
            } else if (subs[1].Size() == 0) {
                if (subs[0].Size() != 0)  groups.push_back(subs[0]);
            } else {
                work.push_back(subs[0]);
                work.push_back(subs[1]);
            }
            if (subs[0].Size() == 0 && subs[1].Size() == 0) {
                groups.push_back(subs[2]);
            } else {
                if (subs[2].Size() > 0) {
        
                    auto ws = subs[2].Split();
                    PrintGroups("ws2", ws);
                    work.insert(work.end(), ws.begin(), ws.end());
                    DEBUG_local_printf("Connected %d\n", subs[2].Connected());
                }
            }
        }
        work.pop_front();
    }
    
    return groups;
}

auto LocalPhaser::Combine(const std::vector<Group> &groups) -> std::vector<Group> {

    struct Item {
        size_t i;
        std::array<int, 3> stat;
    };

    std::vector<Item> items;
    for (size_t i = 0; i < groups.size(); ++i) {
        items.push_back({ i, groups[i].StatDiffs(opts_.min_support_rate, opts_.min_support_count)});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        auto s0 = a.stat[0] + a.stat[2];
        auto s1 = b.stat[0] + b.stat[2];
        if (s0 == 0 || s1 == 0) {
            return s0 > s1;
        } else {
            //return a.stat[0] < b.stat[0] || (a.stat[0] == b.stat[0] && a.stat[2] > b.stat[2]);
            return a.stat[2] - 2*a.stat[0] > b.stat[2] - 2*b.stat[0] ;
        }
    });

    PrintGroups("before combine", groups);
    for (size_t i = 0; i < items.size(); ++i) {
        DEBUG_local_printf("sort -- %zd %d %d %d\n", items[i].i,items[i].stat[0],items[i].stat[1],items[i].stat[2]);
    }

    std::vector<Group> combine ({ groups[items[0].i]});

    for (size_t i = 1; i < items.size(); ++i) {
        auto& it = groups[items[i].i];
        DEBUG_local_printf("test combine\n");
        if (combine[0].IsConsistentWithGroup0(it, opts_)) {
            combine[0].Merge(it);
        } else {
            bool merged = false;
            for (size_t ic = 1; ic < combine.size(); ++ic) {
                if (combine[ic].IsConsistentWith(it, opts_)) {
                    combine[ic].Merge(it);
                    merged = true;
                }
            }
            if (!merged) {
                combine.push_back(it);
            }
        }
    }

    PrintGroups("after combine", combine);
    return combine;
}

auto LocalPhaser::Combine2(const std::vector<Group> &groups) -> std::vector<Group> {

    struct Item {
        size_t i;
        std::array<int, 3> nstat;
        std::array<int, 3> vstat;
    };

    std::vector<Item> items;
    for (size_t i = 0; i < groups.size(); ++i) {
        items.push_back({ i, groups[i].StatDiffs(opts_.min_support_rate),
                             groups[i].StatDiffs(opts_.min_support_rate, opts_.min_support_count)});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.vstat[2] > b.vstat[2] || (a.vstat[2] == b.vstat[2] && a.nstat[2] > b.nstat[2]);
    });

    PrintGroups("before combine", groups);
    for (size_t i = 0; i < items.size(); ++i) {
        DEBUG_local_printf("sort -- %zd %d %d %d\n", items[i].i,items[i].vstat[0],items[i].vstat[1],items[i].vstat[2]);
    }

    size_t n0 = items.size();
    for (size_t i = 0; i < items.size(); ++i) {
        if (groups[items[i].i].queries.find(target_in_queires_) != groups[items[i].i].queries.end()) {
            n0 = i;
            break;
        }
    }

    assert(n0 < items.size());

    std::vector<size_t> nodes { n0 };
    std::vector<std::array<int,2>> dis;
    for (size_t i = 0; i < items.size(); ++i) {
        if (i == n0) continue;
        auto& it = groups[items[i].i];
        auto d1 = groups[items[n0].i].Distance(it, opts_.min_support_rate, opts_.min_support_count);
        auto d2 = groups[items[n0].i].Distance(it, opts_.min_support_rate);
        dis.push_back({d1[0] - d1[2], d2[0]-d2[2]});
        DEBUG_local_printf("sort2 -- %zd %d %d %d, (%d,%d), (%d,%d)\n", items[i].i, items[i].vstat[0],items[i].vstat[1],items[i].vstat[2], d1[0],d1[2], d2[0],d2[2]);

    }


    std::vector<size_t> others;
    if (dis.size() > 0) {

        size_t max_pos = std::max_element(dis.begin(), dis.end(), [](const std::array<int,2> &a, const std::array<int,2>& b) {
            return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
        }) - dis.begin();

        if (max_pos >= n0) {
            max_pos ++;
        }
    
        nodes.push_back(max_pos);

        for (size_t i = 0; i < items.size(); ++i) {
            if (i != nodes[1] && i != nodes[0]) {
                others.push_back(i);
            }
        }

    }
    

    std::vector<Group> combine;
    
    for (auto i : nodes) {
        combine.push_back(groups[items[i].i]);
    } 

    PrintGroups("nodes", combine);
    for (auto i : others) {
        auto& it = groups[items[i].i];

        std::vector<std::array<int, 3>> scores;
        for (size_t ic = 0; ic < combine.size(); ++ic) {
            auto ndis = combine[ic].Distance(it, opts_.min_support_rate);
            auto vdis = combine[ic].Distance(it, opts_.min_support_rate, opts_.min_support_count);
            // auto v2dis = combine[ic].Distance(it, opts_.min_support_rate, opts_.min_support_count/2);
            auto max_distance = combine[ic].MaxDistance(it, opts_.min_support_rate, opts_.min_support_count);
            DEBUG_local_printf("mmm %zd(%zd),%zd(%zd): %d %d %f : %d %d %d---\n", 
                i,it.Size(), ic, combine[ic].Size(), vdis[0], vdis[2], opts_.min_inconsist_rate*(vdis[0] + vdis[2]), 
                ndis[0], ndis[2], max_distance);
            if ((vdis[0] == 0 || vdis[0] < (opts_.min_inconsist_rate*(vdis[0] + vdis[2]) && max_distance < opts_.min_consist_count)
                || vdis[0] < opts_.min_consist_count) &&
                (ndis[0] == 0 || ndis[0] < opts_.min_inconsist_rate*(ndis[0] + ndis[2]))) {
                double s = ndis[0] + ndis[2] > 0 ? (ndis[2] - ndis[0]) *1.0 / (ndis[0]+ndis[2]) : 0;
                scores.push_back({(int)ic, (int)(s*100), ndis[2]});
                DEBUG_local_printf("- mmm %f \n", s);

            }
        }

        std::sort(scores.begin(), scores.end(), [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
            return a[1] > b[1] || (a[1] == b[1]  && a[2] > b[2]) || (a[1] == b[1] && a[2] == b[2] && a[0] < b[0]);
        });

        if (scores.size() > 0 ) {
            combine[scores[0][0]].Merge(it);
        } else {
            combine.push_back(it);
        }
    }

    // first closest groups
    // auto get_score = [this](const Group& g) {
    //     auto s = g.StatDiffs(opts_.min_support_rate, opts_.min_support_count);
    //     return s[2] - s[0]; 
    // };
    // auto best_score = get_score(combine[0]);
    // size_t best = 0;
    // for (size_t i = 1; i < combine.size(); ++i) {
    //     auto s = get_score(combine[i]);
    //     if (s > best_score) {
    //         best = i;
    //         best_score = s;
    //     }
    // }
    // std::swap(combine[0], combine[best]);
    
    for (size_t i = 0; i < combine.size(); ++i) {
        if (combine[i].queries.find(target_in_queires_) != combine[i].queries.end()) {
            if (i != 0) {
                std::swap(combine[0], combine[i]);
            }
            break;
        }
    }


    PrintGroups("after combine", combine);
    
    return combine;
}

auto LocalPhaser::Combine3(const std::vector<Group> &groups) -> std::vector<Group> {

    struct Item {
        size_t i;
        std::array<int, 3> nstat;
        std::array<int, 3> vstat;
        std::array<int, 3> vstat2;
        double Score() const {
            std::array<double, 3> stat;
            stat[0] = vstat[0]*1.0 + vstat2[0]*0.5 + nstat[0]*0.1;
            stat[1] = vstat[1]*1.0 + vstat2[1]*0.5 + nstat[1]*0.1;
            stat[2] = vstat[2]*1.0 + vstat2[2]*0.5 + nstat[2]*0.1;
            return stat[0] + stat[2] == 0 ? 0.0 : (stat[2]-stat[0]) / (stat[0]+stat[2]);
        }
    };

    std::vector<Item> items;
    for (size_t i = 0; i < groups.size(); ++i) {
        items.push_back({ i, groups[i].StatDiffs(opts_.min_support_rate),
                             groups[i].StatDiffs(opts_.min_support_rate, opts_.min_support_count),
                             groups[i].StatDiffs(opts_.min_support_rate, opts_.min_support_count/2)});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.Score() > b.Score();
    });

    PrintGroups("before combine", groups);

    // find the group that contains the target
    size_t n0 = std::find_if(items.begin(), items.end(), [this, &groups](const Item& a) {
        return groups[a.i].queries.find(target_in_queires_) != groups[a.i].queries.end();
    }) - items.begin();
    assert(n0 < items.size());

    DEBUG_local_printf("Node0: %zd (%zd,%zd): %d %d\n", 
                    n0, items[n0].i, groups[items[n0].i].Size(), items[n0].vstat[2] , items[n0].vstat[0]);

    std::vector<Group> combine;
    
    combine.push_back(groups[items[n0].i]);
    std::unordered_set<size_t> c0nodes { items[n0].i };
    for (size_t i = 0; i < items.size(); ++i) {
        auto& it = groups[items[i].i];
        DEBUG_local_printf("CombineToNode0: %zd (%zd,%zd): %d %d\n", 
                        i,items[i].i,it.Size(), items[i].vstat[2] , items[i].vstat[0]);
        if (i == n0) continue;
        if (items[i].vstat[2] > items[i].vstat[0]) {
            
            auto ndis = combine[0].Distance(it, opts_.min_support_rate);
            auto vdis = combine[0].Distance(it, opts_.min_support_rate, opts_.min_support_count);
            auto vdis2 = combine[0].Distance(it, opts_.min_support_rate, opts_.min_support_count/2);
            
            DEBUG_local_printf("cmb0 %zd(%zd): %d %d %f : %d, %d, %d %d\n", 
                            i,it.Size(), vdis[0], vdis[2], opts_.min_inconsist_rate*(vdis[0] + vdis[2]), 
                            vdis2[0], vdis2[2],
                            ndis[0], ndis[2]);

            if ((vdis[0] == 0 || vdis[0] < std::min<double>(opts_.min_inconsist_count, opts_.min_inconsist_rate*(vdis[0] + vdis[2]))) &&
                (vdis2[0] == 0 || vdis2[0] < opts_.min_inconsist_rate*(vdis2[0] + vdis2[2]))) {

                DEBUG_local_printf("Combine0 add %zd\n", i);
                combine[0].Merge(groups[items[i].i]);
                c0nodes.insert(items[i].i);
            }
        
        } else {
            break;
        }
    }

    std::vector<std::array<int,3>> others;
    for (size_t i = 0; i < items.size(); ++i) {
        if (c0nodes.find(items[i].i) != c0nodes.end()) continue;
        auto& it = groups[items[i].i];
        auto d1 = combine[0].Distance(it, opts_.min_support_rate, opts_.min_support_count);
        auto d2 = combine[0].Distance(it, opts_.min_support_rate);
        others.push_back({(int)items[i].i, d1[0], d2[0]});
    }

    std::sort(others.begin(), others.end(), [](const std::array<int,3> &a, const std::array<int,3>& b) {
        return a[1] > b[1] || (a[1] == b[1] && a[2] > b[2]);
    });

    PrintGroups("nodes", combine);
    for (auto &io : others) {
        auto& it = groups[io[0]];

        std::vector<std::array<int, 3>> scores;
        for (size_t ic = 0; ic < combine.size(); ++ic) {
            auto ndis = combine[ic].Distance(it, opts_.min_support_rate);
            auto vdis = combine[ic].Distance(it, opts_.min_support_rate, opts_.min_support_count);
            DEBUG_local_printf("mmm (%zd,%zd),%zd(%zd): %d %d %f : %d %d\n", 
                io[0],it.Size(), ic, combine[ic].Size(), vdis[0], vdis[2], opts_.min_inconsist_rate*(vdis[0] + vdis[2]), 
                ndis[0], ndis[2]);
            if ((vdis[0] == 0 || (vdis[0] < opts_.min_inconsist_rate*(vdis[0] + vdis[2]) && vdis[0] < opts_.min_inconsist_count*2)
                || vdis[0] < opts_.min_consist_count) &&
                (ndis[0] == 0 || ndis[0] < opts_.min_inconsist_rate*(ndis[0] + ndis[2]))) {

                double s = ndis[0] + ndis[2] > 0 ? (ndis[2] - ndis[0]) *1.0 / (ndis[0]+ndis[2]) : 0;
                scores.push_back({(int)ic, (int)(s*10000), ndis[2]});
                DEBUG_local_printf("- mmm %f \n", s);

            }
        }

        std::sort(scores.begin(), scores.end(), [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
            return a[1] > b[1] || (a[1] == b[1]  && a[2] > b[2]) || (a[1] == b[1] && a[2] == b[2] && a[0] < b[0]);
        });

        if (scores.size() > 0 ) {
            combine[scores[0][0]].Merge(it);
        } else {
            combine.push_back(it);
        }
    }

    PrintGroups("after combine", combine);
    
    return combine;
}

void LocalPhaser::FindConsistent(const Group& g0) {
    for (auto q : g0.queries) {
        if (q->query->id == target_->id) continue;

        if (g0.IsConsistentWithGroup0(*q, opts_)) {
                consistent_[ReadOffset::Make(*(target_->o))].push_back(GetMapItem(*target_, *(q->query)));
        }
    }
}


void LocalPhaser::FindInconsistent(const Group& g0, const Group& gx) {
    DEBUG_local_printf("FindInconsistent\n");
    for (auto q : gx.queries) {
        if (q->query->id == target_->id) continue;
        
        if (g0.IsInconsistentWithGroup0(*q, gx, opts_)) {
            inconsistent_[ReadOffset::Make(*(target_->o))].push_back( GetMapItem(*target_, *(q->query)));
        }
    }
}

LocalPhaser::Group::Group(const Query* q) : queries({q}), vars(q->vars.size()) {
    for (size_t i = 0; i < q->vars.size(); ++i) {
        auto v = q->vars[i];
        if (v != 0) {
            vars[i].accu = v;
            vars[i].cov = 1;
        }
    }
    assert(vars.size() > 0);
}

LocalPhaser::Group::Group(const std::unordered_set<const Query*>& qs)
 : queries(qs.begin(), qs.end()), vars((*qs.begin())->vars.size()) {

    for (auto q : qs) {
        for (size_t i = 0; i < q->vars.size(); ++i) {
            auto v = q->vars[i];
            if (v != 0) {
                vars[i].accu += v;
                vars[i].cov += 1;
            }
        }

    }
    assert(vars.size() > 0);
}

LocalPhaser::Group::Group(const std::vector<Query>& qs)
 : vars(qs.begin()->vars.size()) {

    for (const auto &q : qs) {
        queries.insert(&q);
        for (size_t i = 0; i < q.vars.size(); ++i) {
            auto v = q.vars[i];
            if (v != 0) {
                vars[i].accu += v;
                vars[i].cov += 1;
            }
        }

    }
    assert(vars.size() > 0);
}




std::array<int,3> LocalPhaser::Group::StatDiffs() const {
    std::array<int,3> ds = {0, 0, 0};
    std::for_each(vars.begin(), vars.end(), [&ds](const Loc &v) {
        ds[v.GetType()+1] += 1;
    });
    return ds;
}

std::array<int,3> LocalPhaser::Group::StatDiffs(double support_rate) const {
    std::array<int,3> ds = {0, 0, 0};
    std::for_each(vars.begin(), vars.end(), [&ds, support_rate](const Loc &v) {
        ds[v.GetNotableType(support_rate)+1] += 1;
    });
    return ds;
}

std::array<int,3> LocalPhaser::Group::StatDiffs(double support_rate, int support) const {
    std::array<int,3> ds = {0, 0, 0};
    std::for_each(vars.begin(), vars.end(), [&ds, support_rate, support](const Loc &v) {
        ds[v.GetValidType(support_rate, support)+1] += 1;
    });
    return ds;
}

void LocalPhaser::Group::Merge(const Group &b) {
    assert(b.vars.size() > 0);
    if (vars.size() == 0) vars.assign(b.vars.size(), Loc());
    assert(vars.size() == b.vars.size());
    for (size_t i = 0; i < vars.size(); ++i) {
        vars[i].Merge(b.vars[i]);
    }

    queries.insert(b.queries.begin(), b.queries.end());
}

void LocalPhaser::Group::Merge(const Query* q) {
    Merge(Group(q));
}

void LocalPhaser::Group::Remove(const Query* q) {

    for (size_t i = 0; i < vars.size(); ++i) {
        //vars[i].Remove({q->vars[i], 1});
        if (q->vars[i] != 0) {
            vars[i].cov -= 1;
            vars[i].accu -= q->vars[i];
        }
        assert(vars[i].cov >= 0);
    }

    queries.erase(q);
}


bool LocalPhaser::Group::IsConsistentWithGroup0(const Group &b, const PhsOptions::PhaserOptions &opts) const {

    std::array<int,3> count0 = { 0, 0, 0};
    std::array<int,3> count1 = { 0, 0, 0};
    
    for (size_t i = 0; i < vars.size(); ++i) {
        auto &av = vars[i];
        auto &bv = b.vars[i];
        count0[av.GetNotableType(opts.min_support_rate) * bv.GetNotableType(opts.min_support_rate) + 1] ++;
        count1[av.GetValidType(opts) * bv.GetValidType(opts) + 1] ++;
    }

    DEBUG_local_printf("IsConsistentWithGroup0: (%d,%d,%d), (%d,%d,%d)\n", count0[0], count0[1], count0[2], count1[0], count1[1], count1[2]);
    return count0[2] >= std::max<int>((count0[0] + count0[2])*opts.min_consist_rate, opts.min_consist_count) &&
//        count1[0] <= std::min<int>((count1[0] + count1[2])*opts.min_inconsist_rate, opts.min_inconsist_count) &&
        count1[0] <= (count1[0] + count1[2])*opts.min_inconsist_rate &&
        DistanceFromTarget() <= b.DistanceFromTarget();
}

bool LocalPhaser::Group::IsConsistentWith(const Group &b, const PhsOptions::PhaserOptions &opts) const {


    std::array<int,3> count0 = { 0, 0, 0};
    std::array<int,3> count1 = { 0, 0, 0};
    
    for (size_t i = 0; i < vars.size(); ++i) {
        auto &av = vars[i];
        auto &bv = b.vars[i];
        count0[av.GetNotableType(opts.min_support_rate) * bv.GetNotableType(opts.min_support_rate) + 1] ++;
        count1[av.GetValidType(opts) * bv.GetValidType(opts) + 1] ++;
    }

    DEBUG_local_printf("IsConsistentWith: (%d,%d,%d), (%d,%d,%d)\n", count0[0], count0[1], count0[2], count1[0], count1[1], count1[2]);

    return count0[2] >= std::max<int>((count0[0] + count0[2])*opts.min_consist_rate, opts.min_consist_count) &&
//        count1[0] <= std::min<int>((count1[0] + count1[2])*opts.min_inconsist_rate, opts.min_inconsist_count);
        count1[0] <= (count1[0] + count1[2])*opts.min_inconsist_rate;
}

std::array<int, 3> LocalPhaser::Group::Distance(const Group& b, double support_rate) const {

    std::array<int,3> count = { 0, 0, 0};
    
    for (size_t i = 0; i < vars.size(); ++i) {
        auto &av = vars[i];
        auto &bv = b.vars[i];
        count[av.GetNotableType(support_rate) * bv.GetNotableType(support_rate) + 1] ++;
    }
    return count;
}

std::array<int, 3> LocalPhaser::Group::Distance(const Group& b, double support_rate, int support) const {

    std::array<int,3> count = { 0, 0, 0};
    
    for (size_t i = 0; i < vars.size(); ++i) {
        auto &av = vars[i];
        auto &bv = b.vars[i];
        count[av.GetValidType(support_rate, support) * bv.GetValidType(support_rate, support) + 1] ++;
    }
    return count;
}

int LocalPhaser::Group::MaxDistance(const Group& b, double support_rate, int support) const {

    int max_distance = 0;
    int count = 0;
    
    for (size_t i = 0; i < vars.size(); ++i) {
        auto &av = vars[i];
        auto &bv = b.vars[i];
        auto t = av.GetValidType(support_rate, support) * bv.GetValidType(support_rate, support);
        if (t == -1) {
            count ++;
        } else if (t == 1) {
            if (max_distance < count) {
                max_distance = count;
            }
            count = 0;
        }
    }
    return max_distance;
}

bool LocalPhaser::Group::IsInconsistentWithGroup0(const Group &b, const PhsOptions::PhaserOptions &opts, const std::vector<int>& coverage) const {

    auto count = Distance(b, opts.min_support_rate, opts.min_support_count);
    auto max_distance = MaxDistance(b, opts.min_support_rate, opts.min_support_count);
    DEBUG_local_printf("IsInconsistentWithGroup0: (%d, %d, %d, %d)\n", count[0], count[1], count[2], max_distance);
    return count[0] >= opts.min_inconsist_count && (count[0] >= opts.min_inconsist_rate*(count[0] + count[2]));
        //count[0] >= opts.min_inconsist_count*2 ;
}


bool LocalPhaser::Group::IsConsistentWithGroup0(const Query &q, const PhsOptions::PhaserOptions &opts) const {
    assert(vars.size() == q.vars.size());

    std::array<int,3> count = {0, 0, 0};
    for (size_t i = 0; i < q.vars.size(); ++i) {
        auto g0t = vars[i].GetValidType(opts.min_support_rate, opts.min_support_count);
        auto qt = q.vars[i];
        count[1 + g0t * qt] += 1;
    }
    return count[2] >= opts.min_consist_count && 
           count[2] >= opts.min_consist_rate * (count[0] + count[2]);
}

bool LocalPhaser::Group::IsInconsistentWithGroup0(const Query& q, const Group &b, const PhsOptions::PhaserOptions &opts) const {
    assert(vars.size() == q.vars.size());

    std::array<int, 3> count = {0, 0, 0}; // all, diff

    // int distance = 0; distance of two snp

    for (size_t i = 0; i < q.vars.size(); ++i) {
        int at = vars[i].GetValidType(opts);
        int bt = b.vars[i].GetValidType(opts);

        //if (at * bt == -1 && q.vars[i] != 0) {
        if (at * bt != 0 && q.vars[i] != 0) {
            auto t = at * q.vars[i];
            count[t + 1] += 1;

            
        }
    }
    
    DEBUG_local_printf("IsInconsistentWithGroup0(reads): (%d, %d, %d)\n", count[0], count[1], count[2]);
    return count[0] >= opts.min_inconsist_count && (count[0] >= opts.min_inconsist_rate * (count[0] + count[2]));// ||
           //count[0] >= opts.min_inconsist_count*2 ;
}



bool LocalPhaser::Group::Connected() const {
    struct Item {
        const Query* q;
        std::array<size_t, 2> rs;
    };

    std::vector<Item> items;
    for (auto q : queries) {
        std::array<size_t, 2> rs {q->vars.size(), q->vars.size()};
        for (size_t i = 0; i < q->vars.size(); ++i) {
            if (q->vars[i] != 0) {
                if (rs[0] == q->vars.size()) {
                    rs[0] = i;
                }
                rs[1] = i;
            }
        }
        items.push_back({q, rs});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.rs[0] < b.rs[0];
    });

    size_t end = items[0].rs[1];
    std::vector<size_t> split_points;

    for (size_t i = 1; i < items.size(); ++i) {
        if (items[i].rs[0] > end) {
            split_points.push_back(i);
        } else {
            end = std::max(end, items[i].rs[1]);
        }
    }



    return split_points.size() > 0;
}

auto LocalPhaser::Group::Split() const  -> std::vector<Group> {
    
    struct Item {
        const Query* q;
        std::array<size_t, 2> rs;
    };

    std::vector<Item> items;
    for (auto q : queries) {
        std::array<size_t, 2> rs {q->vars.size(), q->vars.size()};
        for (size_t i = 0; i < q->vars.size(); ++i) {
            if (q->vars[i] != 0) {
                if (rs[0] == q->vars.size()) {
                    rs[0] = i;
                }
                rs[1] = i;
            }
        }
        items.push_back({q, rs});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.rs[0] < b.rs[0] || (a.rs[0] == b.rs[0] && a.rs[1] < b.rs[1]);
    });

    size_t end = items[0].rs[1];
    std::vector<size_t> split_points;

    for (size_t i = 1; i < items.size(); ++i) {
        DEBUG_local_printf("items(%zd): %zd-%zd\n", i, items[i].rs[0], items[i].rs[1]);
        if (items[i].rs[0] > end) {
            split_points.push_back(i);
            DEBUG_local_printf("add point\n");
            end = items[i].rs[1];
        } else {
            end = std::max(end, items[i].rs[1]);
        }
    }
    split_points.push_back(items.size());
    DEBUG_local_printf("split_points size: %zd\n", split_points.size());

    std::vector<Group> gs;

    size_t curr = 0;
    for (size_t i = 0; i < split_points.size(); ++i) {
        assert(curr < split_points[i]);
        gs.push_back(Group(items[curr].q));
        curr++;
        while (curr < split_points[i]) {
            gs.back().Merge(items[curr].q);
            curr ++;
        }
    }
    return gs;
}

void LocalPhaser::CorrectTarget(const std::vector<Group>& groups, const std::vector<size_t>& inconsist) {

    std::vector<std::array<size_t,2>> ranges ;

    for (auto i : inconsist) {
        std::vector<Group::Centroid> centroids = { groups[0].GetCentroid(), groups[i].GetCentroid() };
        centroids[0].SetRange2(groups[0], centroids[1], groups[i], opts_);
        DEBUG_local_printf("CorrectTarget: range(%zd): %zd - %zd\n", i, centroids[0].range[0], centroids[0].range[1]);
        ranges.push_back(centroids[0].range);
    }
    std::sort(ranges.begin(), ranges.end(), [](const std::array<size_t, 2>&a, const std::array<size_t, 2>&b) {
        return a[1] - a[0] > b[1] - b[0];
    });

    DEBUG_local_printf("CorrectTarget: final range: %zd - %zd\n", ranges[0][0], ranges[0][1]);

    std::vector<std::vector<std::array<int,2>>> accu_vars;
    accu_vars.push_back(groups[0].StatVars(mapper));
    for (auto i : inconsist) {
        accu_vars.push_back(groups[i].StatVars(mapper));
    }

    std::unordered_set<size_t> valid_position;
    for (const auto &r : ranges) {
        if (r[1] - r[0] >= (size_t)opts_.min_inconsist_count) {
            for (size_t i = r[0]; i < r[1]; ++i) {
                valid_position.insert(i);
            }
        }
        break;
    }

    DEBUG_local_printf("Correct: ");
    auto& tvars = const_cast<std::vector<std::array<int,5>>&>(target_->vars);
    for (size_t i = 0; i < accu_vars[0].size(); ++i) {
        auto& v = accu_vars[0][i];
        auto& tv = tvars[i];

        auto choice = 0;
        //if (i >= ranges[0][0] && i < ranges[0][1]) {
        if (valid_position.find(i) != valid_position.end()) {
            std::vector<int> choices;
            Group::Loc s0(v[0] - v[1], v[0]+v[1]);
            choice = s0.GetValidType(opts_.min_support_rate, opts_.min_support_count);
            choices.push_back(choice);

            for (size_t ix = 1; ix < accu_vars.size(); ++ix) {
                auto& vx = accu_vars[ix][i];
                Group::Loc sx(vx[0] - vx[1], vx[0]+vx[1]);
                auto vchoice = sx.GetValidType(opts_.min_support_rate, opts_.min_support_count);
                DEBUG_local_printf("choice(%d,%d)", choice, vchoice);
                choices.push_back(vchoice);
            }

            bool paired = false;
            for (size_t i = 0; i < choices.size(); ++i) {
                for (size_t j = i+1; j < choices.size(); ++j) {
                    if (choices[i] * choices[j] == -1) {
                        paired = true;
                        break;
                    }
                }
            }
            if (paired) {
                choice = s0.GetValidType(opts_.min_support_rate, opts_.min_support_count/2);
            } else {
                choice = 0;
            }

        }

        if (choice == -1) {
            tv[4] = 1;
        } else if (choice == 1) {
            tv[4] = 0;
        } else {
            if (!opts_.using_vcf) {
                tv[4] = -1;
            }
        }

        DEBUG_local_printf("(%zd|%d,%d,%d->%d), ", i, v[0],v[1], tv[3], tv[4]);
    }
    DEBUG_local_printf("\n");
}

void LocalPhaser::CorrectTarget(const std::vector<Group>& groups) {

    std::vector<std::array<size_t,2>> ranges ;

    for (size_t i = 1; i < groups.size(); ++i) {
        std::vector<Group::Centroid> centroids = { groups[0].GetCentroid(), groups[i].GetCentroid() };
        centroids[0].SetRange2(groups[0], centroids[1], groups[i], opts_);
        DEBUG_local_printf("CorrectTarget: range(%zd): %zd - %zd\n", i, centroids[0].range[0], centroids[0].range[1]);
        ranges.push_back(centroids[0].range);
    }
    std::sort(ranges.begin(), ranges.end(), [](const std::array<size_t, 2>&a, const std::array<size_t, 2>&b) {
        return a[1] - a[0] > b[1] - b[0];
    });

    std::unordered_set<size_t> valid_position;
    for (const auto &r : ranges) {
        if (r[1] - r[0] >= (size_t)opts_.min_inconsist_count) {
            for (size_t i = r[0]; i < r[1]; ++i) {
                valid_position.insert(i);
            }
        }
        break;
    }

    if (ranges.size() > 0) {
        DEBUG_local_printf("CorrectTarget: final range: %zd - %zd\n", ranges[0][0], ranges[0][1]);

        std::vector<std::vector<std::array<int,2>>> accu_vars;
        accu_vars.push_back(groups[0].StatVars(mapper));
        for (size_t i = 1; i < groups.size(); ++i) {
            accu_vars.push_back(groups[i].StatVars(mapper));
        }

        DEBUG_local_printf("Correct0: ");
        auto& tvars = const_cast<std::vector<std::array<int,5>>&>(target_->vars);
        for (size_t i = 0; i < accu_vars[0].size(); ++i) {
            auto& v = accu_vars[0][i];
            auto& tv = tvars[i];

            auto choice = 0;
            //if (i >= ranges[0][0] && i < ranges[0][1]) {
            if (valid_position.find(i) != valid_position.end()) {
                Group::Loc s0(v[0] - v[1], v[0]+v[1]);
                choice = s0.GetValidType(opts_.min_support_rate, opts_.min_support_count);
            }

            if (choice == -1) {
                tv[4] = 1;
            } else if (choice == 1) {
                tv[4] = 0;
            } else {
                if (!opts_.using_vcf) {
                    tv[4] = -1;
                }
            }

            DEBUG_local_printf("(%zd|%d,%d,%d->%d), ", i,v[0],v[1], tv[3], tv[4]);
        }
        DEBUG_local_printf("\n");
    } else {
        // DEBUG_local_printf("ClearTarget\n");
        // ClearTarget();
    }



}

void LocalPhaser::ClearTarget() {
    auto& tvars = const_cast<std::vector<std::array<int,5>>&>(target_->vars);
    for (size_t i = 0; i < tvars.size(); ++i) {
        auto& tv = tvars[i];
        tv[4] = -1;
    }

}

PhaseItem LocalPhaser::GetMapItem(const ReadInfo& target, const ReadInfo& query) const {
    int strand = target.o->SameDirect() == query.o->SameDirect() ? 0 : 1;

    int offset = 0;
    int count = 0;

    size_t it = 0; 
    size_t iq = 0;
    while (it < target.vars.size() && iq < query.vars.size()) {
        int ictg = target.vars[it][0];
        int jctg = query.vars[iq][0];

        if (ictg < jctg) {
            for (; it < target.vars.size(); ++it) {
                if (target.vars[it][0] >= jctg) break;
            }
        } else if (ictg > jctg) {
            for (; iq < query.vars.size(); ++iq) {
                if (query.vars[iq][0] >= ictg) break;
            }
        } else {
            if (target.vars[it][2] != -1 && query.vars[iq][2] != -1) {
                count ++;
                if (strand == 0) {
                    offset += target.vars[it][1] - query.vars[iq][1];
                } else {
                    offset += target.vars[it][1] + query.vars[iq][1];
                }
                //break;
            }
            it++;
            iq++;
        }
    }

    assert(count > 0);
    return {query.o->a_.id, strand, offset / count};

}

void LocalPhaser::PrintGroups(const std::string& msg, const std::vector<Group>& groups) const {
    
    DEBUG_local_printf("------- %s\n", msg.c_str());
    for (size_t i = 0; i < groups.size(); ++i) {
        DEBUG_local_printf("group %zd(%zd): ", i, groups[i].queries.size());
        for (auto q : groups[i].queries) {
            DEBUG_local_printf("%s, ", rd_store_.QueryNameById(q->query->id).c_str());
        }
        //auto diffs = groups[i].StatDiffs(opts_.min_support_rate, opts_.min_support_count);
        auto diffs = groups[i].StatDiffs();
        DEBUG_local_printf("\n    diffs  =(%d,%d,%d):", diffs[0], diffs[1], diffs[2]);
        for (auto v : groups[i].vars) {
            DEBUG_local_printf("(%d,%d),", v.accu, v.cov);
        }
        DEBUG_local_printf("\n");
    }
    fflush(stdout);
}

auto LocalPhaser::MergeGroups(std::vector<Group> &groups, double similary, double coverage, int diff) -> std::vector<Group> {
    std::sort(groups.begin(), groups.end(), [](const Group &a, const Group &b) {
        auto sa = a.StatDiffs();
        auto sb = b.StatDiffs();
        return sa[0]+sa[2] > sb[0]+sb[2];
    });

    std::unordered_set<size_t> done ;
    while (true) {
        size_t donesize = done.size();

        std::vector<size_t> retains;
        for (size_t i = 0; i < groups.size(); ++i) {
            if (done.find(i) == done.end()) {
                auto& qvar = groups[i].vars;

                bool merged = false;
                for (auto ir : retains) {
                    auto tvar = groups[ir].vars;
                    std::array<size_t, 3> nt_stat {0, 0, 0};
                    std::array<size_t, 3> vt_stat {0, 0, 0};
                    std::array<size_t, 3> cov { 0, 0, 0 };
                    for (size_t iv = 0; iv < tvar.size(); ++iv) {
                        int qnt = qvar[iv].GetNotableType(opts_.min_support_rate);
                        int qvt = qvar[iv].GetValidType(opts_.min_support_rate, opts_.min_support_count);
                        int tnt = tvar[iv].GetNotableType(opts_.min_support_rate);
                        int tvt = tvar[iv].GetValidType(opts_.min_support_rate, opts_.min_support_count);
                        if (qnt != 0) cov[0] ++;
                        if (tnt != 0) cov[1] ++;
                        if (qnt != 0 && tnt != 0) cov[2] ++;
                  
                        nt_stat[qnt*tnt + 1] ++;
                        vt_stat[qvt*tvt + 1] ++;
                    }

                    DEBUG_local_printf("test merge: (%zd,%zd,%zd), (%zd,%zd,%zd), (%zd,%zd,%zd)\n",
                        nt_stat[0],nt_stat[1],nt_stat[2], vt_stat[0],vt_stat[1],vt_stat[2], cov[0],cov[0],cov[0]);
                    double qcov = cov[0] > 0 ? cov[2]*1.0 / cov[0] : 0.0;
                    double tcov = cov[1] > 0 ? cov[2]*1.0 / cov[1] : 0.0;
                    double sim = nt_stat[0] + nt_stat[2] > 0 ? nt_stat[2]*1.0 / (nt_stat[0]+nt_stat[2]) : 0.0;
                    if (sim >= similary && (qcov >= coverage || tcov >= coverage ) && nt_stat[0] < (size_t)diff &&
                       (vt_stat[0] == 0 || vt_stat[0] < opts_.min_inconsist_rate*(vt_stat[0] + vt_stat[2]))) {
                        groups[ir].Merge(groups[i]);
                        done.insert(i);
                        merged = true;
                        DEBUG_local_printf("merge: %zd <- %zd\n", ir, i);
                        break;
                    }
                }

                if (!merged) retains.push_back(i);
            }
        }

        //printf("done %zd %zd\n", done.size(), retains.size());
        if (done.size() == donesize) break;
    }
    std::vector<Group> newgroups;
    for (size_t i = 0; i < groups.size(); ++i) {
        if (done.find(i) == done.end()) {
            newgroups.push_back(groups[i]);
        }
    }

    //printf("new %zd %zd\n", groups.size(), newgroups.size());
    DEBUG_local_printf("%zd %zd\n", groups.size(), newgroups.size());
    return newgroups;
}

}   // namespace fsa {
