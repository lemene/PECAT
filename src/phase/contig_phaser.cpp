#include "contig_phaser.hpp"

#include "../read_store.hpp"
#include "../overlap_store.hpp"

#include "phs_options.hpp"
#include "phs_dataset.hpp"
#include "local_phaser.hpp"

namespace fsa {

ContigPhaser::ContigPhaser(Seq::Id ctg, const PhsDataset& dataset, PhsOptions &opts) 
 : ctg_(ctg), dataset_(dataset), rd_store_(dataset_.rd_store_), ol_store_(dataset_.ol_store_), opts_(opts) {

    // 
    Phase();
}

void ContigPhaser::Phase() {

    variants_.assign(GetReadStore().GetSeqLength(ctg_), Variant());

    read_infos_ = CollectReads(ctg_);

    if (opts_.vcf_fname_.empty()) {
        FindVariantsInContig(ctg_, read_infos_, variants_);
        ConfirmVariants(variants_);
        ScanContig();
    } else {
        GetVariantsFromSnps(read_infos_, variants_);
    }
    
    FindVariantsInReads(read_infos_, variants_);
    ClassifyReads(read_infos_);
}

std::unordered_map<ReadOffset, ReadInfo> ContigPhaser::CollectReads(Seq::Id ctg) {
    std::unordered_map<ReadOffset, ReadInfo> read_infos;

    std::unordered_map<Seq::Id, std::vector<const Overlap*>> rgs;
    for (size_t i=0; i < ol_store_.Size(); ++i) {
        const auto &o = ol_store_.Get(i);
//        auto th = opts_.filter_opts_.MaxOverhang(o.a_.len);
        
        if (o.b_.id == ctg) {
//            if (th < 0 || o.a_.start < th || (o.a_.len - o.a_.end) < th) {            
                rgs[o.a_.id].push_back(&o);
//            }
        }
    }

    for (auto& rs : rgs) {
        std::sort(rs.second.begin(), rs.second.end(), [](const Overlap* a, const Overlap *b) {
            double aidt = a->AlignedLength()*a->identity_;
            double bidt = b->AlignedLength()*b->identity_;
            return aidt > bidt || (aidt == bidt && a->b_.start > b->b_.start);
        });
        
        std::vector<const Overlap*> added;
        for (size_t i = 0; i < rs.second.size() && (int)i < opts_.max_count; ++i) {
            const auto o = rs.second[i];

            bool covered = false;
            for (auto ad : added) {
                if (o->b_.id == ad->b_.id) {
                    auto rr = std::min<int>(o->b_.end, ad->b_.end) - std::max<int>(o->b_.start, ad->b_.start);
                    covered = rr > (o->b_.end - o->b_.start) / 2 ;
                    if (covered) break;
                }
            }

            if (!covered) {
                ReadOffset rid = ReadOffset::Make(*o);
                read_infos[rid] = {o->a_.id, o};
                added.push_back(o);
            }
        }
        
    }
    return read_infos;
}


void ContigPhaser::FindVariantsInContig(Seq::Id c, const std::unordered_map<ReadOffset, ReadInfo>& read_infos, std::vector<Variant> &vars) {
    const int C = opts_.snp_match_length_;
    const int variant_type = opts_.variant_type_;
    
    //const auto &ctg = rd_store_.GetSeq(c);  
    for (auto &iter : read_infos) {
        const Overlap& o = *(iter.second.o);
        assert(o.b_.id == c);
        
        const auto &rd = rd_store_.GetSeq(o.a_.id); 
        const auto &ctg = rd_store_.GetSeq(c); 
        size_t ctg_off = 0;
        size_t rd_off = 0;
        for (const auto &d : o.detail_) {
            switch (d.type) {
            case 'M':
                if (variant_type & VARIANT_TYPE_M) {
                    if (d.len >= C) {
                        for (int i=C/2; i<d.len-C/2; ++i) {
                            size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+i : o.b_.end-ctg_off-i-1;
                            size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+i : o.a_.end-rd_off-i-1;

                            //uint8_t ctg_b = ctg[ctg_i];
                            uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : (3 - rd[rd_i]);
                            
                            bool s = true;
                            for (int ii = i - C/2; ii < i+C/2+1; ++ii) {
                                if (ii != i) {
                                    size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+ii : o.b_.end-ctg_off-ii-1;
                                    size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+ii : o.a_.end-rd_off-ii-1;

                                    uint8_t ctg_b = ctg[ctg_i];
                                    uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : (3 - rd[rd_i]);
                                    if (ctg_b != rd_b) {
                                        s = false;
                                        break;
                                    }

                                }

                            }
                            if (s) vars[ctg_i].IncM(rd_b);
                        }
                    }
                }
                ctg_off += d.len;
                rd_off += d.len;
                break;

            case 'D':
                if (variant_type & VARIANT_TYPE_D) {
                    for (int i=0; i<d.len; ++i) {
                        size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+i : o.b_.end-ctg_off-i-1;
                        vars[ctg_i].IncD();
                    }
                }
                ctg_off += d.len;
                break;

            case 'I':
                if (variant_type & VARIANT_TYPE_I) {
                    size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+0 : o.b_.end-ctg_off-0-1;
                    size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+0 : o.a_.end-rd_off-0-1;
                    uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : 3 - rd[rd_i];
                    vars[ctg_i].IncI(rd_b);
                }
                rd_off += d.len;
                break;

            case 'S':
            case 'H':
                break;
            case '=':
            default:
                LOG(ERROR)("Not support cigar type '%c'.", d.type);
            }
        }

        // 计算Coverage
        for (auto i = o.b_.start; i < o.b_.end; ++i) {
            vars[i].IncCov();
        }
        
    }
}


static bool HasHomopolymer(size_t n, const std::vector<int>& bs, size_t start, size_t end) {
    size_t count = 1;
    for (size_t i = start+1; i < end; ++i) {
        if (bs[i] == bs[i-1]) {
            count++;
        } else {
            count = 1;
        }
        if (count >= n) break;
    }
    return count >= n;
}

void ContigPhaser::ConfirmVariants(std::vector<Variant> &vars) {

    const auto &ctg = rd_store_.GetSeq(ctg_);
    assert(ctg.Size() == vars.size());
    // homopolymers
    for (size_t i=0; i<vars.size(); ++i) {
        vars[i].Comfirm(opts_.cov_opts_);

        if (vars[i].Valid()) {
            // 
            if (i < 5 || i + 5 >= vars.size()) {
                vars[i].Disable();
                continue;
            }
            const int N = 7;
            std::vector<int> subctg(N);
            for (size_t j = 0; j < N; ++j) {
                assert(i-N/2+j < ctg.Size());
                subctg[j] = ctg[i-N/2 + j];
            }

            if (HasHomopolymer(4, subctg, 0, subctg.size())) {
                vars[i].Disable();
                continue;
            }
            if (vars[i].var[0] < 4 && vars[i].var[0] != subctg[N/2]) {
                int b = subctg[N/2];
                subctg[N/2] = vars[i].var[0];
                if (HasHomopolymer(4, subctg, 0, subctg.size())) {
                    vars[i].Disable();
                }
                subctg[N/2] = b;
                continue;
            }
            
            if (vars[i].var[1] < 4 && vars[i].var[1] != subctg[N/2]) {
                int b = subctg[N/2];
                subctg[N/2] = vars[i].var[1];
                if (HasHomopolymer(4, subctg, 0, subctg.size())) {
                    vars[i].Disable();
                }
                subctg[N/2] = b;
                continue;
            }
        }
    }

    // std::vector<size_t> density(vars.size(), 0);
    // const size_t R = 10000;
    // auto end = vars.size() >= R ? vars.begin()+R : vars.end();
    // size_t count = std::count_if(vars.begin(), end, [](const Variant &a) {
    //     return a.Valid();
    // });

    // std::fill(density.begin(), density.begin() + (size_t)(end - vars.begin()), count);

    // for (size_t i = 0; i + R < vars.size(); ++i) {
    //     if (vars[i].Valid()) {
    //         count--;
    //         assert(count >= 0);
    //     }
    //     if (vars[i+R].Valid()) {
    //         count++;
    //     }
    //     density[i+R] = count;
    //     for (size_t j = i+1; j < i + R; ++j) {
    //         if (density[j] < count) density[j] = count;
    //     }

    // }
    // const size_t C = 3;
    // for (size_t i = 0; i < vars.size(); ++i) {
    //     if (vars[i].Valid() && density[i] < C) {
    //         vars[i].Disable();
    //     }
    // }
}

void ContigPhaser::GetVariantsFromSnps(const std::unordered_map<ReadOffset, ReadInfo>& read_infos, std::vector<Variant> &vars) {
    auto snp = dataset_.snps_.find(ctg_);

    for (auto &iter : read_infos) {
        const Overlap& o = *(iter.second.o);
        assert(o.b_.id == ctg_);
        
        // 计算Coverage
        for (auto i = o.b_.start; i < o.b_.end; ++i) {
            vars[i].IncCov();
        }
    }

    if (snp != dataset_.snps_.end()) {
        for (const auto& s : snp->second) {
            if (opts_.cov_opts_.Effective(vars[s.first].counts[9])) {
                vars[s.first].var[0] = s.second[0];
                vars[s.first].var[1] = s.second[1];
            }
        }
    }
}


void ContigPhaser::FindVariantsInReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos, const std::vector<Variant> &vars) {
    const int C = opts_.snp_match_length_;
    const int variant_type = opts_.variant_type_;
    for (auto &iter : read_infos) {
        const Overlap& o = *(iter.second.o);

        std::unordered_map<int, std::array<int,5>> cand_vars;

        const auto &rd = rd_store_.GetSeq(o.a_.id); 
        size_t ctg_off = 0;
        size_t rd_off = 0;
        for (const auto &d : o.detail_) {

            switch (d.type) {
            case 'M':
                if (variant_type & VARIANT_TYPE_M) {
                    if (d.len >= C) {
                        for (int i=C/2; i<d.len-C/2; ++i) {
                            size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+i : o.b_.end-ctg_off-i-1;
                            size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+i : o.a_.end-rd_off-i-1;

                            uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : 3 - rd[rd_i];

                            if (vars[ctg_i].Valid()) {
                                if (vars[ctg_i].AtM(rd_b)) {
                                    cand_vars[ctg_i] = {(int)ctg_i, (int)rd_i, rd_b, vars[ctg_i].Offset(rd_b), -1 };
                                    //iter.second.vars.push_back({(int)ctg_i, (int)rd_i, rd_b, vars[ctg_i].Offset(rd_b), -1 });
                                } 
                            }
                        }
                    }
                }
                ctg_off += d.len;
                rd_off += d.len;
                break;

            case 'D':
                if (variant_type & VARIANT_TYPE_D) {
                    for (int i=0; i<d.len; ++i) {
                        size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+i : o.b_.end-ctg_off-i-1;
                        size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+i : o.a_.end-rd_off-i-1;
                        if (vars[ctg_i].AtD()) {
                            //iter.second.vars.push_back({(int)ctg_i, (int)rd_i, (uint8_t)8});
                            cand_vars[ctg_i] = {(int)ctg_i, (int)rd_i, (uint8_t)8, vars[ctg_i].Offset(8), -1};
                        }
                    }
                }
                ctg_off += d.len;
                break;

            case 'I':
                if (variant_type & VARIANT_TYPE_I) {
                    size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+0 : o.b_.end-ctg_off-0-1;
                    size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+0 : o.a_.end-rd_off-0-1;
                    uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : 3 - rd[rd_i];
                    if (vars[ctg_i].AtI(rd_b)) {
                        //iter.second.vars.push_back({(int)ctg_i, (int)rd_i, rd_b+4});
                        cand_vars[ctg_i] = {(int)ctg_i, (int)rd_i, rd_b+4, vars[ctg_i].Offset(rd_b+4), -1};
                    }
                }
                rd_off += d.len;
                break;

            case 'S':
            case 'H':
                break;
            case '=':
            default:
                LOG(ERROR)("Not support cigar type '%c'.", d.type);
            }
            
        }

        for (size_t i = o.b_.start; i < (size_t)o.b_.end; ++i) {
            if (vars[i].Valid()) {
                auto it = cand_vars.find(i);
                if (it != cand_vars.end()) {
                    iter.second.vars.push_back(it->second);
                } else {
                    iter.second.vars.push_back({int(i), -1, -1, -1, -1});
                }
            }
        }
 
    }


}



std::unordered_map<ReadOffset, std::unordered_set<ReadOffset>> ContigPhaser::GroupReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos, int shared_variants) {
    std::vector<const ReadInfo*> sorted_read_infos;
    sorted_read_infos.reserve(read_infos.size());
    for (const auto & i : read_infos) {
        sorted_read_infos.push_back(&i.second);
    }
    
    std::sort(sorted_read_infos.begin(), sorted_read_infos.end(), [](const ReadInfo* a, const ReadInfo* b) {
        return a->o->b_.start < b->o->b_.start || (a->o->b_.start == b->o->b_.start && a->o->b_.end < b->o->b_.end);
    });

    std::unordered_set<ReadOffset> discards;
    std::unordered_map<ReadOffset, std::unordered_set<ReadOffset>> groups;
    for (size_t i=0; i<sorted_read_infos.size(); ++i) {
        
        auto &iri = *sorted_read_infos[i];
        ReadOffset iid = ReadOffset::Make(*sorted_read_infos[i]->o);
        if (discards.find(iid) != discards.end()) continue;   

        for (size_t j=i+1; j < sorted_read_infos.size(); ++j) {

            auto &jri = *sorted_read_infos[j];
            ReadOffset jid = ReadOffset::Make(*sorted_read_infos[j]->o);
            if (discards.find(jid) != discards.end()) continue;   

            if (jri.o->b_.start <= iri.o->b_.end) {
                std::array<int, 2> count {0, 0};
                size_t ii = 0; 
                size_t ji = 0;
                while (ii < iri.vars.size() && ji < jri.vars.size()) {
                    int ictg = iri.vars[ii][0];
                    int jctg = jri.vars[ji][0];

                    if (ictg < jctg) {
                        for (; ii < iri.vars.size(); ++ii) {
                            if (iri.vars[ii][0] >= jctg) break;
                        }
                    } else if (ictg > jctg) {
                        for (; ji < jri.vars.size(); ++ji) {
                            if (jri.vars[ji][0] >= ictg) break;
                        }
                    } else {
                        if (iri.vars[ii][2] != -1 && jri.vars[ji][2] != -1) {
                            count[0] ++;
                        }
                        ii++;
                        ji++;
                    }
                }

                if (count[0] >= shared_variants) {
                    if (dataset_.HasAva() && !dataset_.QueryAva(iri.id, jri.id)) continue;

                        groups[jid].insert(iid);
                        groups[iid].insert(jid);             
                    // if (IsValidPair(*iri.o, *jri.o) && IsValidPair(*jri.o, *iri.o)) {
                    //     groups[jid].insert(iid);
                    //     groups[iid].insert(jid);
                    // } else {
                    //     printf("removej(%zd,%zd) %s-%s\n",i, j, dataset_.QueryNameById(iri.id).c_str()
                    //                                            , dataset_.QueryNameById(jri.id).c_str());
                    //     printf("%s\n%s\n", iri.o->ToM4Line().c_str(), jri.o->ToM4Line().c_str());
                    // }

                }

            }
        }
    }
    return groups;
}

bool ContigPhaser::IsValidPair(const Overlap &irio, const Overlap& jrio) {
    auto ioh = opts_.filter_opts_.MaxOverhang(irio.a_.len);
    if (irio.a_.start > ioh) {
        auto r = irio.MappingToTarget<2>({0, irio.a_.start});
        if (r[0] > r[1]) std::swap(r[0], r[1]);

        auto s = std::max(r[0], jrio.b_.start);
        auto e = std::max(r[1], jrio.b_.end);
        if (e > s && e - s > ioh) {
            return false;
        }
    }

    if (irio.a_.len - irio.a_.end > ioh) {
        auto r = irio.MappingToTarget<2>({irio.a_.end, irio.a_.len});
        if (r[0] > r[1]) std::swap(r[0], r[1]);

        auto s = std::max(r[0], jrio.b_.start);
        auto e = std::max(r[1], jrio.b_.end);
        if (e > s && e - s > ioh) {
            return false;
        }

    }

    return true;
}

void ContigPhaser::ClassifyReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos) {
    auto groups = GroupReads(read_infos, opts_.phase_opts_.min_count); // TODO PHS

    std::vector<std::pair<ReadOffset, int>> read_list(read_infos.size());
    std::transform(read_infos.begin(), read_infos.end(), read_list.begin(), [](const std::pair<ReadOffset, ReadInfo>& a) ->std::pair<ReadOffset, int> {
        return {a.first, a.second.o->AlignedLength()};
    });

    std::sort(read_list.begin(), read_list.end(), [&read_infos](const std::pair<ReadOffset, int> &a,  const std::pair<ReadOffset, int> &b) {
        return a.second > b.second;
    });

    auto add_map_items = [](std::unordered_map<ReadOffset, std::vector<PhaseItem>> &total, const std::unordered_map<ReadOffset, std::vector<PhaseItem>> &sub) {
        for (const auto &i : sub) {
            auto& it = total[i.first];
            it.insert(it.end(), i.second.begin(), i.second.end());
        }
    };

    for (int it = 0; it < opts_.phase_opts_.number_of_iteration; ++it) {
        consistent_.clear();
        inconsistent_.clear();


        for ( const auto& i : read_list) {
            auto &ri = read_infos[i.first];
            for (auto &v : ri.vars) {
                v[4] = v[3];
            }
        }

        std::atomic<size_t> index { 0 };
        std::mutex mutex;
        auto combine_func = [this, &mutex, add_map_items](const std::unordered_map<ReadOffset, std::vector<PhaseItem>>& consistent, 
                                           const std::unordered_map<ReadOffset, std::vector<PhaseItem>>& inconsistent) {
            std::lock_guard<std::mutex> lock(mutex);
            add_map_items(inconsistent_, inconsistent);
            add_map_items(consistent_, consistent);
 
        };

        auto work_func = [this, &index, &groups, &read_list, &read_infos, combine_func, add_map_items, it](size_t thread_id) {
            std::unordered_map<ReadOffset, std::vector<PhaseItem>> consistent;
            std::unordered_map<ReadOffset, std::vector<PhaseItem>> inconsistent;
            std::unordered_set<int> vsnps;
            if (thread_id > 0) opts_.curr_thread_size.fetch_add(1);

            for (size_t i = index.fetch_add(1); i < read_list.size(); i = index.fetch_add(1)) {
                //if (dataset_.QueryNameById(read_list[i].first.id) != "978867" ) continue;
                auto ig = groups.find(read_list[i].first);
           
                if (ig == groups.end() || ig->second.size() == 0) continue;

                LocalPhaser phaser(opts_.phase_opts_, ig->first, ig->second, read_infos, rd_store_, it+1 != opts_.phase_opts_.number_of_iteration);
                phaser.Run();
                add_map_items(consistent, phaser.GetConsistent());
                add_map_items(inconsistent, phaser.GetInconsistent());
                
            }
            combine_func(consistent, inconsistent);
            if (thread_id > 0) opts_.curr_thread_size.fetch_add(-1);

        };

        size_t theads_size = 0;
        AutoThreadRun(work_func, [this, &theads_size]() {
            return opts_.curr_thread_size < (size_t)opts_.thread_size_;
        });

        for ( const auto& i : read_list) {
            auto &ri = read_infos[i.first];
            for (auto &v : ri.vars) {
                v[3] = v[4];
            }
        }
        // for ( const auto& i : inconsistent_) {
        //     auto &ri = read_infos[i.first];
        //     for (auto &v : ri.vars) {
        //         v[3] = v[4];
        //     }
        // }
        
    }
}



void ContigPhaser::ScanContig() {    
    const auto &ctg = rd_store_.GetSeq(ctg_);  
    assert(ctg.Size() == variants_.size());

    auto get_kmers_freq = [](const DnaSeq& c, size_t s, size_t K, size_t W) {
        std::unordered_map<size_t, size_t> kmer_counts;

        size_t kmer = 0;
        for (size_t i = 0; i < K; ++i) {
            kmer += c[s+i] << 2*i;
        }
        kmer_counts[kmer] += 1;
        for (size_t i = 1; i < W; ++i) {
            kmer >>= 2;
            kmer += c[s+i+K-1] << ((K-1)*2);
            kmer_counts[kmer] += 1;
        }

        std::vector<size_t> counts(kmer_counts.size());
        std::transform(kmer_counts.begin(), kmer_counts.end(), counts.begin(), [](const std::pair<size_t, size_t>& a){
            return a.second;
        });

        std::sort(counts.begin(), counts.end(), [](size_t a, size_t b) { return a > b; });
        size_t sum = 0;
        for (size_t i = 0; i < std::min<size_t>(3, counts.size()); ++i) {
            sum += counts[i];
        }

        return sum * 1.0 / W;
    };

    const size_t W = 100;
    const size_t K = 5;
    size_t curr = 0;
    while (curr + W + K < ctg.Size()) {
        auto freq = get_kmers_freq(ctg, curr, K, W);
        if (freq >= 0.66) {
            for (size_t i = 0; i < W; ++i) {
                variants_[curr+i].Disable();
            }
            curr += W;
        } else if (0.66 > freq && freq >= 0.5) {
            curr += 1;
        } else if (0.5 >= freq && freq >= 0.25) {
            curr += W / 2;
        } else {
            curr += W;
        }
    }
}

void ContigPhaser::DumpReadInfos(std::ostream& of) const {
    for (auto &i : read_infos_) {
        auto &ri = i.second;
        assert(ri.o->b_.id == ctg_);

        int offset = ri.o->SameDirect() ? ri.o->a_.start - ri.o->b_.strand : ri.o->a_.start + ri.o->b_.strand;
        for (auto &v : ri.vars) {
            if (v[1] >= 0) {
                offset = ri.o->SameDirect() ? v[0] - v[1] : v[0] + v[1] ;
                break;
            }
        }

        of  << dataset_.QueryNameById(ctg_) << " " << dataset_.QueryNameById(ri.o->a_.id) << " " 
            << (ri.o->SameDirect() ? 0 : 1) << " " << offset << " ";

        for (auto &v : ri.vars) {
            int choice = v[3];
            if (choice == 0 || choice == 1) {
                of << " " <<  v[0] << "|" << int(variants_[v[0]].var[choice]) << "|" << v[1] << "|" << v[2];
            } else {
                of << " " <<  v[0] << "|" << -1 << "|" << v[1] << "|" << v[2];
            }
        }
        of << "\n";
    }
} 

void ContigPhaser::DumpVariants(std::ostream& of) const {
    const auto & cname = dataset_.rd_store_.QueryNameById(ctg_);
    for (size_t i=0; i<variants_.size(); ++i) {
        if (true || variants_[i].Valid()) {
            of << cname << " ";
            of << i;

            for (size_t j=0; j<10; ++j) {
                of << " " << variants_[i].counts[j];
            }
            of << " " << variants_[i].Valid() << "\n";
        }
    }
}
 
void ContigPhaser::DumpInconsistent(std::ostream& of) const {
    for (const auto &i : inconsistent_) {
        of << dataset_.rd_store_.QueryNameById(i.first.id) << ":";
        for (auto p : i.second) {
            if (p.id >= 0) {
                of << dataset_.rd_store_.QueryNameById(p.id) << "|" << p.strand << "|" << p.offset << ",";
            } else {
                of << p.id << "|0|0,";
            }
        }
        of << "\n";
    }
}

void ContigPhaser::DumpConsistent(std::ostream& of) const {
    for (const auto &i : consistent_) {
        if (dataset_.GetOverlapSize(i.first.id) > 1) {

            of << dataset_.rd_store_.QueryNameById(i.first.id) << ":";
            for (auto p : i.second) {
                if (p.id >= 0) {
                    of << dataset_.rd_store_.QueryNameById(p.id) << "|" << p.strand << "|" << p.offset << ",";
                } else {
                    of << p.id << "|0|0,";
                }
            }
            of << "\n";

        }
    }
}
}