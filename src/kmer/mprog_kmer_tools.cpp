#include "mprog_kmer_tools.hpp"

#include "read_store.hpp"

#include "kmer_counter.hpp"
#include "ranked_kmers.hpp"
#include "minimizer_counter.hpp"
#include "minimizer_store.hpp"
#include "minimizer_graph.hpp"
#include "mkseq_store.hpp"


namespace fsa {


    
void Program_Count::Running() {
    assert(!ifname_.empty());

}

template<typename C>
std::vector<std::string> split_string(const std::string &str, C check) {
    auto is_not_split_point = [check](char c) { return !check(c); };
    auto is_split_point = [check](char c) { return check(c); };

    std::vector<std::string> substrs;
    auto begin = std::find_if(str.begin(), str.end(), is_not_split_point);

    while (begin != str.end()) {
        auto end = std::find_if(begin, str.end(), is_split_point);
        substrs.push_back(std::string(begin, end));
        begin = std::find_if(end, str.end(), is_not_split_point);
    }

    return substrs;
}

void Program_Gap::Running() {
    assert(!ifname_.empty());
    std::mutex mutex;
    
    // This function is a placeholder for testing the CountKmer functionality
    LOG(INFO)("Testing CountKmer with k=%d", k_);

    RankedKmers rkmers(split_string(ranked_kmer_fnames_, [](char c) { return c == ';'; }), std::vector<double>(), 20);

    ReadStore rd_store;
    rd_store.Load(ifname_);
    KmerCounter kc(k_);

    std::map<size_t, size_t> all_gaps;
    
    auto combine = [&mutex, &all_gaps](std::unordered_map<size_t, size_t> &gaps) {
        std::lock_guard<std::mutex> lock(mutex);

        for (auto &i : gaps) {
            all_gaps[i.first] += i.second;
        }
    };
 
    std::atomic<size_t> count {0};
    std::atomic<size_t> count1 {0};

    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &kc, &index, &rkmers, combine, &count, &count1, this](size_t _) {
        std::unordered_map<size_t, size_t> gap;
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            auto kmers = kc.CountAll(rd_store.GetSeq(i));
            std::vector<size_t> err_pos(kmers.size());
            for (size_t i = 0; i < kmers.size(); ++i) {
                auto b = (rkmers.GetRank(kmers[i][0] < kmers[i][1] ? kmers[i][0] : kmers[i][1]) >= 0) ? 0 : 1;
                err_pos[i] = b;
            }
            for (size_t i = 0; i + w_ < err_pos.size(); ++i) {
                auto s = std::accumulate(err_pos.begin()+i, err_pos.begin()+i+w_, 0);
                if (s == 0) {
                    count.fetch_add(1);
                } else {
                    count1.fetch_add(1);
                }
            }
        }
        //combine(gap);
    };

    MultiThreadRun(thread_size_, work_func);

    LOG(INFO)("SSS: %zd %zd", count.load(), count1.load());
}



void Program_Bin::Running() {

    auto spec_kmer_fname = SplitStringByChar(specific_, ';');

    // loading specific kmers
    std::vector<KmerStoreUsingVector> spec_kmers;
    for (const auto& fn : spec_kmer_fname) {
        spec_kmers.push_back(KmerStoreUsingVector(fn));
    }

    // check
    for (size_t i = 1; i < spec_kmers.size(); ++i) {
        assert(spec_kmers[i].K() == spec_kmers[i-1].K());
    }
    size_t k = spec_kmers[0].K();

    ReadStore rd_store;
    rd_store.Load(ifname_);

    std::mutex mutex;
    GzFileWriter ofile(ofname_);
    auto save_clear_infos = [&](std::ostringstream& oss) {
        std::lock_guard<std::mutex> lock(mutex);
        ofile << oss.str();
        oss.str("");
    };

    std::atomic<size_t> index { 0 };
    auto work_func = [&](int threadid) {
        std::ostringstream  oss;
        size_t curr = index.fetch_add(1);
        while (curr < rd_store.Size()) {
            const auto& item = rd_store.GetSeq(curr);
            auto count = CountKmers(k, *(item.ToString()), spec_kmers);
            std::vector<double> rate(count.size(), 0);
            for (size_t i = 0; i < rate.size(); ++i) {
                rate[i] = count[i] * 1.0 / spec_kmers[i].Size();
            } 

            std::vector<size_t> idxs(count.size());
            std::iota(idxs.begin(), idxs.end(), 0);

            std::make_heap(idxs.begin(), idxs.end(), [&rate](size_t a, size_t b) { return rate[a] < rate[b]; });
            std::vector<size_t> type = { idxs[0] };
            size_t mx_count = count[idxs[0]];
            double mx_rate = rate[idxs[0]];
            for (size_t ih = 0; ih < idxs.size(); ++ih) {
                std::pop_heap(idxs.begin(), idxs.begin()+idxs.size()-ih);
                if (mx_rate - rate[idxs[0]] > rtol_ && mx_count - count[idxs[0]] > atol_) {
                    break;
                }
                type.push_back(idxs[0]);
            }


            oss << rd_store.QueryNameById(curr) << " " << item.Size() << " ";
            for (size_t c : count) {
                oss << c << " ";
            }
            oss << type[0];
            for (size_t i = 1; i < type.size(); ++i) {
                oss << "_" << type[i];
            }
            oss << "\n";

            if (oss.tellp() >= 10*1024*1024) {
                save_clear_infos(oss);
            }
            curr = index.fetch_add(1);
            if (curr % 10000 == 0) {
                LOG(INFO)("Curr %zd/%zd", curr, rd_store.Size());
            }

        }
        save_clear_infos(oss);
    };

    LOG(INFO)("Classify reads");
    MultiThreadRun((size_t)thread_size_, work_func);
    LOG(INFO)("End classify reads");
}


std::vector<size_t> Program_Bin::CountKmers(size_t k, const std::string& seq, const std::vector<KmerStoreUsingVector>& kmers) {
    std::vector<size_t> count(kmers.size()+1, 0);
    KmerCounter kc (k);
    auto kseq = kc.CountAll(DnaSeq(seq));

    for (size_t i = 0; i < kseq.size(); ++i) {
        auto kmin = std::min(kseq[i][0], kseq[i][1]);
        for (size_t i = 0; i < kmers.size(); ++i) {
            if (kmers[i].Count(kmin) > 0) count[i] ++;
        }
    }
    return count;
}


void Program_Graph::Running() {
    assert(!ifname_.empty());
    std::mutex mutex;
    
    std::unordered_map<KmerId, size_t> all_kmers;

    LOG(INFO)("Loading all reads: %s", ifname_.c_str());
    ReadStore rd_store;
    rd_store.Load(ifname_);

    LOG(INFO)("Loading RankedKmers: %s", ranked_kmer_fnames_.c_str());
    RankedKmers rkmers(SplitStringByChar(ranked_kmer_fnames_, ';'), std::vector<double>(), 20);
    LOG(INFO)("Loaded RankedKmers: size=%zd", rkmers.RankSize());
    
    MkseqStore mkseqs(k_, w_);
    mkseqs.Build(rd_store, &rkmers, thread_size_);
    LOG(INFO)("Minimizer size: %zd, Unique size: %zd", mkseqs.GetMinimizerStore().Size(), mkseqs.GetMinimizerStore().UniqueSize());

    MinimizerGraph graph;
    graph.Build(mkseqs);
    graph.Save("graph.csv");
    mkseqs.Save("mkseqs.txt");

    for (size_t i = 0; i < mkseqs.Size(); ++i) {
        const auto& ks = mkseqs.Get(i);
        std::unordered_set<uint64_t> neighbors;
        for (size_t j = 0; j < ks.Size(); ++j) {
            auto nb0 = graph.Neighbor(ks.Get(j).hash, 0, 3);
            neighbors.insert(nb0.begin(), nb0.end());
            auto nb1 = graph.Neighbor(ks.Get(j).hash, 1, 3);
            neighbors.insert(nb1.begin(), nb1.end());
        }
        LOG(INFO)("get neighbor: size = %zd", neighbors.size());
    }

    //MinimizerStore mkmer_store;

    // auto combine = [&mutex, &mkmer_store](const MinimizerStore& mkmers) {
    //     std::lock_guard<std::mutex> lock(mutex);
    //     mkmer_store.Add(mkmers);
    // };


    // std::atomic<size_t> index {0};
    // auto work_func = [&rd_store, &index, this, combine, &rkmers](size_t _) {
    //     MinimizerCounter mc(k_, w_);
    //     MinimizerStore mkmers;
    //     for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
    //         if (rd_store.GetSeqLength(i) < w_) continue; 
    //         auto kmers = mc.Count(rd_store.GetSeq(i), rkmers);
    //         mkmers.Add(kmers);
    //     }
    //     combine(mkmers);
    // };

    // MultiThreadRun(thread_size_, work_func);
    // mkmer_store.BuildIndex();
 
}

void Program_Histo::Running() {
    assert(!ifname_.empty());
    assert(!ofname_.empty());
    assert(!freq_fname_.empty());

    auto kset = KmerStoreUsingMap(freq_fname_);
    KmerCounter kc(kset.K());

    ReadStore rd_store;
    rd_store.Load(ifname_);

    std::ofstream ofile(ofname_);
    
    
    std::mutex mutex;
    auto save_clear_infos = [&mutex, &ofile](std::ostringstream& oss) {
        std::lock_guard<std::mutex> lock(mutex);
        ofile << oss.str();
        oss.str("");
    };
 

    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &kset, &kc, &index, save_clear_infos, this](size_t _) {
        std::ostringstream oss;
        LOG(INFO)("Thread %zd start %zd", _, rd_store.Size());
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            if (rd_store.GetSeqLength(i) < 1000) continue;
            auto kmers = kc.CountCanon(rd_store.GetSeq(i));

            std::map<KmerId, size_t> kfreq;     // ordered
            for (auto k : kmers) {
                kfreq[kset.Count(k)] ++; // count frequency in the
            }
            oss << rd_store.QueryNameById(i) << " ";
            for (auto &f : kfreq) {
                oss << f.first << "-" << f.second << " ";
            }   
            oss << "\n";
            if (oss.tellp() >= 10*1024*1024) {
                save_clear_infos(oss);
            }
        }
        save_clear_infos(oss);
    };


    MultiThreadRun(thread_size_, work_func);
}

void Program_FreqFreq::Running() {
    assert(!ifname_.empty());

    std::unique_ptr<KmerStoreUsingMap> kset;
    std::unique_ptr<GzFileWriter> local;
    std::unique_ptr<GzFileWriter> global;
    if (!freq_fname_.empty() && !global_fname_.empty()) {
        kset = std::make_unique<KmerStoreUsingMap>(freq_fname_, thread_size_);
        LOG(INFO)("KK: %zd, %zd", kset->K(), k_);
        assert(kset->K() == (size_t)k_);
        global = std::make_unique<GzFileWriter>(global_fname_);
    }

    if (!local_fname_.empty()) {
        local = std::make_unique<GzFileWriter>(local_fname_);
    }

    ReadStore rd_store;
    rd_store.Load(ifname_);
    
    std::mutex mutex;
    auto save_clear_infos = [&mutex, &local, &global](std::ostringstream& ossl, std::ostringstream& ossg) {
        std::lock_guard<std::mutex> lock(mutex);
        
        if (local) {
            local->Write(ossl.str());
        }
        if (global) {
            global->Write(ossg.str());
        }
        ossl.str("");
        ossg.str("");
    };
 
    auto dump_local = [&local, &rd_store](size_t i, const std::vector<KmerId>& kseqs, std::ostream& oss) {
        if (local == nullptr) return;

        std::unordered_map<KmerId, size_t> kcount;
        for (auto k : kseqs) {
            kcount[k]++;
        }

        auto mx = std::max_element(kcount.begin(), kcount.end(), 
            [](const decltype(kcount)::value_type &a,const decltype(kcount)::value_type &b) {
                return a.second < b.second;
        });
        
        std::vector<size_t> histo(mx->second+1, 0);
        for (const auto& it : kcount) {
            histo[it.second] += 1;
        }

        oss << rd_store.QueryNameById(i);
        for (size_t i = 0; i < histo.size(); i++) {
            if (histo[i] > 0) {
                oss  << " " << i << "-" << histo[i];
            }
        }
        oss << "\n";
    };

    auto dump_global = [&global, &rd_store, &kset](size_t i, const std::vector<KmerId>& kmers, std::ostream& oss) {
        if (global == nullptr) return;

        std::map<KmerId, size_t> kfreq;     // ordered
        for (auto k : kmers) {
            kfreq[kset->Count(k)] ++; // count frequency in the
        }

        oss << rd_store.QueryNameById(i) ;
        for (auto &f : kfreq) {
            oss << " " << f.first << "-" << f.second;
        }   
        oss << "\n";
    };

    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &index, save_clear_infos, this, &dump_local, &dump_global](size_t _) {
        LOG(INFO)("Thread %zd start %zd", _, rd_store.Size());
        const size_t MAX_BLCOK_SIZE = 10*1024*1024;
        KmerCounter kc(k_);
        std::ostringstream ossl;
        std::ostringstream ossg;
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            if (rd_store.GetSeqLength(i) < 1000) continue;
            auto kmers = kc.CountCanon(rd_store.GetSeq(i));

            dump_local(i, kmers, ossl);
            dump_global(i, kmers, ossg);

            if (ossl.tellp() >= MAX_BLCOK_SIZE || ossg.tellp() >= MAX_BLCOK_SIZE) {
                save_clear_infos(ossl, ossg);
            }
        }
        save_clear_infos(ossl, ossg);
    };


    MultiThreadRun(thread_size_, work_func);
}



void Program_SegFreq::Running() {
    std::mutex mutex_gen;
    std::mutex mutex_comb;
    
    GzFileReader ifile(ifname_);
    
    size_t k = GetKmerLength(ifname_);
    assert(start_ >= 0 && (size_t)(start_ + len_) <= k && len_  < 6);

    std::vector<std::atomic_uint32_t> seg_freqs(1ULL << (len_*2));
    
    auto work_func = [&mutex_gen, &seg_freqs, &ifile, this](size_t _) {
        LineInBlock line_in_block(ifile, 10000000, &mutex_gen);
    
        KmerCounter kc(len_);
        
        std::string line;
        for (auto valid = line_in_block.GetLine(line); valid; valid = line_in_block.GetLine(line)) {
            auto items = SplitStringBySpace(line);

            auto p = kc.CountForward(items[0].substr(start_, len_))[0];
            seg_freqs[p]++;
        }
    };
    
    MultiThreadRun(thread_size_, work_func);

    std::ofstream ofile(ofname_);
    for (size_t i = 0; i < seg_freqs.size(); ++i) {
        ofile << KmerId2String(i, len_) << "\t" << seg_freqs[i] << "\n";
    }
}

void Program_Verify::Running() {
    assert(!freq_fname0_.empty() && !freq_fname0_.empty() && !ofname_.empty());

    auto kset1 = KmerStoreUsingMap(freq_fname1_);

    std::mutex mutex_gen;
    std::mutex mutex_comb;
    std::ofstream ofile(ofname_);
 
    GzFileReader in(freq_fname0_);

    auto combine_func = [&mutex_comb, &ofile](const std::unordered_map<std::string, size_t>& okmers) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        for (const auto& k : okmers) {
            ofile << k.first << "\t" << k.second << "\n";
        }
    };
        
    auto work_func = [&in, &mutex_gen, combine_func, &kset1](size_t _) {
        LineInBlock line_in_block(in, 1000000, &mutex_gen);
        std::unordered_map<std::string, size_t> okkmers;
        KmerCounter kc(kset1.K());
        std::string line;
        double E = 0.9;
        while (line_in_block.GetLine(line)) {
            auto items = SplitStringBySpace(line);
            size_t freq0 = std::stoul(items[1]);
            for (auto c : {'A', 'C', 'G', 'T'}) {
                auto ks = kc.CountCanon(DnaSeq(items[0] + c)) ;
                assert(ks.size() == 1);
                int freq1 = kset1.Count(ks[0]);
                double pred1 = freq0 * 0.6 + 1;
                if (freq0 > 1)
                LOG(INFO)("%s %c %zd %zd %.02f", items[0].c_str(), c, freq0, freq1, pred1);
                if (freq0 > 100 || freq1 >= pred1) {
                    okkmers[items[0]] = freq0;
                    break;
                }
            }
        }
        combine_func(okkmers);

    };


    MultiThreadRun(thread_size_, work_func);
}



void Program_Test::Running() {
    assert(!ifname_.empty());
    // Implement the test logic here
    LOG(INFO)("Running test with input file: %s", ifname_.c_str());
    
    //Test_CountingKmer();
    //Test_CountingMinimizer();
    //CountMinimizers();
    auto kset = KmerStoreUsingMap(ifname_, thread_size_);
    LOG(INFO)("KmerStore size: %zd", kset.Size());
}

void Program_Test::Test_CountingKmer() {
    // This function is a placeholder for testing the CountKmer functionality
    LOG(INFO)("Testing CountKmer with k=%d", k_);

    ReadStore rd_store;
    rd_store.Load(ifname_);
    KmerCounter kc(k_);
    
    std::atomic<size_t> count {0};

    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &kc, &index, &count](size_t _) {
        std::unordered_map<KmerId, size_t> mkmers;
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            auto kmers = kc.CountAll(rd_store.GetSeq(i));
            count.fetch_add(kmers.size());
            
            for (size_t i = 0; i < kmers.size(); ++i) {
                LOG(INFO)("%zd %08X %08X", i, kmers[i][0],kmers[i][0]);
            }
        }
    };

    MultiThreadRun(thread_size_, work_func);
    LOG(INFO)("Total k-mers counted: %zd", count.load());

}


void Program_Test::Test_CountingMinimizer() {
    // This function is a placeholder for testing the CountKmer functionality
    LOG(INFO)("Testing CountKmer with k=%d", k_);

    RankedKmers rkmers(split_string(ranked_kmer_fnames_, [](char c) { return c == ';'; }), std::vector<double>(), 20);
    RankedKmers rkmers2(std::vector<std::string>(), std::vector<double>(), 20);
    ReadStore rd_store;
    rd_store.Load(ifname_);
    MinimizerCounter mc(k_, w_);
    
    std::atomic<size_t> count {0};
    std::atomic<size_t> count1 { 0 };
    std::atomic<size_t> count2 { 0 };
    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &rkmers, &mc, &index, &count, &count1, &count2, &rkmers, &rkmers2](size_t _) {
        std::unordered_map<KmerId, size_t> mkmers;
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            if (rd_store.GetSeqLength(i) >= 500) {  // TODO parameter
                auto mmer = mc.Count(rd_store.GetSeq(i), rkmers);
                auto mmer2 = mc.Count(rd_store.GetSeq(i), rkmers2);
                std::unordered_set<uint64_t> kkk1;
                for (auto&m : mmer) {
                    kkk1.insert(m.kmer);
                }
                std::unordered_set<uint64_t> kkk2;
                for (auto&m : mmer2) {
                    kkk2.insert(m.kmer);
                }
                size_t c = 0;
                for (auto ik : kkk2)  {
                    if (kkk1.find(ik) != kkk1.end()) {
                        c ++;
                    }
                } 
                count.fetch_add(kkk2.size());
                count1.fetch_add(c);
            }
        }
    };

    MultiThreadRun(thread_size_, work_func);
    LOG(INFO)("Total minimizers counted: %zd %zd %zd", count.load(), count1.load(), count2.load());
}

void Program_Test::CountMinimizers() {
    std::mutex mutex;
    
    std::unordered_map<KmerId, size_t> all_kmers;

    ReadStore rd_store;
    rd_store.Load(ifname_);
    RankedKmers rkmers2(std::vector<std::string>(), std::vector<double>(), 20);
    
    
    auto combine = [&mutex, &all_kmers](const std::unordered_map<KmerId, size_t> mkmers) {
        std::lock_guard<std::mutex> lock(mutex);

        for (auto &i : mkmers) {
            all_kmers[i.first] += i.second;
        }
    };

    std::atomic<size_t> count {0};

    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &index, &count, this, &rkmers2, combine](size_t _) {
        thread_local MinimizerCounter mc(k_, w_);
        std::unordered_map<KmerId, size_t> mkmers;
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            if (rd_store.GetSeqLength(i) < w_) continue; // TODO parameter
            auto kmers = mc.Count(rd_store.GetSeq(i), rkmers2);
            count.fetch_add(kmers.size());
            
            for (size_t i = 0; i < kmers.size(); ++i) {
                mkmers[kmers[i].kmer] += 1;
            }
        }
        combine(mkmers);
    };

    MultiThreadRun(thread_size_, work_func);
    LOG(INFO)("Total minimizers counted: %zd", count.load());
    for (auto &i : all_kmers) {
        printf("Kmer %016llX %zd\n", i.first, i.second);
    }
}


}