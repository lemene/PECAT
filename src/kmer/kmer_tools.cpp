#include "kmer_tools.hpp"

#include "kmer_counter.hpp"
#include "minimizer_counter.hpp"

#include "read_store.hpp"
#include "ranked_kmers.hpp"

namespace fsa {


    
void Program_Count::Running() {
    assert(!ifname_.empty());

}

void Program_Test::Running() {
    assert(!ifname_.empty());
    // Implement the test logic here
    LOG(INFO)("Running test with input file: %s", ifname_.c_str());
    
    //Test_CountingKmer();
    Test_CountingMinimizer();
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



}