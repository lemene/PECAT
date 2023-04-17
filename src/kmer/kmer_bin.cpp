#include "kmer_bin.hpp"

#include <cassert>
#include <deque>

#include "../file_io.hpp"
#include "../utility.hpp"
#include "../read_store.hpp"

namespace fsa {

ArgumentParser KmerBin::GetArgumentParser() {
    ArgumentParser ap("fsa_kmer_tools", "tools about reads", "1.0");
    ap.AddNamedOption(maternal_, "maternal", "kmers in maternal haplotype");
    ap.AddNamedOption(paternal_ , "paternal", "kmers in paternal haplotype");
    ap.AddNamedOption(offspring_ , "offspring", "kmers in paternal haplotype");
    ap.AddNamedOption(thread_size_, "thread_size", "number of threads");
    ap.AddNamedOption(ifname_, "ifname", "input file name");
    ap.AddNamedOption(ofname_, "ofname", "output file name");
    ap.AddNamedOption(output_directory_, "output_directory", "output directory");
    ap.AddNamedOption(th_count_, "count", "");
    ap.AddNamedOption(th_rate_, "rate", "");
    return ap;
}

void KmerBin::Running() {

    auto patkmers = LoadKmers1(paternal_);
    auto matkmers = LoadKmers1(maternal_);
    auto offkmers = LoadKmers1(offspring_);

    size_t k = CheckKmerSet(patkmers, matkmers, offkmers);

    ReadStore rd_store;
    rd_store.Load(ifname_);

    struct InfoItem {
        Seq::Id id;
        size_t len;
        size_t pat_only;
        size_t mat_only;
        size_t found;       // 找到的kmer
        int type; // 1 paternal. -1 maternal 0 ambigours
    };

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
        thread_local std::vector<InfoItem> infos;
        size_t curr = index.fetch_add(1);
        while (curr < rd_store.Size()) {
            const auto& item = rd_store.GetSeq(curr);
            auto count = CountKmers(k, *(item.ToString()), patkmers, matkmers, offkmers);
            std::array<double, 2> thresholds { std::max<double>(th_count_, count[0]*th_rate_), std::max<double>(th_count_, count[1]*th_rate_)};
            
            InfoItem info;
            info.id = curr;
            info.len = item.Size();
            info.pat_only = count[0];
            info.mat_only = count[1];
            info.found = count[2];
            info.type = count[0]*1.0 / patkmers.Size() > (count[1]*1.0 + thresholds[1]) / matkmers.Size() ? 1 :
                        (count[0]*1.0+thresholds[0]) / patkmers.Size() < count[1]*1.0 / matkmers.Size() ? -1 : 0;

            oss << rd_store.QueryNameById(curr) << " " 
                  << info.len << " "
                  << info.pat_only << " "
                  << info.mat_only << " "
                  << info.found << " " 
                  << info.len - k + 1 - info.found << " "
                  << info.type << "\n";

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

KmerBin::KmerSet0 KmerBin::LoadKmers0(const std::string &fname) {
    KmerSet0 kmers;
    std::mutex mutex_gen;
    std::mutex mutex_comb;

    kmers.k = GetKmerLength(fname);

    const size_t block_size = 1000;
    GzFileReader reader(fname);
    auto generate_func = [&mutex_gen, &reader](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        return reader.GetLines(lines);
    };

    auto combine_func = [&mutex_comb, &kmers](std::unordered_map<KmerId, int>& ks) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        kmers.kmers.insert(ks.begin(), ks.end());
        ks.clear();
    };

    auto work_func = [block_size, generate_func, combine_func, &reader, this](size_t id) {
        std::vector<std::string> lines(block_size);
        std::unordered_map<KmerId, int> ks;

        size_t sz = generate_func(lines);
        while (sz > 0) {
            
            for (size_t i=0; i<sz; ++i) {
                auto items = SplitStringBySpace(lines[i]);
                ks[KmerStringToId(items[0])] = std::stoi(items[1]);
            }
            
            combine_func(ks);
            sz = generate_func(lines);
        }
    };
    
    MultiThreadRun(std::min<size_t>(thread_size_, 3), work_func);

    LOG(INFO)("Load %zd kmers(k=%zd) from %s", kmers.kmers.size(), kmers.k, fname.c_str());
    return kmers;
}

KmerBin::KmerSet KmerBin::LoadKmers1(const std::string &fname) {
    KmerSet kmers;

    kmers.k = GetKmerLength(fname);

    GzFileReader reader(fname);
    std::string line;
    while (reader.GetLine(line)) {
        auto items = SplitStringBySpace(line);
        kmers.kmers.push_back({KmerStringToId(items[0]),std::stoi(items[1])});
    }

    LOG(INFO)("Load %zd kmers(k=%zd) from %s", kmers.kmers.size(), kmers.k, fname.c_str());
    std::sort(kmers.kmers.begin(), kmers.kmers.end(), [](const KmerItem& a, const KmerItem &b) { return a.kmer < b.kmer; });
    LOG(INFO)("Sort %zd kmers(k=%zd) from %s", kmers.kmers.size(), kmers.k, fname.c_str());
    kmers.BuildIndex();
    return kmers;
}

std::array<size_t, 3> KmerBin::CountKmers(size_t k, const std::string& seq, 
        const KmerSet& patkmers, const KmerSet& matkmers,const KmerSet& offkmers) {
    
    std::array<size_t, 3> count {0, 0, 0};
    for (size_t i=0; i+k <= seq.size(); ++i) {
        std::string s(seq.begin()+i, seq.begin()+i+k);
        std::for_each(s.begin(), s.end(), [](char &c) { c = ::toupper(c); });
        std::string s_rc = Seq::ReverseComplement(s);

        auto kid = KmerStringToId(s);
        auto vkid = KmerStringToId(s_rc);


        if (patkmers.Find(kid) || patkmers.Find(vkid)) count[0] ++;
        if (matkmers.Find(kid) || matkmers.Find(vkid)) count[1] ++;
        if (offkmers.Find(kid) || offkmers.Find(vkid)) count[2] ++;
        if (thread_size_ == 1) {
            if (patkmers.Find(kid)) printf("%zd 1, %s\n", i, s.c_str());
            if (patkmers.Find(vkid)) printf("%zd 1, %s\n", i, s_rc.c_str());
            if (matkmers.Find(kid)) printf("%zd -1, %s\n", i, s.c_str());
            if (matkmers.Find(vkid)) printf("%zd -1, %s\n", i, s_rc.c_str());
        }
    }
    return count;
}
size_t KmerBin::CheckKmerSet(const KmerSet& patkmers, const KmerSet& matkmers, const KmerSet& offkmers) const {

    size_t k = 0;
    if (!offkmers.Empty()) {
        k = offkmers.k;
        if (!patkmers.Empty() && k != patkmers.k) {
            LOG(ERROR)("k in pateral is not equal the one in offspring");
        }
        if (!matkmers.Empty() && k != matkmers.k) {
            LOG(ERROR)("k in materal is not equal the one in offspring");
        }
    } else {
        if (!patkmers.Empty() && !matkmers.Empty()) {
            k = patkmers.k;
            if (patkmers.k !=  matkmers.k) {
                LOG(ERROR)("k in materal is not equal the one in pateral");
            }
        } else {
            LOG(ERROR)("pateral or maternal is empty");
        }
    }
    return k;
}

KmerBin::KmerId KmerBin::KmerStringToId(const std::string& str) {
    static DnaSerialTable table;
    KmerId id = 0;
    for (auto c : str) {
        id = (id << 2) + table[c];
    }
    return id;    
}

size_t KmerBin::GetKmerLength(const std::string &fname) {
    GzFileReader reader(fname);
    if (reader.Valid()) {
        auto line = reader.GetNoEmptyLine();
        if (!line.empty()) {
            auto items = SplitStringBySpace(line);
            return items[0].size();

        }

    }
    return 0;
}

void KmerBin::KmerSet::BuildIndex() {
    index.assign(1024, {kmers.size(),0});

    for (size_t i = 0; i<kmers.size(); ++i) {
        size_t idx = kmers[i].kmer >> (k*2 - 10);

        if (index[idx][0] > i) {
            index[idx][0] = i;
        }

        if (index[idx][1] < i+1) {
            index[idx][1] = i+1;
        }
    }

}

bool KmerBin::KmerSet1::Find(KmerId kid) const {
    auto se = index[kid >> (k*2 - 10)];
    size_t s = se[0]; 
    size_t e = se[1];
    //printf("s e %zd %zd\n", s, e);

    while (s < e) {
        size_t m = (s+e) / 2;
        if (kmers[m].kmer == kid) {
            return true;
        } else if (kmers[m].kmer < kid) {
            s = m+1;
        } else {
            e = m;
        }
    }
    return false;
}
} // namespace fsa
