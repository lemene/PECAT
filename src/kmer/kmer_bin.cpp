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


std::array<size_t, 3> KmerBin::CountKmers(size_t k, const std::string& seq, 
        const KmerSet& patkmers, const KmerSet& matkmers,const KmerSet& offkmers) {
    
    std::array<size_t, 3> count {0, 0, 0};
    KmerCounter kc (k);
    auto kseq = kc.CountAll(DnaSeq(seq));

    // meryl ACTG
    for (size_t i = 0; i < kseq.size(); ++i) {
        auto kid = kseq[i][0];
        auto vkid = kseq[i][1];
        auto kmin = std::min(kid, vkid);
        if (patkmers.Find(kmin)) count[0] ++;
        if (matkmers.Find(kmin)) count[1] ++;
        if (offkmers.Find(kmin)) count[2] ++;
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


} // namespace fsa
