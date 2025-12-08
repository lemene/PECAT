#include "mkseq_store.hpp"

#include <atomic>
#include <mutex>

#include "read_store.hpp"
#include "utility.hpp"
namespace fsa {



void MkseqStore::Build(const ReadStore& rd_store, RankedKmers* rks, size_t thread_size) {
    assert(rd_store_ == nullptr);
    rd_store_ = &rd_store;

    std::mutex mutex;

    mkseqs_.assign(rd_store.Size(), MinimizerStore::Seq());
    LOG(INFO)("read_store %zd %zd", mkseqs_.size(), rd_store.Size());
    auto combine = [&mutex, this](const std::vector<std::vector<Minimizer>>& mkmers, const std::vector<size_t> ids) {
        std::lock_guard<std::mutex> lock(mutex);

        assert(mkmers.size() == ids.size());

        for (size_t i = 0; i < ids.size(); ++i) {
            assert(ids[i] < mkseqs_.size());
            mkseqs_[ids[i]] = mkmer_store_.Add(mkmers[i]);
        }
    };

    
    std::atomic<size_t> index {0};
    auto work_func = [&rd_store, &index, this, &rks, combine](size_t _) {
        thread_local MinimizerCounter mc(k_, w_);
        std::vector<std::vector<Minimizer>> mkmers;
        std::vector<size_t> ids;
        for (size_t i = index.fetch_add(1); i < rd_store.Size(); i = index.fetch_add(1)) {
            if (rd_store.GetSeqLength(i) < w_ + k_) continue; // TODO parameter
            auto kmers = mc.Count(rd_store.GetSeq(i), rks);
            mkmers.push_back(kmers);
            ids.push_back(i);
            
            if (mkmers.size() > 1000) {
                combine(mkmers, ids);
                mkmers.clear();
                ids.clear();
            }
        }
        combine(mkmers, ids);
        LOG(INFO)("combine %zd", mkmer_store_.Size());

    };
    
    MultiThreadRun(thread_size, work_func);
    
    LOG(INFO)("read_store %zd %zd", mkseqs_.size(), rd_store.Size());
    LOG(INFO)("BuildIndex");
    mkmer_store_.BuildIndex();
    LOG(INFO)("END BuildIndex");
}

void MkseqStore::Save(const std::string& fname) const {
    std::ofstream ofile(fname);
    for (size_t i = 0; i < rd_store_->Size(); i++) {
        ofile << rd_store_->QueryNameById(i) << ' ';
        for (size_t j = 0; j < mkseqs_[i].Size(); ++j) {
            auto &m = mkseqs_[i].Get(j);
            ofile << m.pos << '-' << KmerId2String(m.kmer, k_) << ' ';
        }
        ofile << '\n';
    }
}

const std::string& MkseqStore::QueryNameById(size_t seq_id) {
    return rd_store_->QueryNameById(seq_id);
}
}