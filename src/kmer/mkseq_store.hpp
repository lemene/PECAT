#pragma once

#include <cassert>
#include <vector>

#include "minimizer_store.hpp"

namespace fsa {

class ReadStore;
class RankedKmers;

class MkseqStore {

public:
    MkseqStore(size_t k=19, size_t w=1000) : k_(k), w_(w) {}
    void Build(const ReadStore& rd_store, RankedKmers* rks=nullptr, size_t threads=1);
    size_t Size() const { return mkseqs_.size(); }
    const MinimizerStore::Seq& Get(size_t i) const { return mkseqs_[i]; }
    const MinimizerStore& GetMinimizerStore() const { return mkmer_store_; }
    void Save(const std::string& fname) const;
    size_t K() const { return k_; }
    size_t W() const { return w_; }
    const Minimizer& GetMinimizer(size_t seq_id, size_t postion) const { 
        assert(seq_id < mkseqs_.size());
        return mkseqs_[seq_id].Get(postion); 
    }
protected:
    size_t k_ { 19 };
    size_t w_ { 1000 };
    const ReadStore* rd_store_;
    std::vector<MinimizerStore::Seq> mkseqs_;
    MinimizerStore mkmer_store_;
};

}