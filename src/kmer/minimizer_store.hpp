#pragma once

#include "minimizer_counter.hpp"

namespace fsa {


class MinimizerStore {
public:
    class Seq {
    public:
        Seq(const MinimizerStore* store, size_t s, size_t e) : store_(store), range_({s, e}) {
        }
        Seq() {}
        size_t Size() const { return range_[1]-range_[0]; }
        //const Minimizer& operator [](size_t i) { return store_
        const Minimizer& Get(size_t i) const { 
            assert(i < Size());
            return store_->Get(range_[0] + i); 
        }
    protected:
        const MinimizerStore* store_ {nullptr};
        std::array<size_t, 2> range_ {0, 0};
    };
public:
    MinimizerStore() { }

    size_t Size() const { return mkmers_.size(); }
    size_t UniqueSize() const { return unique_mkmers_.size(); }

    Seq Add(const std::vector<Minimizer>& mkmers);
    void Add(const MinimizerStore& mkmers) { Add(mkmers.mkmers_); }
    void BuildIndex();
    const Minimizer& Get(size_t i) const { return mkmers_[i]; }
protected:
    std::vector<Minimizer> mkmers_;
    std::vector<size_t> index_;
    std::unordered_map<uint64_t, std::array<size_t,2>> unique_mkmers_;
};

}