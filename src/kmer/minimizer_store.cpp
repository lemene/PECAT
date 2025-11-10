#include "minimizer_store.hpp"

#include <numeric>  // iota
namespace fsa {

auto MinimizerStore::Add(const std::vector<Minimizer>& mkmers) -> Seq {
    assert(index_.size() == 0);
    size_t s = mkmers_.size();
    mkmers_.insert(mkmers_.end(), mkmers.begin(), mkmers.end());
    return Seq(this, s, mkmers_.size());
}

void MinimizerStore::BuildIndex() {
    index_.assign(mkmers_.size(), 0);
    std::iota(index_.begin(), index_.end(), 0);
    LOG(INFO)("SS %zd %zd", index_.size(), mkmers_.size());

    std::sort(index_.begin(), index_.end(), [this](size_t a, size_t b) {
        return mkmers_[a].rank <  mkmers_[b].rank || 
              (mkmers_[a].rank == mkmers_[b].rank && mkmers_[a].hash < mkmers_[b].hash);
    });
    LOG(INFO)("SS0");

    std::array<size_t, 2> unique = {0, 0};
    for (size_t i = 1; i < index_.size(); ++i) {
        if (mkmers_[i].hash != mkmers_[unique[0]].hash) {
            unique[1] = i;
            unique_mkmers_[mkmers_[unique[0]].hash] = unique;
            unique[0] = i;
        }
    }
    LOG(INFO)("SS1");
    unique[1] - index_.size();
    unique_mkmers_[mkmers_[unique[0]].hash] = unique;

}

}