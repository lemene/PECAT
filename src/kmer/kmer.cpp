#include "kmer.hpp"

namespace fsa {

void KmerSet::BuildIndex() {
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

bool KmerSet1::Find(KmerId kid) const {
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