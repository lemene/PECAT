#include "ranked_kmers.hpp"

#include <cassert>
#include <thread>
#include <mutex>
#include <sstream>
#include <fstream>
#include "../utility.hpp"
#include "../utils/logger.hpp"

#include "../file_io.hpp"

namespace fsa {
    
struct BaseTable {
    BaseTable() {
        table['A'] = 0;
        table['C'] = 1;
        table['G'] = 2;
        table['T'] = 3;
        table['a'] = 0;
        table['c'] = 1;
        table['g'] = 2;
        table['t'] = 3;
    }
    uint8_t get(uint8_t c) const { return table[c]; }
    uint8_t table[256];
};



struct Str2Kmer {
public:
    Str2Kmer(uint32_t ksize) {
        shift1 = 2 * (ksize - 1);
        mask = (1ULL<<2*ksize) - 1;
    }

    uint64_t operator()(const std::string &str) {
        uint64_t kmer = 0;
        for (auto a : str) {
            kmer = (kmer << 2 | base_table.get(a)) & mask;           // forward k-mer

        }
        return kmer;
    }

    uint64_t shift1;
    uint64_t mask;
    BaseTable base_table;
};

class BlockReader {
public:
    BlockReader(std::ifstream &ifs, char delim) : ifs_(ifs), delim_(delim) { ifs_.seekg(0); }

    size_t get(char* block, size_t bsize) {
        size_t bindex = 0;
        if (buf_size_ > 0) {
            assert(buf_size_ < bsize);
            std::copy(buf_, buf_+buf_size_, block);
            bindex = buf_size_;
            buf_size_ = 0;
        }

        ifs_.read(block+bindex, bsize-bindex);
        bindex += ifs_.gcount();

        if (bindex == bsize) {
            size_t end = 0;
            for (size_t i = 0; i < bindex; ++i) {
                if (block[bindex-i-1] == delim_) {
                    end = bindex - i; 
                    break;
                }
            }
            //printf("%zd >= %zd && %zd > %zd - %zd\n", bindex, end, max_buf_size_, bindex, end); fflush(stdout);
            assert(bindex >= end && max_buf_size_ > bindex - end);
            std::copy(block+end, block+bindex, buf_);
            buf_size_ = bindex - end;
            bindex = end; 

        } 
        return bindex;
    }

protected:
    const static size_t max_buf_size_ = 1000;
    char buf_[max_buf_size_];
    size_t buf_size_ { 0 };
    std::ifstream &ifs_; 
    char delim_;
};


auto RankedKmers::Group::Create(const std::string &fname, double level, size_t threads) -> Group* {
    std::ifstream ifs(fname);
    auto count = CountLinesInFile(ifs, ActualThreads(threads, 20));
    auto kmer_size = GetKmerSizeFromFile(ifs);
    LOG(INFO)("Loading RankedKmers from %s: size=%zd, k=%d", fname.c_str(), count, kmer_size);

    if (count > 10000000) {
        return new GroupWithBloom(ifs, level, count, kmer_size, threads);
    } else {
        return new GroupWithSet(ifs, level, count, kmer_size, threads);
    }
}

uint32_t RankedKmers::Group::GetKmerSizeFromFile(std::ifstream& ifs) {
    ifs.clear();   ifs.seekg(0);
    std::string kmer;
    uint64_t freq;
    if (ifs >> kmer >> freq) {
        return kmer.size();
    } else {
        return 0;
    }
}

template<typename C>
void RankedKmers::Group::LoadToX(std::ifstream& ifs, size_t threads, C add_to_x) {
    ResetStream(ifs);

    std::mutex mutex_gen;
    std::mutex mutex_comb;

    Str2Kmer str2kmer(kmer_size_);
    BlockReader reader(ifs, '\n');

    auto gen_func = [&reader, &mutex_gen](char* block, size_t bsize) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        return reader.get(block, bsize);
    };

    auto comb_func = [add_to_x, &mutex_comb](const std::vector<uint64_t>& kmers) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        add_to_x(kmers);
    };

    auto work_func = [&gen_func, &comb_func, &str2kmer](size_t id) {
        size_t batch_size = 10000;
        std::vector<uint64_t> batch;

        char block[1*1024*1024];
        size_t bmax = sizeof(block);
        for (size_t bsize = gen_func(block, bmax); bsize > 0; bsize = gen_func(block, bmax)) {
            std::istringstream iss(std::string(block, block+bsize));
            
            std::string kmer;
            uint64_t freq;
            while (iss >> kmer >> freq) {
                batch.push_back(str2kmer(kmer));
            }

            if (batch.size() >= batch_size) {
                comb_func(batch);
                batch.clear();
            }
        }
        comb_func(batch);
        batch.clear();

    };

    MultiThreadRun(ActualThreads(threads, 4), work_func);
}

RankedKmers::GroupWithBloom::GroupWithBloom(std::ifstream& ifs, double level, size_t count, size_t kmer_size, size_t threads) 
 : Group(level, count, kmer_size) {

    Load(ifs, threads);
}

bool RankedKmers::GroupWithBloom::Load(std::ifstream& ifs, size_t threads) {
 
    ResetStream(ifs);

    if (count_ > 0) {
        kmer_size_ = GetKmerSizeFromFile(ifs);
        bloom_ = std::shared_ptr<bloom_filter>(MakeBloom(count_));
        LoadToX(ifs, threads, [this](const std::vector<uint64_t> kmers) {
            bloom_->insert(kmers.begin(), kmers.end());
        });
    } 
}

std::shared_ptr<bloom_filter> RankedKmers::GroupWithBloom::MakeBloom(size_t count) {

    //set up bloom filter
    bloom_parameters parameters;
    parameters.projected_element_count = std::max<size_t>(count, (uint64_t)1000);
    parameters.false_positive_probability = 0.0001; 
    parameters.maximum_number_of_hashes = 2;
    assert(!(!parameters));
    parameters.compute_optimal_parameters();
    return std::shared_ptr<bloom_filter>(new bloom_filter(parameters));
}

RankedKmers::GroupWithSet::GroupWithSet(std::ifstream& ifs, double level, size_t count, size_t kmer_size, size_t threads) 
 : Group(level, count, kmer_size) {
    Load(ifs, threads);
}

bool RankedKmers::GroupWithSet::Load(std::ifstream& ifs, size_t threads) {
    LoadToX(ifs, threads, [this](const std::vector<uint64_t> kmers) {
        kmers_.insert(kmers.begin(), kmers.end());
    });
    return true;
}

RankedKmers::RankedKmers(const std::vector<std::string> &fnames, const std::vector<double> &weights, size_t threads) {
    assert(fnames.size() == weights.size() || weights.size() == 0);

    for (size_t i = 0; i < fnames.size(); ++i) {
        double wt = weights.size() == 0 ? 1.0 : weights[i];
        groups_.push_back(Group::Create(fnames[i], wt, threads));
    }

    assert(Check() && "TODO kmer");
}


}
