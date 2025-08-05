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


size_t count_lines_in_file(std::ifstream& ifs, size_t threads) {
    ifs.seekg(0);
    std::string line;
    std::mutex mutex;
    std::atomic<size_t> count { 0 };

    auto generate_func = [&ifs, &mutex](char* buf, size_t bufsize) -> size_t {
        std::lock_guard<std::mutex> lock(mutex);

        ifs.read(buf, bufsize);
        return ifs.gcount();
    };
    auto worker_func = [&count, &generate_func](size_t id) {
        char buf[1024*1000];
        size_t bsize = generate_func(buf, sizeof(buf));
        size_t cnt = 0;
        while (bsize > 0) {
            for (char* p = buf; p < buf+bsize; ++p) {
                if (*p == '\n') cnt ++;
            }
            bsize = generate_func(buf, sizeof(buf));
        }
        count.fetch_add(cnt);
    };

    MultiThreadRun(threads, worker_func);
    return count.load();
}

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


RankedKmers::Group::Group(const std::string &fname, double level, size_t threads) {
    level_ = level;
    load(fname, threads);
}

bool RankedKmers::Group::load(const std::string& fname, size_t threads) {
 
    std::ifstream ifs(fname);

    count_ = count_lines_in_file(ifs, actual_threads(threads, 20));
    if (count_ > 0) {
        kmer_size_ = get_kmer_size_from_file(ifs);
        LOG(INFO)("%zd kmers(k=%u) in file %s", count_, kmer_size_, fname.c_str());
        bloom_ = std::shared_ptr<bloom_filter>(make_bloom(count_));
        load_kmers_to_bloom(ifs, threads, *bloom_.get());
    } else {
        LOG(INFO)("Empty file: %s", fname.c_str());
    }
	LOG(INFO)("[M::%s] load kmers %s", __func__, fname.c_str());
}

std::shared_ptr<bloom_filter> RankedKmers::Group::make_bloom(size_t count) {

    //set up bloom filter
    bloom_parameters parameters;
    parameters.projected_element_count = std::max<size_t>(count, (uint64_t)1000);
    parameters.false_positive_probability = 0.001; 
    parameters.maximum_number_of_hashes = 2;
    assert(!(!parameters));
    parameters.compute_optimal_parameters();
    return std::shared_ptr<bloom_filter>(new bloom_filter(parameters));
}

uint32_t RankedKmers::Group::get_kmer_size_from_file(std::ifstream& ifs)  const {
    ifs.clear();   ifs.seekg(0);
    std::string kmer;
    uint64_t freq;
    if (ifs >> kmer >> freq) {
        return kmer.size();
    } else {
        return 0;
    }
}

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

    
void RankedKmers::Group::load_kmers_to_bloom(std::ifstream& ifs, size_t threads, bloom_filter& bloom) {
    ifs.clear();  ifs.seekg(0);

    std::mutex mutex_gen;
    std::mutex mutex_comb;

    assert(kmer_size_ > 0);
    Str2Kmer str2kmer(kmer_size_);
    BlockReader reader(ifs, '\n');

    auto gen_func = [&reader, &mutex_gen](char* block, size_t bsize) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        return reader.get(block, bsize);
    };

    auto comb_func = [&bloom, &mutex_comb](const std::vector<uint64_t>& kmers) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        bloom.insert(kmers.begin(), kmers.end());
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

    MultiThreadRun(actual_threads(threads, 4), work_func);
}

RankedKmers::RankedKmers(const std::vector<std::string> &fnames, const std::vector<double> &weights, size_t threads) {
    assert(fnames.size() == weights.size() || weights.size() == 0);

    for (size_t i = 0; i < fnames.size(); ++i) {
        
        double wt = weights.size() == 0 ? 1.0 : weights[i];
        groups_.push_back(Group(fnames[i], wt, threads));
    }

    assert(check() && "TODO kmer");
}


}
