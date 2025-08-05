#pragma once

#include <stdint.h>
#include <vector>
#include <memory>
#include <unordered_set>

#include "../utils/logger.hpp"

#include "../thirdparty/bloom/bloom_filter.hpp"

namespace fsa {
class RankedKmers {
public:
    class Group {
    public:
        Group(const std::string &fname, double level, size_t threads=1);
        bool contains(uint64_t kmer) const { return bloom_->contains(kmer); }
        uint32_t get_kmer_size() const { return kmer_size_; }
    protected:
        bool load(const std::string& fname, size_t threads);
        std::shared_ptr<bloom_filter> make_bloom(size_t count);
        static size_t actual_threads(size_t threads, size_t max_threads)  { return std::min<size_t>(threads, max_threads); }
        uint32_t get_kmer_size_from_file(std::ifstream& ifs) const;
        void load_kmers_to_bloom(std::ifstream& ifs, size_t threads, bloom_filter& bloom);
    protected:
        double level_ { 1.0 };
        std::shared_ptr<bloom_filter> bloom_;
        std::unordered_set<uint64_t> kkk_;
        uint32_t kmer_size_ { 0 };
        size_t count_ { 0 };
    };

public:

    RankedKmers(const std::vector<std::string> &fnames, const std::vector<double> &weights, size_t threads=1);

    size_t RankSize() const { return groups_.size(); }
    int GetRank(uint64_t kmer) const {
        for (size_t i = 0; i < groups_.size(); ++i) {
            if (groups_[i].contains(kmer)) {
                return (int)i;
            }
        }
        return -1;
    }
    
protected:
    bool check() const {
        // TODO 检查kmer size是否匹配。
        for (auto &kg : groups_) {
            if (kg.get_kmer_size() == 0) return false;
        }
        return true;
    }
protected:
    std::vector<Group> groups_;
};


}

