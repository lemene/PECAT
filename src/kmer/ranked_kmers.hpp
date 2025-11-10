#pragma once

#include <stdint.h>
#include <vector>
#include <memory>
#include <unordered_set>
#include <fstream>

#include "../utils/logger.hpp"

#include "../thirdparty/bloom/bloom_filter.hpp"

namespace fsa {
class RankedKmers {
public:
    class Group {
    public:
        static Group* Create(const std::string &fname, double level, size_t threads=1);
        virtual bool Contains(uint64_t kmer) const = 0;
        virtual uint32_t GetKmerSize() const { return kmer_size_; }
        virtual ~Group() {}

    protected:
        static uint32_t GetKmerSizeFromFile(std::ifstream& ifs);
        static size_t ActualThreads(size_t threads, size_t max_threads)  { return std::min<size_t>(threads, max_threads); }
        static void ResetStream(std::ifstream& ifs) { ifs.clear(); ifs.seekg(0); }

        template<typename C>
        void LoadToX(std::ifstream& ifs, size_t threads, C add_to_x);
    protected:
        Group(double level, size_t count, size_t kmer_size) : level_(level), count_(count), kmer_size_(kmer_size) {}
        double level_ { 1.0 };
        size_t count_ { 0 };
        uint32_t kmer_size_ { 0 };

    };

    class GroupWithBloom : public Group {
    public:
        GroupWithBloom(std::ifstream& ifs, double level, size_t count, size_t kmer_size, size_t threads=1);
        virtual ~GroupWithBloom() {}
        virtual bool Contains(uint64_t kmer) const { return bloom_->contains(kmer); }
    protected:
        bool Load(std::ifstream& ifs, size_t threads);
        std::shared_ptr<bloom_filter> MakeBloom(size_t count);
    protected:
        std::shared_ptr<bloom_filter> bloom_;
    };

    class GroupWithSet : public Group {
    public:
        GroupWithSet(std::ifstream& ifs, double level, size_t count, size_t kmer_size, size_t threads=1);
        virtual ~GroupWithSet() {}
        virtual bool Contains(uint64_t kmer) const { return kmers_.find(kmer) != kmers_.end(); }
    protected:
        bool Load(std::ifstream& ifs, size_t threads);
    protected:
        std::unordered_set<uint64_t> kmers_;
    };

public:
    RankedKmers(const std::vector<std::string> &fnames, const std::vector<double> &weights, size_t threads=1);
    ~RankedKmers() {
        for (auto g : groups_) {
            delete g;
        }
    }

    size_t RankSize() const { return groups_.size(); }
    int GetRank(uint64_t kmer) const {
        for (size_t i = 0; i < groups_.size(); ++i) {
            if (groups_[i]->Contains(kmer)) {
                return (int)i;
            }
        }
        return -1;
    }
    
protected:
    bool Check() const {
        // TODO 检查kmer size是否匹配。
        for (auto &g : groups_) {
            if (g->GetKmerSize() == 0) return false;
        }
        return true;
    }
protected:
    std::vector<Group*> groups_;
};


}

