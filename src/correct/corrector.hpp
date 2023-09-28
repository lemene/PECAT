#pragma once


#include <bitset>
#include <fstream>
#include <array>
#include <vector>
#include <unordered_set>
#include "align/alignment.hpp"
#include "../utils/logger.hpp"
#include "overlap_store.hpp"

namespace fsa {

class BitSet {
public:
    static const size_t block_size = 256;

    void reset() { 
        for (auto &b : bits) b.reset();
    }
    void set(size_t i, bool torf) {
        //LOG(INFO)("SSET %d %d", i, torf);
        while (i >= block_size*bits.size()) {
            bits.push_back(std::bitset<block_size>());
        }
        auto ib = i / block_size;
        auto off = i % block_size;
        //LOG(INFO)("off %d %d, %zd", ib, off, bits.size());
        bits[ib].set(off, torf);
        //LOG(INFO)("cSET");
    }

    bool operator [](size_t i) const {
        //LOG(INFO)("S[]S");
        if (i < block_size*bits.size()) {
            auto ib = i / block_size;
            auto off = i % block_size;
            return bits[ib][off];
        } else {
            return false;
        }
        //LOG(INFO)("S[]C");
    }
    BitSet operator &(const BitSet& b) const {
        //LOG(INFO)("SSS");
        BitSet c;
        if (bits.size() >= b.bits.size()) {
            c.bits = b.bits;
            for (size_t i = 0; i < c.bits.size(); ++i) {
                c.bits[i] &= bits[i];
            }
        } else {
            c.bits = bits;
            for (size_t i = 0; i < c.bits.size(); ++i) {
                c.bits[i] &= b.bits[i];
            }
        }
        
        //LOG(INFO)("ccc");
        return c;
    }
    size_t count() const {
        size_t cnt = 0;
        for (const auto &b : bits) {
            cnt += b.count();
        }
        return cnt;
    }
    std::string to_string() const {
        std::string r;
        for (const auto &b : bits) {
            r += b.to_string();
        }
        return r;
    }

    std::vector<std::bitset<block_size>> bits;
};

//typedef std::bitset<500> MyBitSet;
typedef BitSet MyBitSet;

class Corrector {
public:
    
    struct Loc {
        static Loc Invalid() { return { -1, 0, -1}; }

        Loc(int c=-1, int r=-1, int b=-1) : col(c), row(r), base(b) { }
        bool operator == (const Loc &a) const{
            return col == a.col && row == a.row && base == a.base;
        }
        
        bool operator != (const Loc &a) const{
            return !(*this == a);
        }
        bool operator < (const Loc &a) const {
            return col < a.col || (col == a.col && (row < a.row || (row == a.row && base < a.base)));
        }

        int col;      
        int row ;
        int base ;
    };

    struct Link {
        bool operator == (const Link &a) const{
            return prev == a.prev;
        }
        bool operator == (const Loc &a) const{
            return prev == a;
        }
        bool operator < (const Link &a) const {
            return prev < a.prev;
        }

        void Reset(const Loc& p, int id) {
            prev = p;
            count = 1;
            seqs.reset();
            seqs.set(id, true);
        }
        
        Loc prev {-1, -1, -1};
        size_t count {0};
        MyBitSet seqs;
        //double w;       // weight
    };

    struct Tag {
        bool operator == (const Tag &a) const{
            return curr == a.curr && prev == a.prev;
        }
        bool operator != (const Tag &a) const{
            return !(*this == a);
        }
        bool operator < (const Tag &a) const {
            return curr < a.curr || (curr == a.curr && prev < a.prev);
        }
        Loc curr;
        Loc prev;
        int id; // seq id
    };

    
    virtual void SetParameter(const std::string &name, const std::string &v) = 0;
    virtual void SetParameter(const std::string &name, double v) = 0;

    virtual void Build(const DnaSeq& target, const std::array<size_t,2> &range, const std::vector<Alignment> &aligned) = 0;
    virtual void Clear() = 0;
    virtual void Consensus() = 0;
    virtual const std::string& GetSequence() const = 0;
    virtual const std::string& GetQuality() const = 0;
    virtual void SaveReadInfos(std::ostream &os, int tid, const class ReadStore& rs) const = 0;

 
};

class OverlapGrouper {
public:
    OverlapGrouper(OverlapStore& ols) : ol_store_(ols) { }

    void BuildIndex(size_t thread_size, const std::unordered_set<int>& read_ids);

    void ClusterReads(const std::vector<Seq::Id>& reads);

    class Group {
    public:
        Group(Seq::Id i) : id(i) {}
        bool Empty() const { return ols.size() == 0; }
        size_t Size() const { return index.size(); }
        size_t Size(size_t i) const { return index[i][1] - index[i][0]; }
        const Overlap* Get(size_t i, size_t j) { return ols[index[i][0]+j]; }

        void Sort(double opt_ohwt);
        std::vector<double> GetWeight(double opt_ohwt);

        Seq::Id id;
        std::vector<const Overlap*> ols;
        std::vector<std::array<size_t, 2>> index;
    };
    Group Get(int id);

    struct Index {
        std::array<int, 2> by_qurey; 
        std::array<int, 2> by_target;
    };
//protected:
    OverlapStore& ol_store_;
    std::vector<const Overlap*> sorted_;
    std::unordered_map<int, Index> index_;
};


} // namespace fsa {
    