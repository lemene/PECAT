#pragma once


#include <bitset>
#include <fstream>
#include <array>
#include <vector>
#include "align/alignment.hpp"


namespace fsa {


class Corrector {
public:
    static const int MAX_COV = 440;
    
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
        std::bitset<MAX_COV> seqs;
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
    virtual void SetParameter(const std::string &name, int v) = 0;
    virtual void SetParameter(const std::string &name, double v) = 0;

    virtual void Build(const DnaSeq& target, const std::array<size_t,2> &range, const std::vector<Alignment> &aligned) = 0;
    virtual void Clear() = 0;
    virtual void Consensus() = 0;
    virtual const std::string& GetSequence() const = 0;
    virtual const std::string& GetQuality() const = 0;
    virtual void SaveReadInfos(std::ostream &os, int tid, const class ReadStore& rs) const = 0;

 
};


} // namespace fsa {
    