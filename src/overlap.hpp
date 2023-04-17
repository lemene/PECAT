#ifndef FSA_OVERLAP_HPP
#define FSA_OVERLAP_HPP

#include <string>
#include <array>

#include "sequence.hpp"

namespace fsa {

class Overlap {
public:
    struct Filter {
        Filter() {}
        Filter(const std::string& strs) { From(strs); }
        void From(const std::string& strs);
        std::string ToString() const ;
        bool ValidQuery(const Overlap &ol) const;
        bool Valid(const Overlap &ol, int replen=0, bool bothside=true) const;
        double Identity(const Overlap &ol) const { return large_indel_ > 0 ? IdentityIgnoringIndel(ol, large_indel_) : ol.identity_; }
        double IdentityIgnoringIndel(const Overlap &ol, int indel) const;
        int MaxOverhang(int len) const;

 
        int min_length { 0 };                   // l
        double min_identity { 0 };              // i
        int min_aligned_length { 0 };          // alr
        double min_aligned_rate { 1.0 };          // al
        int max_overhang { 10000000 };           // oh          
        double max_overhang_rate { 1.0 };       // ohr
        int min_accept_aligned_length { 10000000 };       // aal
        double min_accept_aligned_rate { 1.0 };       // aalr
        int large_indel_ { 0 };                           // sid
    };

public:
	struct Read {
		int id;
		int strand;
		int start;
		int end;
		int len;
	};

    struct Detail {
        int len : 24;
        char type;
    };
public:


    enum class Loc {
        Left, Right, Contained, Equal, Containing, Abnormal
    };
public:
    Overlap() = default;

    std::string ToM4Line() const;
    bool FromM4Line(const std::string &line);
    
    bool CheckEnd(int error) const;

    bool SameDirect() const { return a_.strand == b_.strand; }
    size_t AlignedLength() const { return (a_.end - a_.start + b_.end - b_.start) / 2;  }
    
    size_t IdealAlignedLength() const { 
        if (SameDirect()) {
            return std::min(a_.start, b_.start) + AlignedLength() + std::min(a_.len-a_.end, b_.len-b_.end);
        } else {
            return std::min(a_.start, b_.len-b_.end) + AlignedLength() + std::min(a_.len-a_.end, b_.start);
        }
    }

    bool IsContaining(int err) const;
    bool IsContained(int err) const;
    bool IsAbnormal(int err) const;
    bool IsProper(int err) const;
    
    Loc Location(int err) const;
    Loc Location(Seq::Id id, int err) const ;
    Loc Location(const Read& r, int err) const { return Location(r.id, err);}

    static std::array<int, 2> Overhang(const Overlap &o, Loc loc);

    std::array<int, 2> Overhang() const;
    std::array<int, 2> Overhang2() const;

    template<int N>
    std::array<int, N>  MappingToSource(const std::array<int, N>& pos) const {
        return Mapping<N>(a_, b_, pos);
    }

    template<int N>
    std::array<int, N>  MappingToTarget(const std::array<int, N>& pos) const {
        return Mapping<N>(b_, a_, pos);
    }

    template<int N>
    std::array<int, N> MappingTo(const Read &a, const std::array<int, N> &pos) const {
        if (&a == &a_) { 
            return MappingToSource<N>(pos);
        } else if (&a == &b_) {
            return MappingToTarget<N>(pos);
        } else {
            assert("a is not valid parameter");
            return std::array<int, N>();
        }
    }

    template<int N>
    static std::array<int, N> Mapping(const Read& a, const Read& b, const std::array<int, N> &bpos);
    
    static bool IsConsistent(const Overlap &ab, const Overlap &ac, const Overlap &bc, int err);

    static int ReverseStrand(int s) { return s == 0 ? 1 : 0;  }
    static Loc ReverseLocation(Loc loc, bool direct=true);

    
    bool Extend(int overhang);
    
    const Read& GetRead(Seq::Id id) const { 
        assert(id == a_.id || id == b_.id); 
        return a_.id == id ? a_ : b_;
    }

    const Read& GetOtherRead(Seq::Id id) const {
        assert(id == a_.id || id == b_.id); 
        return a_.id != id ? a_ : b_;
    }


    int Extension(Seq::Id id, int end) const;

	Read a_;    // query
	Read b_;    // target
	double identity_;
    std::vector<Detail> detail_;
    mutable long long attached {0};   // used by other, and avoid building map (Overlap -> info)
};


template<int N>
std::array<int, N> Overlap::Mapping(const Read& a, const Read& b, const std::array<int, N> &bpos) {
    std::array<int, N> apos;
    if (a.strand == b.strand) {
        for (size_t i=0; i<N; ++i) {
            if (bpos[i] < b.start) {
                apos[i] = a.start + (bpos[i] - b.start);
            } else if (bpos[i] >= b.start && bpos[i] < b.end) {
                apos[i] = a.start + (long long)(bpos[i] - b.start) * (a.end - a.start) / (b.end - b.start);
            } else {
                assert(bpos[i] >= b.end);
                apos[i] = a.end + (bpos[i] - b.end);
            }
        }
    }
    else {
        for (size_t i=0; i<N; ++i) {
            if (bpos[i] > b.end) {
                apos[i] = a.start - (bpos[i] - b.end);
            } else if (bpos[i] <= b.end && bpos[i] > b.start) {
                apos[i] = a.start - (long long)(bpos[i]-b.end) * (a.end - a.start) / (b.end - b.start);
            } else { 
                assert(bpos[i] <= b.start);
                apos[i] = a.end - (bpos[i] - b.start);
            }
        }
    }
    return apos;
}

} // namespace fsa {

#endif  // FSA_OVERLAP_HPP
