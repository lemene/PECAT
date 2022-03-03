#pragma once

#include <vector>
#include <array>
#include <unordered_map>

#include "asm_options.hpp"
#include "../overlap.hpp"

namespace fsa {

    struct RdReason {
        enum Type {
            RS_OK = 0,
            RS_CONTAINED,
            RS_COVERAGE,
            RS_COVERAGE_TYPE,
            RS_OVERHANG,
            RS_UNKNOWN
        };
        bool IsOk() const { return type == RS_OK; }
        static RdReason Ok() { return RdReason(Type::RS_OK); }
        static RdReason Contained(int id) { return RdReason(Type::RS_CONTAINED, id); }
        static RdReason Coverage(int low, int high) { return RdReason(Type::RS_COVERAGE, low, high); }
        static RdReason Coverage(const std::array<int, 2>& a) { return RdReason(Type::RS_COVERAGE, a[0], a[1]); }
        static RdReason CoverageType(int t) { return RdReason(Type::RS_COVERAGE_TYPE, t); }
        static RdReason Overhang() { return RdReason(Type::RS_OVERHANG, 0, 0); }

        static const char* ToString(Type t) {
            switch(t) {
                case RS_OK: return "OK";
                case RS_CONTAINED: return "Contained";
                case RS_COVERAGE: return "Coverage";
                case RS_COVERAGE_TYPE: return "CoverageType";
                case RS_OVERHANG: return "Overhang";
                case RS_UNKNOWN: default: return "Unknown";
            }
        }

        const char* ToString() const { return ToString(type); }
        Type type;
        std::array<int, 2> sub;

protected:
        RdReason(Type t=RS_UNKNOWN, int s0=0, int s1=0) { type = t; sub[0] = s0;  sub[1] = s1;}
    };

    struct OlReason {
        enum Type {
            RS_OK = 0,
            RS_SIMPLE,
            RS_DUPLICATE,
            RS_FILTERED_READ,
            RS_LOCAL,
            RS_CONSISTENCY,
            RS_CONSISTENCY1,
            RS_CONTIG,
            RS_UNKNOWN
        };
        static OlReason Ok() { return OlReason(RS_OK); }
        static OlReason Simple() { return OlReason(RS_SIMPLE); }
        static OlReason Duplicate() { return OlReason(RS_DUPLICATE); }
        static OlReason FilteredRead(int id) { return OlReason(RS_FILTERED_READ, id); }
        static OlReason Consistency(int id) { return OlReason(RS_CONSISTENCY, id); }
        static OlReason Consistency1(int id) { return OlReason(RS_CONSISTENCY1, id); }
        static OlReason Local(int t, int v) { return OlReason(RS_LOCAL, t*10000+v); }
        static OlReason Contig() { return OlReason(RS_CONTIG); }
        
        OlReason(Type t=RS_UNKNOWN, int p0=0, int p1=0) : type(t), sub{p0, p1} { }

        static const char* ToString(Type t) {
            switch(t) {
                case RS_OK:                 return "OK";
                case RS_SIMPLE:             return "Simple";
                case RS_DUPLICATE:          return "Duplicate";
                case RS_FILTERED_READ:      return "FilteredRead";
                case RS_LOCAL:              return "Local";
                case RS_CONSISTENCY:        return "Consistency";
                case RS_CONSISTENCY1:       return "Consistency1";
                case RS_CONTIG:             return "Contig";
                case RS_UNKNOWN: default:   return "Unknonw";
            }

        }
        const char* ToString() { return ToString(type); }

        Type type;
        std::array<int, 2> sub;
    };


    struct IdentityInfo {
        float identity;
        int32_t len;
        int32_t start;
        int32_t end;
    };

    struct OverhangInfo {
        std::array<int16_t, 2> overhang;    // left and right;
        int32_t len;
        int32_t alen;   // aligned len
    };


struct ReadStatInfo {
    int id {-1};
    double overhang_r_threshold {0.0};
    double overhang_l_threshold {0.0};
    std::vector<double> identity_threshold;
    std::array<int, 2>  minmax_coverage {{-1, -1}};
    int len {-1};
    size_t count {0};
    std::array<int, 3> coverage {{0, 0, 0}};
    int covtype { 0 };
    RdReason filtered {RdReason::Ok()};
    std::array<int,2> cliff {{-1, -1}};
    double IdentityThreshold(int start, int end) const ;
    bool CheckIdentity(const Overlap& o) const ;
    bool CheckOverhang(const Overlap& o) const ;
    
    std::vector<IdentityInfo> identities;
    std::vector<OverhangInfo> overhangs;

    void Stat(int id, AsmOptions& opts);
};


class StatReadInfo {
public:
    StatReadInfo(std::unordered_map<Seq::Id, ReadStatInfo> &ri) 
     : read_infos_(ri) {}
    StatReadInfo() 
     : StatReadInfo(default_read_infos_) {}

    size_t Size() const { return read_infos_.size(); }
    void Add(Overlap& o);
    void Merge(StatReadInfo &sri);
    void Clear();
    std::unordered_map<Seq::Id, ReadStatInfo> &read_infos_;
    std::unordered_map<Seq::Id, ReadStatInfo> default_read_infos_;
};


} // namespace fsa