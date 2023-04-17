#pragma once

#include <array>
#include <vector>

#include "../sequence.hpp"
#include "../overlap.hpp"
#include "../utils/misc.hpp"

namespace fsa {

struct CoverageOptions : public StringOptions {
    CoverageOptions(const std::string& strs) { From(strs); }
    void From(const std::string& strs);
    std::string ToString() const;


    bool Effective(int c) const { return c >= valid_range[0] && c <= valid_range[1]; }
    int Threshold(int c) const {
        double slope = (rate[1]*range[1] - rate[0]*range[0]) / (range[1] - range[0]);
        if (c < range[0])       return c* rate[0];
        else if (c < range[1])  return (c - range[0]) * slope + rate[0]*range[0];
        else                    return c * rate[1];
    }

    std::array<int,2> range {{10, 100}};
    std::array<double,2> rate {{0.4, 0.20}};
    std::array<int, 2> valid_range {{10, 1000}};
    double slope;
};

struct ContigInfo {

};

struct Variant {
    void IncM(uint8_t v) { counts[v]++; }
    void IncD()          { counts[8]++; }
    void IncI(uint8_t v) { counts[v+4]++;}
    void IncCov()        { counts[9]++; }
    bool Valid() const  { return var[0] != var[1]; }
    bool AtM(uint8_t v) const { return v == var[0] || v == var[1]; }
    bool AtD() const { return var[0] == 8 || var[1] == 8; }
    bool AtI(uint8_t v) const { return v+4 == var[0] || v+4 == var[1];}
    void Comfirm(const CoverageOptions& covopts);
    void Disable() { var[0] = 0; var[1] = 0; }
    int Offset(uint8_t b) const {
        return var[0] == b ? 0 : (var[1] == b ? 1 : -1);
    }
    
    uint16_t counts[10] = {0};
    int8_t var[2] = {0};
};

struct ReadOffset {
    ReadOffset(Seq::Id i=-1, int o0=0, int o1=1) : id(i), offset0(o0), offset1(o1) {}
    bool operator == (const ReadOffset& r) const { return id == r.id && offset0 == r.offset0 && offset1 == r.offset1; }
    static ReadOffset Make(const Overlap& o) { return {o.a_.id, o.a_.start, o.b_.start}; }
    Seq::Id id;
    uint16_t offset0;
    uint16_t offset1;
};

struct ReadInfo {
    Seq::Id id;
    const class Overlap *o;
    std::vector<std::array<int,5>> vars;    // ctg_i, rd_i, rd_b, choice, work

    int GetVar(size_t ctg) const {
        for (const auto &v : vars) {
            if (v[0] == (int)ctg) return v[2];
        }
        return -1;
    }
};

struct SnpSite {
    uint32_t ctg;       
    uint32_t offset;    
    bool operator == (const SnpSite& a) const {
        return a.ctg == ctg && a.offset == offset;
    }
    struct Hash {
        size_t operator() (const SnpSite &a) const {
            return std::hash<uint64_t>()(((uint64_t)a.ctg << 32) + a.offset);
        }
    };
};

struct SnpAllele {
    SnpSite site;
    uint8_t base;
    bool operator == (const SnpAllele& a) const {
        return a.site == site && a.base == base;
    }
    bool operator != (const SnpAllele& a) const {
        return !(*this == a);
    }
    struct Hash {
        size_t operator() (const SnpAllele &a) const {
            return 31 * SnpSite::Hash()(a.site) + a.base;
        }
    };
};
// 
struct PhaseItem {
    Seq::Id id;
    int strand : 8;
    int offset : 24;
};

} // namespace fsa


namespace std {
    template<>
    struct hash<fsa::ReadOffset> {
        size_t operator() (const fsa::ReadOffset& s) const noexcept {
            return hash<uint64_t>()(((uint64_t)s.id << 32) + (s.offset0 << 16) + s.offset1);
        }
    };
}