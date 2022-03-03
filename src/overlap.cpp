#include "overlap.hpp"

#include <sstream>
#include <vector>

#include "./utils/logger.hpp"
#include "utility.hpp"

namespace fsa {

void Overlap::Filter::From(const std::string &str) {
    auto items = SplitStringByChar(str, ':');

    for (auto &i : items) {
        auto kv = SplitStringByChar(i, '=');
        if (kv[0] == "l") {
            min_length = std::stoi(kv[1]);
        } else if (kv[0] == "i") {
            min_identity = std::stoi(kv[1]);
        } else if (kv[0] == "al") {
            min_aligned_length = std::stoi(kv[1]);

        } else if (kv[0] == "alr") {
            min_aligned_rate = std::stod(kv[1]);

        } else if (kv[0] == "oh") {
            max_overhang = std::stoi(kv[1]);
        } else if (kv[0] == "ohr") {
            max_overhang_rate = std::stod(kv[1]);
        } else if (kv[0] == "aal") {
            min_accept_aligned_length = std::stoi(kv[1]);
        } else if (kv[0] == "aalr") {
            min_accept_aligned_rate = std::stof(kv[1]);
        } else if (kv[0] == "ilid") {       // igore large indel
            large_indel_ = std::stoi(kv[1]);
        } else {
            LOG(ERROR)("Unrecoginze filtering option %s", kv[0].c_str());
        }
    }
}

std::string Overlap::Filter::ToString() const  {
    std::ostringstream oss;
    oss.precision(2);
    oss.setf(std::ios::fixed);
    oss << "l=" << min_length
        << ":i=" << min_identity
        << ":al=" << min_aligned_length
        << ":alr=" << min_aligned_rate
        << ":oh="  << max_overhang
        << ":ohr=" << max_overhang_rate
        << ":aal=" << min_accept_aligned_length
        << ":aalr=" << min_accept_aligned_rate
        << ":ilid=" << large_indel_;
    return oss.str();
}


double Overlap::Filter::IdentityIgnoringIndel(const Overlap& o, int indel) const {
    int big_dellen = 0;
    int big_inlen = 0;
    int inlen = 0;
    for (auto &d : o.detail_) {
        if (d.type == 'I') {
            if (d.len >= indel) big_dellen += d.len;
            inlen += d.len;
        } else if (d.type == 'D') {
            if (d.len >= indel) big_inlen += d.len;
        }
    }

    int len = o.b_.end - o.b_.start + inlen;
    int match = len * o.identity_ / 100;
    assert(len > big_dellen + big_inlen);
    return match * 100.0 / (len - big_dellen - big_inlen);
}

bool Overlap::Filter::Valid(const Overlap &ol, int replen, bool bothside) const {
    if (ol.a_.id == ol.b_.id) return false;

    if (ol.a_.len < min_length || ol.b_.len < min_length)   return false;
    if (Identity(ol) < min_identity)                        return false;

    int vlen = std::max<int>(ol.AlignedLength(), replen);
    if (vlen < min_aligned_length) {
        if (min_aligned_rate >= 0) {
            if (vlen < ol.a_.len * min_aligned_rate && 
                vlen < ol.b_.len * min_aligned_rate ) return false;
        }
    }

    if (min_accept_aligned_length >= 0 && vlen >= min_accept_aligned_length)   return true;
    if (min_accept_aligned_rate >= 0.0 && (
         vlen >= ol.a_.len * min_accept_aligned_rate || 
         vlen >= ol.b_.len * min_accept_aligned_rate)) return true;
    
    if (max_overhang >= 0 || max_overhang_rate >= 0) {
        int aoh = 0,  boh = 0;
        if (max_overhang >= 0 && max_overhang_rate >= 0) {
            aoh = std::min<int>(max_overhang, ol.a_.len*max_overhang_rate);
            boh = std::min<int>(max_overhang, ol.b_.len*max_overhang_rate);
        } else {
            aoh = max_overhang >= 0 ? max_overhang : (int)(ol.a_.len*max_overhang_rate);
            boh = max_overhang >= 0 ? max_overhang : (int)(ol.b_.len*max_overhang_rate);
        }

        aoh += replen;
        boh += replen;

        if (bothside) {
            if (ol.SameDirect()) {
                if (ol.a_.start > aoh && ol.b_.start > boh) return false;
                if (ol.a_.len - ol.a_.end > aoh && ol.b_.len - ol.b_.end > boh )  return false;
            } else {
                if (ol.a_.start > aoh && ol.b_.len - ol.b_.end > boh) return false;
                if (ol.a_.len - ol.a_.end > aoh && ol.b_.start > boh)  return false;
            }
        } else {
            if (ol.SameDirect()) {
                if (ol.a_.start > aoh && ol.b_.start > boh && ol.a_.len - ol.a_.end > aoh && ol.b_.len - ol.b_.end > boh)
                    return false;

            } else {
                if (ol.a_.start > aoh && ol.b_.len - ol.b_.end > boh && ol.a_.len - ol.a_.end > aoh && ol.b_.start > boh) 
                    return false;
            }

        }

    }
    return true;
}

bool Overlap::Filter::ValidQuery(const Overlap &ol) const {
    if (ol.a_.id == ol.b_.id) return false;

    if (ol.a_.len < min_length || ol.b_.len < min_length)   return false;
    if (Identity(ol) < min_identity)                        return false;

    if (ol.AlignedLength() < min_aligned_length) {
        if (min_aligned_rate >= 0) {
            if (ol.AlignedLength() < ol.a_.len * min_aligned_rate && 
                ol.AlignedLength() < ol.b_.len * min_aligned_rate ) return false;
        }
    }

    if (min_accept_aligned_length >= 0 && ol.AlignedLength() >= min_accept_aligned_length)   return true;
    if (min_accept_aligned_rate >= 0.0 && (
         ol.AlignedLength() >= ol.a_.len * min_accept_aligned_rate || 
         ol.AlignedLength() >= ol.b_.len * min_accept_aligned_rate)) return true;
    
    if (max_overhang >= 0 || max_overhang_rate >= 0) {
        int aoh = 0,  boh = 0;
        if (max_overhang >= 0 && max_overhang_rate >= 0) {
            aoh = std::min<int>(max_overhang, ol.a_.len*max_overhang_rate);
            boh = std::min<int>(max_overhang, ol.b_.len*max_overhang_rate);
        } else {
            aoh = max_overhang >= 0 ? max_overhang : (int)(ol.a_.len*max_overhang_rate);
            boh = max_overhang >= 0 ? max_overhang : (int)(ol.b_.len*max_overhang_rate);
        }
        if (ol.SameDirect()) {
            if (ol.a_.start > aoh && ol.b_.start > boh) return false;
            if (ol.a_.len - ol.a_.end > aoh && ol.b_.len - ol.b_.end > boh )  return false;
        } else {
            if (ol.a_.start > aoh && ol.b_.len - ol.b_.end > boh) return false;
            if (ol.a_.len - ol.a_.end > aoh && ol.b_.start > boh)  return false;
        }

    }
    return true;
}

std::string Overlap::ToM4Line() const {
    std::ostringstream oss;

    oss << a_.id + 1 << " " << b_.id + 1 << " " << identity_ << " " << -(int)AlignedLength() << " "
        << a_.strand << " " << a_.start << " " << a_.end << " " << a_.len << " "
        << b_.strand << " " << b_.start << " " << b_.end << " " << b_.len;

    return oss.str();
}

bool Overlap::FromM4Line(const std::string &line) {

    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 12) {

        // M4文件的Id就是read在fasta文件的序号。
        a_.id = atoi(items[0].c_str()) - 1;
        b_.id = atoi(items[1].c_str()) - 1;

        identity_ = atof(items[2].c_str());

        a_.strand = atoi(items[4].c_str());
        a_.start = atoi(items[5].c_str());
        a_.end = atoi(items[6].c_str());
        a_.len = atoi(items[7].c_str());

        b_.strand = atoi(items[8].c_str());
        b_.start = atoi(items[9].c_str());
        b_.end = atoi(items[10].c_str());
        b_.len = atoi(items[11].c_str());

        // 调整strand，保证a.strand = 0, 先设置b，再设置a
        b_.strand = a_.strand == b_.strand ? 0 : 1;
        a_.strand = 0;

        return true;
    }
    else {
        return false;
    }
}


bool Overlap::CheckEnd(int error) const {
    int a_start = a_.strand == 0 ? a_.start : a_.len - a_.end;
    int a_end = a_.strand == 0 ? a_.end : a_.len - a_.start;

    int b_start = b_.strand == 0 ? b_.start : b_.len - b_.end;
    int b_end = b_.strand == 0 ? b_.end : b_.len - b_.start;

    return !((a_start > error && b_start > error) ||
        (a_end < a_.len - error && b_end < b_.len - error));
}

bool Overlap::IsConsistent(const Overlap &ab, const Overlap &ac, const Overlap &bc, int err) {

    auto aid = ab.a_.id == ac.a_.id || ab.a_.id == ac.b_.id ? ab.a_.id : ab.b_.id;
    auto bid = ab.GetOtherRead(aid).id;
    auto cid = bc.GetOtherRead(bid).id;
    
    std::array<int, 2> src = {0, ab.GetRead(bid).len};
    std::array<int,2> dst0 = Overlap::Mapping<2>(ab.GetRead(aid), ab.GetRead(bid), src);

    std::array<int,2> dst1 = Overlap::Mapping<2>(ac.GetRead(aid), ac.GetRead(cid), Overlap::Mapping<2>(bc.GetRead(cid), bc.GetRead(bid), src));
    extern int a;
    if (a) printf("iii %d %d %d -- (%d %d) (%d %d)\n", ab.GetRead(aid).len, ab.GetRead(bid).len, ac.GetRead(cid).len, dst0[0], dst0[1], dst1[0], dst1[1]);
    return std::abs(dst0[0]-dst1[0]) <= err && std::abs(dst0[1]-dst1[1]) <= err;
}

Overlap::Loc Overlap::ReverseLocation(Loc loc, bool direct) {
    switch (loc) {
        case Loc::Left: return direct ? Loc::Right : Loc::Left;
        case Loc::Right: return direct ? Loc::Left : Loc::Right;
        case Loc::Equal: return Loc::Equal;
        case Loc::Containing: return Loc::Contained;
        case Loc::Contained: return Loc::Containing;
        case Loc::Abnormal: return Loc::Abnormal;
        default:    { assert(!"never come here"); return Loc::Abnormal; }
    }
}   

Overlap::Loc Overlap::Location(int err) const {
    auto compare = [err](int a, int b) { 
        if (a < b - err) return -1;
        else if (a > b + err) return 1;
        else return 0;
    };

    auto a_start_loc = compare(a_.start, 0);
    auto a_end_loc = compare(a_.end, a_.len);
    auto b_start_loc = compare(b_.start, 0);
    auto b_end_loc = compare(b_.end, b_.len);
    if (SameDirect()) {

        // a: --------->
        // b:       -------->
        if (a_start_loc > 0 && a_end_loc == 0 && b_start_loc == 0 && b_end_loc < 0) {
                return Loc::Left;
        }

        // a:       -------->
        // b:  ---------->
        if (a_start_loc == 0 && a_end_loc < 0 && b_start_loc > 0 && b_end_loc == 0) {
                return Loc::Right;
        }

        // a: -------->
        // b: -------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Equal;
        }

        // a: --------->
        // b:    ---->
        if (a_start_loc >= 0 && a_end_loc <= 0 &&(a_start_loc!=0 || a_end_loc !=0) && b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Containing;
        }

        // a:   ----->
        // b: ---------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc >= 0 && b_end_loc <= 0 && (b_start_loc != 0 || b_end_loc != 0)) {
            return Loc::Contained;
        }
    } else {

        // query :  <---------
        // target:       -------->
        if (a_start_loc == 0 && a_end_loc < 0 && b_start_loc == 0 && b_end_loc < 0) {
                return Loc::Left;
        }

        // query :       <---------
        // target:  ---------->
        if (a_start_loc > 0 && a_end_loc == 0 && b_start_loc > 0 && b_end_loc == 0) {
                return Loc::Right;
        }

        // query : <--------
        // target: -------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Equal;
        }

        // query : <---------
        // target:    ---->
        if (a_start_loc >= 0 && a_end_loc <= 0 &&(a_start_loc != 0 || a_end_loc != 0)&& b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Containing;
        }

        // query :   <----
        // target: ---------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc >= 0 && b_end_loc <= 0 && (b_start_loc != 0 || b_end_loc != 0)) {
            return Loc::Contained;
        }
    }

    return Loc::Abnormal;
}


Overlap::Loc Overlap::Location(Seq::Id id, int err) const {
    assert(id == b_.id || id == a_.id);
    Loc loc = Location(err);
    return id == a_.id ? loc : ReverseLocation(loc, SameDirect());
}

bool Overlap::IsContaining(int err) const {
    Loc loc = Location(err);
    return loc == Loc::Equal || loc == Loc::Containing;
}

bool Overlap::IsContained(int err) const {
    Loc loc = Location(err);
    return loc == Loc::Equal || loc == Loc::Contained;
}


bool Overlap::IsProper(int err) const {
    Loc loc = Location(err);
    return loc == Loc::Left || loc == Loc::Right;
}


std::array<int, 2> Overlap::Overhang(const Overlap &o, Loc loc) {
    std::array<int,2> overhang({-1, -1}); // -1 mean no overhang
    
    if (o.SameDirect()) {
        if (loc == Overlap::Loc::Left) {
            //int oh = std::max(o.b_.start-0, o.a_.len-o.a_.end);
            overhang[0] = o.a_.len - o.a_.end;
            overhang[1] = o.b_.start;
        } else if (loc == Overlap::Loc::Right) {
            //int oh = std::max(o.a_.start-0, o.b_.len-o.b_.end);
            overhang[0] = o.a_.start;
            overhang[1] = o.b_.len - o.b_.end;
        } else if (loc == Overlap::Loc::Contained) {
            //int oh = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
        } else if (loc == Overlap::Loc::Containing) {
            //int oh = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        } else if (loc == Overlap::Loc::Equal) {
            //int oh1 = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            //int oh2 = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            //int oh = std::max(oh1, oh2);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        }
    } else {
        if (loc == Overlap::Loc::Left) {
            //int oh = std::max(o.b_.start-0, o.a_.start-0);
            overhang[0] = o.a_.start;
            overhang[1] = o.b_.start;
        } else if (loc == Overlap::Loc::Right) {
            //int oh = std::max(o.a_.len-o.a_.end, o.b_.len-o.b_.end);
            overhang[0] = o.a_.len - o.a_.end;
            overhang[1] = o.b_.len - o.b_.end;
        } else if (loc == Overlap::Loc::Contained) {
            //int oh = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
        } else if (loc == Overlap::Loc::Containing) {
            //int oh = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        } else if (loc == Overlap::Loc::Equal) {
            //int oh1 = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            //int oh2 = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            //int oh = std::max(oh1, oh2);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        }
    }
    return overhang;   
}


std::array<int, 2> Overlap::Overhang() const {
    std::array<int,2> result;
    int a0 = a_.start - 0;
    int a1 = a_.len - a_.end;
    int b0 = b_.start - 0;
    int b1 = b_.len - b_.end;

    if (SameDirect()) {
        std::vector<int> oh{a0+b1, a1+b0, a0+a1, b0+b1};

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{a0, b1};
        } else if (min_oh == 1) {
            result = std::array<int,2>{a1, b0};
        } else if (min_oh == 2) {
            result = std::array<int,2> {std::max(a0, a1), 0};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {0, std::max(b0, b1)};
        }
    } else {
        std::vector<int> oh {a0+b0, a1+b1, a0+a1, b0+b1};

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{a0, b0};
        } else if (min_oh == 1) {
            result = std::array<int,2>{a1, b1};
        } else if (min_oh == 2) {
            result = std::array<int,2> {std::max(a0, a1), 0};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {0, std::max(b0, b1)};
        }
    }
    return result;
}

std::array<int, 2> Overlap::Overhang2() const {
    std::array<int,2> result;   // 0表示没有判断为overhang，1表示左端为overhnag 2表示右端为overhang 3表示两端为overhang
    int a0 = a_.start - 0;
    int a1 = a_.len - a_.end;
    int b0 = b_.start - 0;
    int b1 = b_.len - b_.end;

    if (SameDirect()) {
        std::vector<int> oh(4);
        oh[0] = a0 + b1;
        oh[1] = a1 + b0;
        oh[2] = a0 + a1;
        oh[3] = b0 + b1;

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{1, 2};
        } else if (min_oh == 1) {
            result = std::array<int,2>{2, 1};
        } else if (min_oh == 2) {
            result = std::array<int,2> {3, 0};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {0, 3};
        }
    } else {
        std::vector<int> oh(4);
        oh[0] = a0 + b0;
        oh[1] = a1 + b1;
        oh[2] = a0 + a1;
        oh[3] = b0 + b1;

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{1, 1};
        } else if (min_oh == 1) {
            result = std::array<int,2>{2, 2};
        } else if (min_oh == 2) {
            result = std::array<int,2> {3, 0};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {0, 3};
        }
    }
    return result;
}



bool Overlap::Extend(int maxoh) {
    if (Location(0) != Overlap::Loc::Abnormal) return true;

    if (Location(maxoh) != Overlap::Loc::Abnormal) {

        if (a_.strand == b_.strand) {
            if (a_.start <= maxoh && b_.start <= maxoh) {
                a_.start = 0;
                b_.start = 0;
            }
            else if (a_.start <= maxoh) {
                b_.start -= a_.start;
                a_.start = 0;
            } 
            else if (b_.start <= maxoh) {
                a_.start -= b_.start;
                b_.start = 0;
            }

            if (a_.end >= a_.len - maxoh && b_.end >= b_.len - maxoh) {
                a_.end = a_.len;
                b_.end = b_.len;
            }
            else if (a_.end >= a_.len - maxoh) {
                b_.end += a_.len - a_.end;
                a_.end = a_.len;
            }
            else if (b_.end >= b_.len - maxoh) {
                a_.end += b_.len - b_.end;
                b_.end = b_.len;
            }

        }
        else {
            if (a_.start <= maxoh && b_.end >= b_.len - maxoh) {
                a_.start = 0;
                b_.end = b_.len;
            }
            else if (a_.start <= maxoh) {
                b_.end += a_.start;
                a_.start = 0;
            }
            else if (b_.end >= b_.len - maxoh) {
                a_.start -= b_.len - b_.end;
                b_.end = b_.len;
            }

            if (b_.start <= maxoh && a_.end >= a_.len - maxoh) {
                b_.start = 0;
                a_.end = a_.len;
            }
            else if (b_.start <= maxoh) {
                a_.end += b_.start;
                b_.start = 0;

            }
            else if (a_.end >= a_.len - maxoh) {
                b_.start -= a_.len - a_.end;
                a_.end = a_.len;
            }
        }

        assert(Location(0) != Overlap::Loc::Abnormal);

        return true;
    } else {
        return false;
    }

 
}


int Overlap::Extension(Seq::Id id, int end) const {
    assert(end == 0 || end == 1);
    auto tr = GetRead(id);
    auto qr = GetOtherRead(id);

    if (SameDirect()) {
        if (end == 0) {
            return std::max(0, qr.start-tr.start);
        } else {
            return std::max(0, (qr.len-qr.end) - (tr.len-tr.end));
        }
    } else {
        if (end == 0) {
            return std::max(0, (qr.len-qr.end)-tr.start);
        } else {
            return std::max(0, qr.start - (tr.len-tr.end));
        }
    }
}

} // namespace fsa {