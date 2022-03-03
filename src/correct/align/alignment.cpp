#include "alignment.hpp"

namespace fsa {

void Alignment::Swap(bool isSameDirect) {
    std::swap(target, query);
    std::swap(target_start, query_start);
    std::swap(target_end, query_end);
    std::swap(aligned_target, aligned_query);

    if (!isSameDirect) {
        std::swap(target_start, target_end);
        target_start = target->Size() - target_start;
        target_end = target->Size() - target_end;

        std::swap(query_start, query_end);
        query_start = query->Size() - query_start;
        query_end = query->Size() - query_end;

        Seq::ReverseComplementInPlace(aligned_target);
        Seq::ReverseComplementInPlace(aligned_query);
 
    }
}

void Alignment::Rearrange(std::string &alq, std::string &alt) {
    
    // 原则: 将alt的base尽量往
    //       在保证alt的前提下，alq尽量前移
    //       alt的长度尽可能短

    // 将mdf尽量向前移动，ref不改变
    //  ...CGX...     ---\    ...CGX... 
    //  ...--C...     .../    ...C--...
    //
    auto move_forward = [](std::string &mdf, const std::string& ref) {
        bool modified = false;
        for (size_t i=1; i<mdf.size(); i++) {
            if (mdf[i-1] == '-' && mdf[i] != '-') {
                int cand = -1;
                for (int j = i-1; j>=0 && mdf[j] == '-'; j--) {
                    if (ref[j] == mdf[i]) {
                        cand = j;
                    }
                }

                if (cand != -1 ) {
                    std::swap(mdf[i], mdf[cand]);
                    modified = true;
                }
            }
        }
        return modified;
    };


    //  qlt: ...-GA...     ---\    ...-GA...   ---\    ...GA... 
    //  alt: ...A--...     .../    ...--A...   ---/    ...-A...  
    //
    auto move_backward = [](const std::string &alq, std::string& alt) {
        bool modified = false;
        for (size_t i=0; i+1<alt.size(); i++) {
            if (alt[i] != '-' && alt[i+1] == '-' && alq[i] == '-') {
                int cand = -1;
                for (size_t j = i+1; j < alt.size() && alt[j] == '-'; j++) {
                    if (alq[j] == alt[i]) {
                        cand = j;
                    }
                }

                if (cand != -1 ) {
                    std::swap(alt[i], alt[cand]);
                    modified = true;
                }
            }
        }
        return modified;
    };

    //
    // alq: ...AAA-...       ---\    ...-AAA...
    // alt: ...---C...       ---/    ...C---...
    //
    auto swap_dash_forward = [](std::string &alq, std::string& alt) {
        bool modified = false;
        for (size_t i=1; i<alt.size(); i++) {
            if (alt[i-1] == '-' && alt[i] != '-' && alq[i] == '-') {
                for (int j = i; j>=1 && alt[j-1] == '-'; j--) {
                    
                    std::swap(alt[j], alt[j-1]);
                    std::swap(alq[j], alq[j-1]);
                    modified = true;
                }
            }

        }
        return modified;
    };

    // 删除dash
    //  alq: ...A-C...       ---\    ...AC...
    //  alt: ...A-C...       ---/    ...AC...
    //
    //  alq: ...-A..A-...       ---\    ...A..A...
    //  alt: ...A-..-A...       ---/    ...A..A...
    auto remove_dash = [] (std::string &alq, std::string &alt) {
        size_t oldsize = alq.size();
        for (size_t i=0; i<alt.size(); ++i) {
            if (alq[i] == alt[i] && alt[i] == '-') {
                alt.erase(alt.begin()+i);
                alq.erase(alq.begin()+i);
                i--;
            }
        }

        for (size_t i=1; i< alt.size(); ++i) {
            if (alt[i-1] == '-' && alq[i] == '-'  && alt[i] == alq[i-1]) {
                alt.erase(alt.begin()+i-1);
                alq.erase(alq.begin()+i);
                i--;
                //std::swap(alt[i-1], alt[i]);
            } else if (alt[i] == '-' && alq[i-1] == '-'  && alt[i-1] == alq[i]) {
                alt.erase(alt.begin()+i);
                alq.erase(alq.begin()+i-1);
                i--;
                //std::swap(alq[i-1], alq[i]);
            }
        }
        return oldsize != alq.size();   // modified
    };

    bool finished = false;
    while (!finished) {
        finished = true;
        finished &= !move_forward(alt, alq);
        finished &= !move_forward(alq, alt);
        finished &= !move_backward(alq, alt);
        finished &= !move_backward(alt, alq);

        finished &= !swap_dash_forward(alq, alt);

        break;

    }

}


void Alignment::Rearrange1(std::string &alq, std::string &alt) {

    for (size_t i=0; i<alq.size() - 1; ++i) {
        if (alt[i] == '-') {

            for (size_t j=i+1; j<alq.size(); ++j) {
                if (alt[j] != '-') {
                    if (alq[i] == alt[j]) {
                        alt[i] = alt[j];
                        alt[j] = '-';
                    }
                    break;
                }
            }
        }

        if (alq[i] == '-') {
            for (size_t j=i+1; j<alq.size(); ++j) {
                if (alq[j] != '-') {
                    if (alt[i] == alq[j]) {
                        alq[i] = alq[j];
                        alq[j] = '-';
                    }
                    break;
                }
            }
        }

    }
}

bool Alignment::TrimEnds(size_t checklen, int stub) {
    if (AlignSize() < checklen) return false;

    int sc = 0;
    size_t as = 0, ts = target_start, qs = query_start;
    size_t ias = 0, its = target_start, iqs = query_start;

    for (ias = 0; ias < checklen; ++ias) {
        if (aligned_query[ias] == aligned_target[ias]) {
            if (aligned_query[ias] != '-') {
                if (sc == 0) {
                    as = ias;
                    qs = iqs;
                    ts = its;
                }
                sc ++;
                if (sc >= stub) break;
            } 
        } else  {
            sc = 0;
        } 

        // increase iqs, its
        if (aligned_query[ias] != '-') iqs++;
        if (aligned_target[ias] != '-') its++;
    }

    int ec = 0;
    size_t ae = 0, te = target_end, qe = query_end;
    size_t iae = aligned_target.size(), ite = target_end, iqe = query_end;

    for (; iae > aligned_target.size() - checklen; iae -- ) {
        if (aligned_query[iae-1] == aligned_target[iae-1]) {
            if (aligned_query[iae-1] != '-') {
                if (ec == 0) {
                    ae = iae;
                    qe = iqe;
                    te = ite;
                }
                ec ++;
                if (ec >= stub) break;
            }
        } else {
            ec = 0;
        }
        
        // decrease iqe, ite
        if (aligned_query[iae-1] != '-') iqe--;
        if (aligned_target[iae-1] != '-') ite--; 
    }


    if (sc >= stub && ec >= stub) {
        target_start = ts;
        query_start = qs;

        target_end = te;
        query_end = qe;

        aligned_target = aligned_target.substr(as, ae-as);
        aligned_query  = aligned_query.substr(as, ae-as);
        return true;
    } else {
        return false;
    }
}

int Alignment::Score() const {
    int score = 0;
    const int ins = 1;
    const int del = 1;
    const int mism = 2; 
    for (size_t i=0; i<AlignSize(); ++i) {
        auto bs = GetAlign(i);
        if (bs[0] != bs[1]) {
            if (bs[0] == '-') {
 
            }
        }
    }
    return score;
}

} // namespace fsa {