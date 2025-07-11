#include "match_info.hpp"

#include <numeric>

#include "utils/logger.hpp"

namespace fsa {


MatchInfo::MatchInfo(const Overlap* ol, const DnaSeq& qseq, const DnaSeq& tseq) : ol_(ol) {
    match_.assign(ol->b_.end - ol->b_.start, InfoItem());

    auto get_base = [](const Overlap::Read &r, const DnaSeq& seq, size_t idx) {
        return r.strand == 0 ? seq[r.start+idx] : (3 - seq[r.end - idx - 1]);
    };
    
    size_t qidx = 0;        // not from ol.b_.start;
    size_t tidx = 0;

    size_t matched = 0;
    size_t matched_len = 0;
    size_t SMALL_INDEL = 50;

    for (const auto &d : ol->detail_) {
        switch (d.type){
        case 'M':
        case '=':
        case 'X':
            for (size_t i = 0; i < (size_t)d.len; ++i) {
                uint8_t cq = get_base(ol->a_, qseq, qidx+i);
                uint8_t ct = tseq[tidx+i + ol->b_.start];
                match_[tidx+i].ref = ct;
                match_[tidx+i].base = cq;
                if (ct == cq) {
                    matched ++;
                }
            }
            qidx += d.len;
            tidx += d.len;
            matched_len += d.len;
            break;
        case 'D':
            for (size_t i = 0; i < (size_t)d.len; ++i) {
                char ct = tseq[tidx+i + ol->b_.start];
                match_[tidx+i].ref = ct;
                match_[tidx+i].base = 4;
            }
            tidx += d.len;
            if (d.len <= SMALL_INDEL) {
                matched_len += d.len;
            }
            break; 
        case 'I':
            insert_.push_back({qidx, qidx+d.len});
            match_[tidx].ins = insert_.size();

            qidx += d.len;
            if (d.len <= SMALL_INDEL) {
                matched_len += d.len;
            }
            break;
        case 'S':
        case 'H':
            break; // do nothing
        default:
            LOG(ERROR)("never come here %c", d.type);
        }
    } 
    matched_identity_ = matched * 1.0 / matched_len;
    CalculateMaxLocalDistance(1000);    // TODO: make it configurable
    for (size_t i = 0; i < match_.size(); ++i) {
        assert(match_[i].ref  == tseq[i + ol->b_.start]);
    }
}

void MatchInfo::CalculateMaxLocalDistance(size_t win_size) {
    if (match_.size() < win_size) {
        LOG(INFO)("%d %d %d", ol_->a_.start, ol_->a_.end, ol_->a_.len);
    }
    win_size = std::min(win_size, match_.size());
    assert(match_.size() >= win_size);

    std::vector<std::array<uint32_t,2>> dist(match_.size(), {0, 0});
    for (size_t i=0; i < match_.size(); ++i) {
        if (match_[i].ref == match_[i].base) {
            dist[i][0] = 1;
        }
        dist[i][1] = 1 + GetInssize(match_[i].ins);
    }

    std::array<uint32_t,2> curr_dist = std::accumulate(dist.begin(), dist.begin() + win_size, std::array<uint32_t,2>({0,0}),
        [](const std::array<uint32_t,2>&a, const std::array<uint32_t,2> &b) -> std::array<uint32_t,2> {
            return {a[0] + b[0], a[1] + b[1]};
        }
    );

    double max_dist = 1 - curr_dist[0] * 1.0 / curr_dist[1];

    for (size_t i = win_size; i < dist.size(); ++i) {
        curr_dist[0] += dist[i][0];
        curr_dist[0] -= dist[i-win_size][0];
        curr_dist[1] += dist[i][1];
        curr_dist[1] -= dist[i-win_size][1];
        double d = 1 - curr_dist[0] * 1.0 / curr_dist[1];
        if (d > max_dist) {
            max_dist = d;
        }
    }
    max_local_distance_ = max_dist;
}

std::vector<std::array<size_t,2>> MatchInfo::LocalDistance(size_t win_size)  const {
    win_size = std::min(win_size, match_.size());
    assert(match_.size() >= win_size);

    std::vector<std::array<size_t,2>> dist(match_.size(), {0, 0});
    std::vector<std::array<size_t,2>> win_dist(match_.size() - win_size + 1, {0, 0});
    for (size_t i=0; i < match_.size(); ++i) {
        if (match_[i].ref == match_[i].base) {
            dist[i][0] = 1;
        }
        dist[i][1] = 1 + GetInssize(match_[i].ins);
    }

    win_dist[0] = std::accumulate(dist.begin(), dist.begin() + win_size, std::array<size_t,2>({0,0}),
        [](const std::array<size_t,2>&a, const std::array<size_t,2> &b) -> std::array<size_t,2> {
            return {a[0] + b[0], a[1] + b[1]};
        }
    );

    for (size_t i = win_size; i < dist.size(); ++i) {
        win_dist[i-win_size+1][0] = win_dist[i-win_size][0] + dist[i][0] - dist[i-win_size][0];
        win_dist[i-win_size+1][1] = win_dist[i-win_size][1] + dist[i][1] - dist[i-win_size][1];
    }

    return win_dist;
}

std::vector<std::array<size_t,2>> MatchInfo::GetHighQualityRegions(size_t win_size, double max_dist, size_t min_clip, size_t min_intv) const {
    std::vector<std::array<size_t,2>> regs;

    auto local_dist = LocalDistance(win_size);
    std::vector<uint8_t> flag(match_.size(), 0);

    std::vector<std::array<size_t,2>> vregs;
    if (LClip() >= min_clip) {
        size_t e = min_intv > match_.size() ? match_.size() :  min_intv;
        vregs.push_back({0, e});
    }
    for (size_t i=0; i < local_dist.size(); ++i) {
        double d = 1 - local_dist[i][0]*1.0 / local_dist[i][1];
        if (d > max_dist) {
            size_t s = i < min_intv ? 0 : i - min_intv;
            size_t e = i + win_size + min_intv > match_.size() ? match_.size() : i + win_size + min_intv;
            vregs.push_back({s, e});
            //LOG(INFO)("vregs %zd-%zd %.02f > %.02f", s,e, d, max_dist);
        }
    }

    if (RClip() >= min_clip) {
        size_t s = min_intv > match_.size() ? 0 :  match_.size() - min_intv;
        vregs.push_back({s, match_.size()});
    }

    size_t s = 0;
    for (auto& vr : vregs) {
        if (s < vr[0]) {
            regs.push_back({s, vr[0]});
        }
        assert (s <= vr[1]);
        s = vr[1];
    }
    if (s < match_.size()) {
        regs.push_back({s, match_.size()});
    }

    return regs;
}

}