#include "sequence_store.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>

#include "file_io.hpp"
#include "./utils/logger.hpp"
#include "./utility.hpp"

namespace fsa {

void BaseStore::Add(const char* str, size_t len) {
    // check block
    assert(len < (size_t)block_size_);
    if (size_ + len >= data_.size() * block_size_ * 4) {
        data_.push_back(std::vector<uint8_t>());
        data_.back().assign(block_size_, 0);
    }

    for (size_t i = 0; i < len; ++i) {
        auto &d = Locate(size_ + i);
        auto off = (size_ + i) % 4;
        d += DnaSeq::s_SerialTable[str[i]] << (off)*2;
    }
    size_ += len;
}


void SequenceStore::Load(const std::string &fname, const std::string &type) {
    LoadFastq(fname);
}


void SequenceStore::LoadFastq(const std::string &fname) {
    std::mutex mutex_gen;
    std::mutex mutex_combine;

    GzFileReader in(fname);

    size_t load_threads = 8;
    size_t thread_size = 8;
    if (thread_size > load_threads) thread_size = load_threads;
    std::atomic<size_t> block_index { 0 };
    
    auto generate_func = [&mutex_gen, &in, &block_index](std::vector<char> &block) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        size_t curr = block_index.fetch_add(1);
        LOG(INFO)("BLOCK_INDEX: %zd, mem=%zd", curr, GetMemoryUsage());

        return in.GetBlock1(block, FindLastFastqInBlock);
    };

    auto combine_func = [&mutex_combine, this](const std::vector<char> &block, size_t block_size, std::vector<std::vector<std::array<size_t, 2>>>& segments) {
        std::lock_guard<std::mutex> lock(mutex_combine);

        // for (auto &segs : segments) {
        //     // std::string name(block.begin()+segs[0][0], block.begin()+segs[0][1]);
        //     // auto id = string_pool_.GetIdByStringUnsafe(name);
        //     // items_.push_back(Sequence());
        //     // items_.back().offset = data_.Size();
        //     // items_.back().len = segs[2][1] - segs[2][0];
        //     // data_.Add(&block[segs[0][0]], items_.back().len);

        // }
  
        segments.clear();
    };



    auto work_func = [generate_func, combine_func, &fname, &in, this](size_t threadid) {
        std::vector<char> block(100000000);
        size_t block_size = generate_func(block);

        std::vector<std::vector<std::array<size_t, 2>>> segments;
        
        while (block_size > 0) {
            for (auto segs = GetFastqSegment(block, block_size); segs.size() == 4; segs = GetFastqSegment(block, block_size, segs.back()[1])) {
                segments.push_back(segs);

            }
            LOG(INFO)("thread=%zd, size=%zd, mem=%zd", threadid, items_.size(), GetMemoryUsage());

            combine_func(block, block_size, segments);
            block_size = generate_func(block);
        }

    
        if (!in.IsEnd()) {
            LOG(WARNING)("No all read in file are loaded: %s", fname.c_str());
        }
    };

    if (in.Valid()) {
        printf("thread_size = %zd", thread_size);
        MultiThreadRun(thread_size, work_func);
    } else {
        LOG(ERROR)("Failed to load file: %s", fname.c_str());
    }
}

std::vector<std::array<size_t, 2>> SequenceStore::GetFastqSegment(const std::vector<char> &block, size_t bsize, size_t start) {

    assert(bsize <= block.size() && start <= bsize);

    std::vector<std::array<size_t,2>> segs;

    auto consume_tochar = [&block, bsize]( size_t &p, char c) {
        while (p < bsize && block[p] != c) p++;
    };

    auto consume_nospace = [&block, bsize]( size_t &p) {
        while (p < bsize && !::isspace(block[p])) p++;
    };

    auto check_char = [&block, bsize](size_t &p, char c) {
        return p < bsize && block[p] == c;
    };

    do {
        // get head 
        segs.push_back({ start, start});
        auto &head = segs.back();
        consume_tochar(head[0], '@');
        if (!check_char(head[0], '@')) break;
        head[0] ++;

        head[1] = head[0];
        consume_nospace(head[1]);

        // get subhead
        segs.push_back({ head[1], head[1]});
        auto &subhead = segs.back();
        consume_tochar(subhead[1], '\n');
        //while (subhead[0] < subhead[1] && ::isspace(block[subhead[0]])) subhead[0]++;
        //while (subhead[1] > subhead[0] && ::isspace(block[subhead[1]])) subhead[1]--;

        // get seq
        segs.push_back({ subhead[1], subhead[1]});
        auto &seq = segs.back();
        if (!check_char(seq[0], '\n')) break;
        seq[0]++;
        seq[1] = seq[0];
        consume_tochar(seq[1], '\n');

        // get quality
        segs.push_back({ seq[1], seq[1]});
        auto &quality = segs.back();
        if (!check_char(quality[0], '\n')) break;
        quality[0] ++; // skip '\n'
        if (!check_char(quality[0], '+')) break;

        consume_tochar(quality[0], '\n');
        if (!check_char(quality[0], '\n')) break;
        quality[0] ++; // skip '\n'
        quality[1] = quality[0];
        consume_nospace(quality[1]);

        return segs;

    } while (false);

    segs.clear();
    return segs;

}


size_t FindLastFastqInBlock(const std::vector<char> &block, size_t bsize) {
    auto check_plus = [](const std::vector<char> &block, size_t bsize) {
        size_t line_no = 0;
        for (size_t i = bsize; i > 0; --i) {
            if (block[i-1] == '\n') {
                line_no ++;
                if (line_no == 2) {
                    // assert(i < bsize);
                    return block[i] == '+';
                }
            }
        }
        return false;
    };

    for (size_t i = bsize; i >= 2; --i) {
        if (block[i-1] == '@' && block[i-2] == '\n') {
            if (check_plus(block, i-2)) {   // skip '\n'
                return i-1;
            }
        }
    }

    return 0;
}

} // namespace fsa {
