#ifndef FSA_OVERLAP_STORE_HPP
#define FSA_OVERLAP_STORE_HPP


#include <array>
#include <vector>
#include <list>
#include <deque>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "overlap.hpp"
#include "./utils/logger.hpp"
#include "utils/string_pool.hpp"
#include "sequence.hpp"
#include "utility.hpp"
#include "file_io.hpp"


namespace fsa {
struct OverlapSet {
    OverlapSet() { 
        AddBlock();
    }
    size_t Size() const { return size; }

    void Insert(const std::vector<Overlap> &ols, size_t sz) {
        for (size_t i=0; i<sz; ) {
            if (overlaps.back().size() == bsize) {
                AddBlock();
            }

            size_t s = std::min<size_t>(sz - i, bsize - overlaps.back().size());
            overlaps.back().insert(overlaps.back().end(), ols.begin()+i, ols.begin()+i+s);
            size += s;
            i += s;
        }
    }

    void Add(const Overlap& ol) {
        if (overlaps.back().size() == bsize) {
            AddBlock();
        }
        overlaps.back().push_back(ol);
        size += 1;
    }

    Overlap& Get(size_t i)  { size_t i0 = i / bsize; size_t i1 = i% bsize; return overlaps[i0][i1]; }
    const Overlap& Get(size_t i) const { size_t i0 = i / bsize; size_t i1 = i% bsize; return overlaps[i0][i1]; }

    void AddBlock() {
        overlaps.push_back(std::vector<Overlap>());
        overlaps.back().reserve(bsize);
    }
    size_t size = 0;
    size_t bsize = 10000000;
    std::vector<std::vector<Overlap>> overlaps;
};

class OverlapStore {
public:
    OverlapStore() : OverlapStore(default_string_pool_) {}
    OverlapStore(StringPool &sp) : string_pool_(sp) {}

    template<typename C=bool (*)(Overlap &o)>
    void Load(const std::string &fname, const std::string &type="", size_t thread_size=1, C check=[](Overlap &o) {return true; });

    template<typename C=bool (*)(Overlap &o)>
    void LoadFast(const std::string &fname, const std::string &type="", size_t thread_size=1, C check=[](Overlap &o) {return true; });

    static std::string DetectFileType(const std::string &fname);
    
    template<typename C = bool(*)(const Overlap &o)>
    void Save(const std::string &fname, const std::string &type="", size_t thread_size=1, C check=[](const Overlap &o){return true;}) const; 

    size_t Size() const { return overlaps_.Size(); }
    Overlap& Get(size_t i)  { return  overlaps_.Get(i); }
    const Overlap& Get(size_t i) const { return  overlaps_.Get(i); }
    // size_t Size() const { return overlaps_.size(); }
    // Overlap& Get(size_t i)  { return  overlaps_[i]; }
    // const Overlap& Get(size_t i) const { return  overlaps_[i]; }


    std::string GetReadName(int id) {
        return string_pool_.QueryStringById(id);
    }

    Seq::Id GetReadId(const std::string &name) {
        return string_pool_.GetIdByString(name);
    }

    std::array<Seq::Id, 2> GetReadIdRange() const;
    const StringPool& GetStringPool() const { return string_pool_; }

    std::unordered_map<int, std::unordered_map<int, Overlap*>> Group() const;
    std::unordered_map<int, std::unordered_map<int, Overlap*>> Group();
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> Group(bool (*better)(const Overlap&, const Overlap&)) const;
    std::unordered_map<int, std::unordered_map<int, Overlap*>> GroupTarget(bool (*better)(const Overlap* a, const Overlap *b));
    std::unordered_map<int, std::unordered_map<int, Overlap*>> GroupQuery(bool (*better)(const Overlap* a, const Overlap *b));
    void Group(std::unordered_map<int, std::unordered_map<int, const Overlap*>>& groups, size_t thread_size) const;
    void Group(std::unordered_map<int, std::unordered_map<int, const Overlap*>>& groups, const std::unordered_set<int>& keys, size_t thread_size) const;
    void GroupTarget(std::unordered_map<Seq::Id, std::unordered_map<Seq::Id, std::vector<const Overlap*>>> &groups, size_t threads) const;


    template<typename F, typename C>
    void LoadFileMt(const std::string &fname, F lineToOl, C check, size_t thread_size=1);

    template<typename F, typename C>
    void LoadFileFast(const std::string &fname, F lineToOl, C check, size_t thread_size);
 
    template<typename L, typename C>
    void SaveFile(const std::string &fname, L toLine, size_t thread_size, C check) const;

    template<typename C>
    void LoadFileTxtMt(const std::string &fname, C check, size_t thread_size);
    template<typename C>
    void LoadFileTxtFast(const std::string &fname, C check, size_t thread_size);

    static bool FromM4Line(const std::string &line, Overlap &o, StringPool::NameId& ni);
    static bool FromM4aLine(const std::string &line, Overlap &o, StringPool::NameId& ni);
    static bool FromPafLine(const std::string &line, Overlap &o, StringPool::NameId& ni);

    static bool FromM4LineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen) {
        replen = 0;
        return FromM4Line(line, o, ni);
    }
    static bool FromM4aLineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen) {
        replen = 0;
        return FromM4aLine(line, o, ni);
    }

    static bool FromPafLineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen);

    static std::string ToM4aLine(const Overlap& o, const  StringPool::NameId& ni);
    static std::string ToM4Line(const Overlap& o, const StringPool::NameId& ni);
    static std::string ToPafLine(const Overlap &o, const StringPool::NameId& ni) ;

    std::string ToM4aLine1(const Overlap& o) const { return ToM4aLine(o,  StringPool::UnsafeNameId(string_pool_)); }
    std::string ToM4Line1(const Overlap& o) const { return ToM4Line(o,  StringPool::UnsafeNameId(string_pool_)); }
    std::string ToPafLine1(const Overlap &o) const { return ToPafLine(o,  StringPool::UnsafeNameId(string_pool_)); }

    std::vector<std::string> GetLineFromFile(const std::string& fname) {
        GzFileReader in(fname);
        std::vector<std::string> lines;
        std::string s = in.GetNoEmptyLine();
        while (!s.empty()) {
            lines.push_back(s);
            s = in.GetNoEmptyLine();
        }
        return lines;
    }
protected:
    //std::deque<Overlap> overlaps_;
    OverlapSet overlaps_;
    std::mutex mutex_overlaps_;

    StringPool &string_pool_;
    StringPool default_string_pool_;

    size_t load_threads = 80;

};

template<typename F, typename C>
void OverlapStore::LoadFileMt(const std::string &fname, F lineToOl, C check, size_t thread_size) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;
    GzFileReader in(fname);
    const size_t block_size = 10000;
    if (thread_size > load_threads) thread_size = load_threads;

    auto generate_func = [&mutex_gen, &in](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        return in.GetLines(lines);
    };

    auto combine_func = [&check, &mutex_comb, this](std::vector<Overlap> &ols, size_t sz, StringPool::TempNameId &ni) {
        if (ni.names_to_ids.size() > 0) {
            auto id2id = string_pool_.MergeNameId(ni);
            ni.names_to_ids.clear();
            for (size_t i=0; i<sz; ++i) {
                auto& o = ols[i];
                o.a_.id = id2id[o.a_.id];
                o.b_.id = id2id[o.b_.id];
            }
        }
        {
            std::lock_guard<std::mutex> lock(mutex_comb);
            for (size_t i=0; i<sz; ++i) {
                auto& o = ols[i];
                if (check(o)) {
                    overlaps_.Add(o);
                }
            }
        }
    };

    auto work_func = [&check, this, lineToOl, block_size, generate_func, combine_func, &fname, &in](size_t id) {
        std::vector<std::string> lines(block_size);
        std::vector<Overlap> ols(block_size);
        thread_local StringPool::TempNameId ni;
        while (true) {
            size_t line_size = generate_func(lines);
            if (line_size > 0) {
                size_t ol_size = 0;
                for (size_t i=0; i<line_size; ++i) {
                    Overlap o;
                    if (lineToOl(lines[i], o, ni)) {
                        ols[ol_size++] = o;
                    } else {
                        LOG(ERROR)("Failed to convert line to overlap \n   %s", lines[i].c_str());
                    }
                }
                combine_func(ols, ol_size, ni);;
                
            } else {
                break;
            }

        }   
        if (!in.IsEnd()) {
            LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
        }
    };

    if (in.Valid()) {
        MultiThreadRun(thread_size, work_func);
    } else {
        LOG(ERROR)("Failed to load file: %s", fname.c_str());
    }
}


// template<typename F, typename C>
// void OverlapStore::LoadFileFast(const std::string &fname, F lineToOl, C check, size_t thread_size) {
//     std::mutex mutex_gen;
//     GzFileReader in(fname);

//     const size_t block_size = 10000;
//     if (thread_size > load_threads) thread_size = load_threads;
    
//     auto generate_func = [&mutex_gen, &in](std::vector<std::string> &lines) {
//         std::lock_guard<std::mutex> lock(mutex_gen);
//         return in.GetLines(lines);
//     };

//     auto combine_func = [&check, this](std::vector<Overlap> &ols, size_t sz) {
//         std::lock_guard<std::mutex> lock(mutex_overlaps_);
//         //overlaps_.insert(overlaps_.end(), ols.begin(), ols.begin()+sz);
//         overlaps_.Insert(ols, sz);
//     };

//     auto work_func = [&check, this, lineToOl, block_size, generate_func, combine_func, &fname, &in](size_t id) {
//         std::vector<std::string> lines(block_size);
//         std::vector<Overlap> ols(block_size);
//         thread_local StringPool::UnsafeNameId ni(string_pool_);
//         while (true) {
//             size_t line_size = generate_func(lines);
            
//             if (line_size > 0) {
//                 size_t ol_size = 0;
//                 for (size_t i=0; i<line_size; ++i) {
//                     Overlap o;
//                     if (lineToOl(lines[i], o, ni)) {
//                         if (check(o)) {
//                             ols[ol_size++] = o;
//                         }
//                     } else {
//                         LOG(ERROR)("Failed to convert line to overlap \n   %s", lines[i].c_str());
//                     }
//                 }
//                 combine_func(ols, ol_size);;
//             } else {
//                 break;
//             }

//         }   
//         if (!in.IsEnd()) {
//             LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
//         }
//     };

//     if (in.Valid()) {
//         MultiThreadRun(thread_size, work_func);
//     } else {
//         LOG(ERROR)("Failed to load file: %s", fname.c_str());
//     }
// }

template<typename F, typename C>
void OverlapStore::LoadFileFast(const std::string &fname, F lineToOl, C check, size_t thread_size) {
    std::mutex mutex_gen;
    GzFileReader in(fname);

    const size_t block_size = 10000;
    if (thread_size > load_threads) thread_size = load_threads;
    
    auto generate_func = [&mutex_gen, &in](std::vector<char> &block) {
        std::lock_guard<std::mutex> lock(mutex_gen);
        return in.GetBlock(block, "\n", 1);
    };

    auto combine_func = [&check, this](std::vector<Overlap> &ols, size_t sz) {
        std::lock_guard<std::mutex> lock(mutex_overlaps_);
        //overlaps_.insert(overlaps_.end(), ols.begin(), ols.begin()+sz);
        overlaps_.Insert(ols, sz);
    };

    auto work_func = [&check, this, lineToOl, block_size, generate_func, combine_func, &fname, &in](size_t id) {
        std::vector<char> block(10000000);
        std::vector<Overlap> ols;
        ols.reserve(10000);

        StringPool::UnsafeNameId ni(string_pool_);
        size_t block_size = generate_func(block);

        std::string line;
        while (block_size > 0) {
            auto curr = block.begin();
            while (curr < block.begin()+block_size) {
                auto next = std::find(curr, block.begin()+block_size, '\n');
                //LOG(INFO)("pos: %d %d, %zd, %c %c", curr-block.begin(), next-block.begin(), block_size, *curr, *next);
                line = std::string(curr, next);
                //printf("line: %s\n", line.c_str());
                Overlap o;
                if (lineToOl(line, o, ni)) {
                    if (check(o)) {
                        ols.push_back(o);
                    }
                } else {
                    LOG(INFO)("line: %d %d, %zd, %c %c", curr-block.begin(), next-block.begin(), block_size, *curr, *next);
                    LOG(ERROR)("Failed to convert line to overlap \n   \"%s\"", line.c_str());
                }
                curr = next;
                if (curr < block.begin()+block_size && *curr == '\n')  curr++;
                // if (id == 0) {
                //     printf("curr %d/%zd\n", curr-block.begin(), block_size);
                // }
            }

            combine_func(ols, ols.size());
            ols.clear();
            block_size = generate_func(block);
        }

    
        if (!in.IsEnd()) {
            LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
        }
    };

    if (in.Valid()) {
        MultiThreadRun(thread_size, work_func);
    } else {
        LOG(ERROR)("Failed to load file: %s", fname.c_str());
    }
}

template<typename L, typename C>
void OverlapStore::SaveFile(const std::string &fname, L toLine, size_t thread_size, C check) const {
    std::ofstream of(fname);

    std::atomic<size_t> index { 0 };
    std::mutex mutex_of;

    auto work_func = [&](int tid) {
        std::ostringstream oss;
        StringPool::UnsafeNameId ni(string_pool_);

        auto flush = [&mutex_of, &of](std::ostringstream &oss) {
            std::lock_guard<std::mutex> lock(mutex_of);
            of << oss.str();
            oss.str("");
        };

        size_t curr = index.fetch_add(1);
        while (curr < Size()) {
            
            const auto &o = Get(curr);
            if (check(o)) {
                oss << (toLine)(o, ni) << "\n";
            }
            if (oss.tellp() > 100000000) {
                flush(oss);
            }
            curr = index.fetch_add(1);
        }
        flush(oss);

    };

    if (of.is_open()) {
        MultiThreadRun(thread_size, work_func);
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
    
    // std::ofstream of(fname);
    // StringPool::UnsafeNameId ni(string_pool_);
    // if (of.is_open()) {
    //     for (size_t i=0; i<Size(); ++i) {
    //         const auto &o = Get(i);
    //         if (check(o)) {
    //             of << (toLine)(o, ni) << "\n";
    //         }
    //     }
    // }
}


template<typename C>
void OverlapStore::Load(const std::string &fname, const std::string &type, size_t thread_size, C check) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "m4" || t == "m4.gz") {
        LoadFileMt(fname, &OverlapStore::FromM4Line, check, thread_size);
    } else if (t == "m4a" || t == "m4a.gz") {
        LoadFileMt(fname, &OverlapStore::FromM4aLine, check, thread_size);
    } else if (t == "paf" || t == "paf.gz") {
        LoadFileMt(fname, &OverlapStore::FromPafLine, check, thread_size);
    } else if (t == "txt") {
        LoadFileTxtMt(fname, check, thread_size);
    } else {
        LOG(ERROR)("Failed to recognize overlap files type: %s", t.c_str());
    }
}

template<typename C>
void OverlapStore::LoadFast(const std::string &fname, const std::string &type, size_t thread_size, C check) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "m4" || t == "m4.gz") {
        LoadFileFast(fname, &OverlapStore::FromM4Line, check, thread_size);
    } else if (t == "m4a" || t == "m4a.gz") {
        LoadFileFast(fname, &OverlapStore::FromM4aLine, check, thread_size);
    } else if (t == "paf" || t == "paf.gz") {
        LoadFileFast(fname, &OverlapStore::FromPafLine, check, thread_size);
    } else if (t == "txt") {
        LoadFileTxtFast(fname, check, thread_size);
    } else {
        LOG(ERROR)("Failed to recognize overlap files type: %s", t.c_str());
    }
}

template<typename C>
void OverlapStore::Save(const std::string &fname, const std::string &type, size_t thread_size, C check) const {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "m4") {        
        SaveFile(fname, &OverlapStore::ToM4Line, thread_size, check);
    } else if (t == "m4a") {
        SaveFile(fname, &OverlapStore::ToM4aLine, thread_size, check);
    } else if (t == "paf") {
        SaveFile(fname, &OverlapStore::ToPafLine, thread_size, check);
    } else {
        LOG(ERROR)("Failed to recognize overlap files type: %s", t.c_str());
    }
}

template<typename C>
void OverlapStore::LoadFileTxtMt(const std::string &fname, C check, size_t thread_size) {
    std::vector<std::string> files = GetLineFromFile(fname);

    for (const auto& f : files) {
        Load(f, "", thread_size, check);
    }    
}


template<typename C>
void OverlapStore::LoadFileTxtFast(const std::string &fname, C check, size_t thread_size) {
    std::vector<std::string> files = GetLineFromFile(fname);

    std::atomic<size_t> index {0};

    auto work_func = [&files, &index, this, check, thread_size](size_t id) {
        size_t curr = index.fetch_add(1);
        while (curr < files.size()) {
            LoadFast(files[curr], "", load_threads, check);
            curr = index.fetch_add(1);
        }
    };
    MultiThreadRun(std::max<size_t>(1, thread_size / load_threads), work_func);
}


} // namespace fsa {

#endif // FSA_OVERLAP_STORE_HPP
