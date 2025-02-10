#pragma once

#include <array>
#include <vector>
#include <list>
#include <deque>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <htslib/sam.h>


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


    std::string QueryNameById(int id) const {
        return string_pool_.QueryStringById(id);
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
    void Group(std::unordered_map<int, std::unordered_map<int, const Overlap*>>& groups, 
               std::unordered_map<int, std::unordered_map<int, std::vector<const Overlap*>>> &dup_groups,
               bool (*better)(const Overlap&, const Overlap&),
               const std::unordered_set<int>& keys, 
               size_t thread_size) const;
    void GroupTarget(std::unordered_map<Seq::Id, std::unordered_map<Seq::Id, std::vector<const Overlap*>>> &groups, size_t threads) const;
    void GroupQuery(std::unordered_map<Seq::Id, std::unordered_map<Seq::Id, std::vector<const Overlap*>>> &groups, size_t threads) const;


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

    template<typename C>
    void LoadFileBam(const std::string &fname, C check, size_t thread_size);
    template<typename C>
    void LoadFileBamFast(const std::string &fname, C check, size_t thread_size);


    template<typename C>
    void LoadFileBam(const std::string &fname, C check, size_t thread_size, const std::vector<Bed>& beds);

    // TODO ugly, only for sam
    void PreLoad(Reader &reader, std::vector<std::string>& done);
    void PreLoad(const std::string &fname);
    void AfterLoad(const std::string &fname);

    static int FromM4Line(const std::string &line, Overlap &o, StringPool::NameId& ni);
    static int FromM4aLine(const std::string &line, Overlap &o, StringPool::NameId& ni);
    static int FromPafLine(const std::string &line, Overlap &o, StringPool::NameId& ni);
    static int FromSamLine(const std::string &line, Overlap &o, StringPool::NameId& ni);

    static int FromM4LineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen) {
        replen = 0;
        return FromM4Line(line, o, ni);
    }
    static int FromM4aLineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen) {
        replen = 0;
        return FromM4aLine(line, o, ni);
    }

    static int FromPafLineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen);

    static int FromSamLineEx(const std::string &line, Overlap &o, StringPool::NameId& ni, int &replen) {
        replen = 0;
        return FromSamLine(line, o, ni);
    }

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

    size_t load_threads = 10;
    static std::unordered_map<int, int> loading_infos_; // TODO for loading sam file

};

template<typename F, typename C>
void OverlapStore::LoadFileMt(const std::string &fname, F lineToOl, C check, size_t thread_size) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;
    GzFileReader in(fname);

    if (thread_size > load_threads) thread_size = load_threads;

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
            //overlaps_.Insert(v_ols, v_ols.size());
            for (size_t i=0; i<sz; ++i) {
                auto& o = ols[i];
                if (check(o)) {
                    overlaps_.Add(o);
                }
            }
        }
    };

    auto work_func = [&check, this, lineToOl, combine_func, &fname, &mutex_gen, &in](size_t id) {
        const size_t max_overlap_size = 1000000;
        const size_t max_block_size = 10000000;

        std::vector<Overlap> ols(max_overlap_size);
        size_t ol_size = 0;

        LineInBlock line_in_block(in, max_block_size, &mutex_gen);
        thread_local StringPool::TempNameId ni;
        std::string line;
        for (bool is_valid = line_in_block.GetLine(line); is_valid; is_valid = line_in_block.GetLine(line)) {
            if (line.size() < 1 || line[0] == '#') continue;
            Overlap o;
            auto r = lineToOl(line, o, ni);
            if (r > 0) {
                ols[ol_size++] = o;
            } else if (r < 0) {
                LOG(ERROR)("Failed to convert line to overlap \n    %s\n    %s", line.c_str(), fname.c_str());
            } else {
                // r == 0 pass
            }
            if (ol_size >= max_overlap_size) {
                combine_func(ols, ol_size, ni);
                ol_size = 0;
            }
        }
        combine_func(ols, ol_size, ni);
        ol_size = 0;

        if (!in.IsEnd()) {
            LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
        }
    };

    if (in.Valid()) {
        PreLoad(fname);
        MultiThreadRun(thread_size, work_func);
        AfterLoad(fname);
    } else {
        LOG(ERROR)("Failed to load file: %s", fname.c_str());
    }
}

template<typename F, typename C>
void OverlapStore::LoadFileFast(const std::string &fname, F lineToOl, C check, size_t thread_size) {
    std::mutex mutex_gen;
    GzFileReader in(fname);

    if (thread_size > load_threads) thread_size = load_threads;

    auto combine_func = [&check, this](std::vector<Overlap> &ols, size_t sz) {
        std::lock_guard<std::mutex> lock(mutex_overlaps_);
        overlaps_.Insert(ols, sz);
    };

    auto work_func = [&check, this, lineToOl, &mutex_gen, combine_func, &fname, &in](size_t id) {
        const size_t max_block_size = 10000000;
        const size_t max_overlap_size = 100000;

        std::vector<Overlap> ols;
        ols.reserve(max_overlap_size);

        StringPool::UnsafeNameId ni(string_pool_);
        LineInBlock line_in_block(in, max_block_size, &mutex_gen);

        std::string line;        
        for (bool is_valid = line_in_block.GetLine(line); is_valid; is_valid = line_in_block.GetLine(line)) {
            if (line.size() < 1 || line[0] == '#') continue;

            Overlap o;
            auto r = lineToOl(line, o, ni);
            if (r > 0) {
                if (check(o)) {
                    ols.push_back(o);
                }
            } else if (r < 0) {
                LOG(ERROR)("Failed to convert line to overlap \n   \"%s\"", line.c_str());
            } else {
                // r == 0 pass
            }

            if (ols.size() >= max_overlap_size) {
                combine_func(ols, ols.size());
                ols.clear();
            }
        }
        if (ols.size() > 0) {
            combine_func(ols, ols.size());
            ols.clear();
        }

        if (!in.IsEnd()) {
            LOG(WARNING)("No all overlaps in file are loaded: %s", fname.c_str());
        }
    };

    if (in.Valid()) {
        PreLoad(fname);
        MultiThreadRun(thread_size, work_func);
        AfterLoad(fname);
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
    } else if (t == "sam" || t == "sam.gz") {
        LoadFileMt(fname, &OverlapStore::FromSamLine, check, thread_size);
    } else if (t == "txt") {
        LoadFileTxtMt(fname, check, thread_size);
    } else if (t == "bam") {
        LoadFileBam(fname, check, thread_size);
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
    } else if (t == "sam" || t == "sam.gz") {
        LoadFileFast(fname, &OverlapStore::FromSamLine, check, thread_size);
    } else if (t == "txt") {
        LoadFileTxtFast(fname, check, thread_size);
    } else if (t == "bam") {
        LoadFileBamFast(fname, check, thread_size);
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

template<typename C>
void OverlapStore::LoadFileBam(const std::string &fname, C check, size_t thread_size) {

    std::mutex mutex;
    std::atomic<size_t> index { 0 };

    auto combine_func = [&mutex, this](std::vector<Overlap> &ols, size_t sz, StringPool::TempNameId &ni) {
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
            std::lock_guard<std::mutex> lock(mutex);
            overlaps_.Insert(ols, sz);
        }
    };
    
    auto work_func = [combine_func, &index, check, this, &fname](int tid) {

        samFile *in = hts_open(fname.c_str(), "r");
        assert(in != nullptr);
        
        bam_hdr_t *header = sam_hdr_read(in);
        assert(header != nullptr);
        auto idx = sam_index_load2(in, fname.c_str(), (fname+".bai").c_str());
        assert(idx != nullptr);

        const size_t max_overlap_size = 100000;
        std::vector<Overlap> ols;
        ols.reserve(max_overlap_size);
        thread_local StringPool::TempNameId ni;

        for (size_t tgtid = index.fetch_add(1); tgtid < header->n_targets; tgtid = index.fetch_add(1)) {

            hts_itr_t* itr = sam_itr_queryi(idx, tgtid, 0, header->target_len[tgtid]);
            bam1_t *b = bam_init1();
            
            while (sam_itr_next(in, itr, b) >= 0) {
                Overlap ol;

                //if ((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FSUPPLEMENTARY)) continue; 
                ol.a_.id = ni.GetIdByName(bam_get_qname(b));
                int tid = b->core.tid;
                if (tid >= 0 && tid < header->n_targets) {
                    ol.b_.id = ni.GetIdByName(header->target_name[tid]);
                } else {
                    continue;
                }

                ol.a_.len = b->core.l_qseq;
                ol.a_.strand = bam_is_rev(b) ? 1 : 0;

                ol.b_.len = header->target_len[tgtid];
                ol.b_.strand = 0;
                ol.b_.start = b->core.pos;

                uint32_t *cigar = bam_get_cigar(b);
                int rpos = 0;
                int qpos = 0;
                int match = 0;
                int clip = 0;
                for(int i=0; i < b->core.n_cigar;++i){
                    int icigar = cigar[i];
                    int n = bam_cigar_oplen(icigar);
                    char t = bam_cigar_opchr(icigar);
                    ol.detail_.push_back({n, t});
                    switch (t) {
                        
                    case '=':
                    case 'M':
                        match += n;
                        rpos += n;
                        qpos += n;
                        break;
                    case 'X':
                        rpos += n;
                        qpos += n;
                        break;
                    case 'I':
                        qpos += n;
                        break;
                    case 'D':
                        rpos += n;
                        break;
                    case 'S':
                    case 'H':
                        clip += n;
                        break;
                    default:
                        //LOG(ERROR)("Not support %c", t);
                        break;
                    }
                }
                ol.b_.end = ol.b_.start + rpos;
                ol.identity_ = match * 2.0 / (qpos + rpos);

                ol.a_.len = clip+qpos;

                assert (b->core.n_cigar >= 1);
                char cigar0 = bam_cigar_opchr(cigar[0]);
                char cigar1 = bam_cigar_opchr(cigar[b->core.n_cigar-1]);
    
                if ((cigar0 == 'S' || cigar0 == 'H') && ol.a_.strand == 0) {
                    ol.a_.start = bam_cigar_oplen(cigar[0]);
                    ol.a_.end = ol.a_.start + qpos;
                } else if ((cigar1 == 'S' || cigar1 == 'H')  && ol.a_.strand == 1) {
                    ol.a_.start = bam_cigar_oplen(cigar[b->core.n_cigar-1]);
                    ol.a_.end = ol.a_.start + qpos;
                } else {
                    ol.a_.start = 0;
                    ol.a_.end = ol.a_.start + qpos;
                }
                assert(0 <= ol.a_.start && ol.a_.start < ol.a_.end && ol.a_.end <= ol.a_.len);
                assert(0 <= ol.b_.start && ol.b_.start < ol.b_.end && ol.a_.end <= ol.b_.len);

                if (check(ol)) {
                    ols.push_back(ol);
                    if (ols.size() >= max_overlap_size) {
                        combine_func(ols, ols.size(), ni);
                        ols.clear();
                    }
                }

            }
            bam_destroy1(b);


        }
        if (ols.size() > 0) {
            combine_func(ols, ols.size(), ni);
            ols.clear();
        }

        bam_hdr_destroy(header);
        hts_close(in);
    };

        
    MultiThreadRun(thread_size, work_func);
}

template<typename C>
void OverlapStore::LoadFileBam(const std::string &fname, C check, size_t thread_size, const std::vector<Bed> &beds) {

    std::mutex mutex;
    std::atomic<size_t> index { 0 };

    auto combine_func = [&mutex, this](std::vector<Overlap> &ols, size_t sz, StringPool::TempNameId &ni) {
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
            std::lock_guard<std::mutex> lock(mutex);
            overlaps_.Insert(ols, sz);
        }
    };
    
    auto work_func = [combine_func, &index, check, this, &fname, &beds](int tid) {

        samFile *in = hts_open(fname.c_str(), "r");
        assert(in != nullptr);
        
        bam_hdr_t *header = sam_hdr_read(in);
        assert(header != nullptr);
        auto idx = sam_index_load2(in, fname.c_str(), (fname+".bai").c_str());
        assert(idx != nullptr);

        const size_t max_overlap_size = 100000;
        std::vector<Overlap> ols;
        ols.reserve(max_overlap_size);
        thread_local StringPool::TempNameId ni;

        for (size_t i = index.fetch_add(1); i < beds.size(); i = index.fetch_add(1)) {
            size_t tgtid = sam_hdr_name2tid(header, beds[i].target.c_str());
            hts_itr_t* itr = sam_itr_queryi(idx, tgtid, beds[i].start, beds[i].end);
            bam1_t *b = bam_init1();
            
            while (sam_itr_next(in, itr, b) >= 0) {
                Overlap ol;

                //if ((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FSUPPLEMENTARY)) continue; 
                ol.a_.id = ni.GetIdByName(bam_get_qname(b));
                int tid = b->core.tid;
                if (tid >= 0 && tid < header->n_targets) {
                    ol.b_.id = ni.GetIdByName(header->target_name[tid]);
                } else {
                    continue;
                }

                ol.a_.len = b->core.l_qseq;
                ol.a_.strand = bam_is_rev(b) ? 1 : 0;

                ol.b_.len = header->target_len[tgtid];
                ol.b_.strand = 0;
                ol.b_.start = b->core.pos;

                uint32_t *cigar = bam_get_cigar(b);
                int rpos = 0;
                int qpos = 0;
                int match = 0;
                int clip = 0;
                for(int i=0; i < b->core.n_cigar;++i){
                    int icigar = cigar[i];
                    int n = bam_cigar_oplen(icigar);
                    char t = bam_cigar_opchr(icigar);
                    ol.detail_.push_back({n, t});
                    switch (t) {
                        
                    case '=':
                    case 'M':
                        match += n;
                        rpos += n;
                        qpos += n;
                        break;
                    case 'X':
                        rpos += n;
                        qpos += n;
                        break;
                    case 'I':
                        qpos += n;
                        break;
                    case 'D':
                        rpos += n;
                        break;
                    case 'S':
                    case 'H':
                        clip += n;
                        break;
                    default:
                        //LOG(ERROR)("Not support %c", t);
                        break;
                    }
                }
                ol.b_.end = ol.b_.start + rpos;
                ol.identity_ = match * 2.0 / (qpos + rpos);

                ol.a_.len = clip+qpos;

                assert (b->core.n_cigar >= 1);
                char cigar0 = bam_cigar_opchr(cigar[0]);
                char cigar1 = bam_cigar_opchr(cigar[b->core.n_cigar-1]);
    
                if ((cigar0 == 'S' || cigar0 == 'H') && ol.a_.strand == 0) {
                    ol.a_.start = bam_cigar_oplen(cigar[0]);
                    ol.a_.end = ol.a_.start + qpos;
                } else if ((cigar1 == 'S' || cigar1 == 'H')  && ol.a_.strand == 1) {
                    ol.a_.start = bam_cigar_oplen(cigar[b->core.n_cigar-1]);
                    ol.a_.end = ol.a_.start + qpos;
                } else {
                    ol.a_.start = 0;
                    ol.a_.end = ol.a_.start + qpos;
                }
                assert(0 <= ol.a_.start && ol.a_.start < ol.a_.end && ol.a_.end <= ol.a_.len);
                assert(0 <= ol.b_.start && ol.b_.start < ol.b_.end && ol.a_.end <= ol.b_.len);

                if (check(ol)) {
                    ols.push_back(ol);
                    if (ols.size() >= max_overlap_size) {
                        combine_func(ols, ols.size(), ni);
                        ols.clear();
                    }
                }

            }
            bam_destroy1(b);


        }
        if (ols.size() > 0) {
            combine_func(ols, ols.size(), ni);
            ols.clear();
        }

        bam_hdr_destroy(header);
        hts_close(in);
    };

        
    MultiThreadRun(thread_size, work_func);
}


template<typename C>
void OverlapStore::LoadFileBamFast(const std::string &fname, C check, size_t thread_size) {

    samFile *in = hts_open(fname.c_str(), "r");
    if (in != nullptr) {
    
        std::unordered_map<Seq::Id, size_t> reflen;
        bam_hdr_t *header = sam_hdr_read(in);
        if (header != nullptr) {
            for (int i = 0; i < header->n_targets; i++) {
                auto id = string_pool_.GetIdByString(header->target_name[i]);
                reflen[id] = header->target_len[i];
            }
        } else {
            LOG(ERROR)("Error reading BAM header: %s", fname.c_str());
        }

        auto idx = sam_index_load2(in, fname.c_str(), (fname+".bai").c_str());


        std::mutex mutex;
        std::mutex mutex_gen;
        std::atomic<size_t> index { 0 };
        auto gen_func = [&mutex_gen,&idx, &header](size_t tgtid) {
            std::lock_guard<std::mutex> lock(mutex_gen);
            hts_itr_t* itr = sam_itr_queryi(idx, tgtid, 0, header->target_len[tgtid]);
            return itr;

        };
        auto combine_func = [&mutex, this](std::vector<Overlap> &ols, size_t sz) {
            std::lock_guard<std::mutex> lock(mutex);
            overlaps_.Insert(ols, sz);
        };
        
        auto work_func = [gen_func, combine_func, &idx, &header, &in, &index, check, this](int tid) {
            
            const size_t max_overlap_size = 100000;
            std::vector<Overlap> ols;
            ols.reserve(max_overlap_size);

            for (size_t tgtid = index.fetch_add(1); tgtid < header->n_targets; tgtid = index.fetch_add(1)) {
                

                //hts_itr_t* itr = sam_itr_queryi(idx, tgtid, 0, header->target_len[tgtid]);
                hts_itr_t* itr = gen_func(tgtid);
                bam1_t *b = bam_init1();
                
                while (sam_itr_next(in, itr, b) >= 0) {
                    Overlap ol;

                    //if ((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FSUPPLEMENTARY)) continue; 
                    ol.a_.id = string_pool_.GetIdByString(bam_get_qname(b));
                    int tid = b->core.tid;
                    if (tid >= 0 && tid < header->n_targets) {
                        ol.b_.id = string_pool_.GetIdByString(header->target_name[tid]);
                    } else {
                        continue;
                    }

                    ol.a_.len = b->core.l_qseq;
                    ol.a_.strand = bam_is_rev(b) ? 1 : 0;

                    ol.b_.len = header->target_len[tgtid];
                    ol.b_.strand = 0;
                    ol.b_.start = b->core.pos;

                    uint32_t *cigar = bam_get_cigar(b);
                    int rpos = 0;
                    int qpos = 0;
                    int match = 0;
                    int clip = 0;
                    for(int i=0; i < b->core.n_cigar;++i){
                        int icigar = cigar[i];
                        int n = bam_cigar_oplen(icigar);
                        char t = bam_cigar_opchr(icigar);
                        ol.detail_.push_back({n, t});
                        switch (t) {
                            
                        case '=':
                        case 'M':
                            match += n;
                            rpos += n;
                            qpos += n;
                            break;
                        case 'X':
                            rpos += n;
                            qpos += n;
                            break;
                        case 'I':
                            qpos += n;
                            break;
                        case 'D':
                            rpos += n;
                            break;
                        case 'S':
                        case 'H':
                            clip += n;
                            break;
                        default:
                            //LOG(ERROR)("Not support %c", t);
                            break;
                        }
                    }
                    ol.b_.end = ol.b_.start + rpos;
                    ol.identity_ = match * 2.0 / (qpos + rpos);

                    ol.a_.len = clip+qpos;

                    assert (b->core.n_cigar >= 1);
                    char cigar0 = bam_cigar_opchr(cigar[0]);
                    char cigar1 = bam_cigar_opchr(cigar[b->core.n_cigar-1]);
        
                    if ((cigar0 == 'S' || cigar0 == 'H') && ol.a_.strand == 0) {
                        ol.a_.start = bam_cigar_oplen(cigar[0]);
                        ol.a_.end = ol.a_.start + qpos;
                    } else if ((cigar1 == 'S' || cigar1 == 'H')  && ol.a_.strand == 1) {
                        ol.a_.start = bam_cigar_oplen(cigar[b->core.n_cigar-1]);
                        ol.a_.end = ol.a_.start + qpos;
                    } else {
                        ol.a_.start = 0;
                        ol.a_.end = ol.a_.start + qpos;
                    }
                    assert(0 <= ol.a_.start && ol.a_.start < ol.a_.end && ol.a_.end <= ol.a_.len);
                    assert(0 <= ol.b_.start && ol.b_.start < ol.b_.end && ol.a_.end <= ol.b_.len);

                    if (check(ol)) {
                        ols.push_back(ol);
                        if (ols.size() >= max_overlap_size) {
                            combine_func(ols, ols.size());
                            ols.clear();
                        }
                    }

                }
                bam_destroy1(b);


            }
            if (ols.size() > 0) {
                combine_func(ols, ols.size());
                ols.clear();
            }

        };

            
            MultiThreadRun(thread_size, work_func);
        bam_hdr_destroy(header);
        hts_close(in);
    }

}

class OverlapGrouper {
public:
    OverlapGrouper(OverlapStore& ols) : ol_store_(ols) { }

    void BuildIndex(size_t thread_size, const std::unordered_set<int>& read_ids);
    void BuildIndex(size_t thread_size, const std::unordered_set<int>& read_ids, int (*compare)(const Overlap&, const Overlap&));
    

    void ClusterReads(const std::vector<Seq::Id>& reads);

    class Group {
    public:
        Group(Seq::Id i) : id(i) {}
        bool Empty() const { return ols.size() == 0; }
        size_t Size() const { return index.size(); }
        size_t Size(size_t i) const { return index[i][1] - index[i][0]; }
        const Overlap* Get(size_t i, size_t j) const { return ols[index[i][0]+j]; }

        void Sort(double opt_ohwt);
        std::vector<double> GetWeight(double opt_ohwt);

        Seq::Id id;
        std::vector<const Overlap*> ols;
        std::vector<std::array<size_t, 2>> index;
    };
    Group Get(int id) const;
    std::vector<const Overlap*> GetRelatedOverlaps(Seq::Id id) const;

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

