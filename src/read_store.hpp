#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "./utils/logger.hpp"
#include "sequence.hpp"
#include "seq_io.hpp"
#include "utils/string_pool.hpp"

namespace fsa {

class ReadStore {
public:
    ReadStore(StringPool &sp) : string_pool_(sp) { offset_ = string_pool_.Size(); }
    ReadStore(StringPool &sp, size_t offset) : string_pool_(sp), offset_(offset) { }
    ReadStore() : ReadStore(default_string_pool_) {}
    struct Item {
        Item(){}
        Item(const std::string &s) : seq(s) {}
        DnaSeq seq;
    };

    StringPool& GetStringPool() { return string_pool_; }
    const StringPool& GetStringPool() const { return string_pool_; }

    Seq::Id QueryIdByName(const std::string &name) const { return string_pool_.QueryIdByString(name);}
    const std::string& QueryNameById(Seq::Id id) const { return string_pool_.QueryStringById(id); }
    void SaveIdToName(const std::string& fname) const { string_pool_.Save(fname); }

    const DnaSeq& GetSeq(Seq::Id id) const { assert(size_t(id) >= offset_); return items_[id-offset_].seq;  }
    const DnaSeq& GetSeq(const std::string &name) { return GetSeq(QueryIdByName(name)); }

    std::string GetSeq(const Seq::Tile& sa);
    size_t GetSeqLength(Seq::Id id) const { return GetSeq(id).Size(); }

    std::array<size_t, 2> GetIdRange() const { return {offset_, offset_ + Size()}; }
    size_t Size() const { return items_.size(); }

    void Load(const std::string &fname, const std::string &type="", bool all=true, const std::unordered_set<Seq::Id>& seqids=std::unordered_set<Seq::Id>());
    void LoadFasta(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadFastq(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadFofn(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadTxt(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) { LoadFofn(fname, all, seqids); }
    
    static std::string DetectFileType(const std::string &fname);
protected:
    void LoadReader(SeqReader& reader, bool all, const std::unordered_set<Seq::Id>& seqids);
    void AddItem(const SeqReader::Item &item, bool all, const std::unordered_set<Seq::Id>& seqids);

protected:
    StringPool default_string_pool_;
    StringPool &string_pool_;   
    mutable std::vector<Item> items_;
    size_t offset_ { 0 };           ///< 序列数据要求连续读取，其序列名称ID号连续，offset为ID起点。 
};


template<typename F> 
void LoadReadFileFasta(const std::string &fname, F action) {
    FastaReader* reader = new FastaReader(fname);
    if (reader->IsValid()) {
        SeqReader::Item item;
        while (reader->Next(item)) {
            assert(!item.head.empty());
            action(item);
        }
        if (!reader->IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}

template<typename F> 
void LoadReadFileFastq(const std::string &fname, F action) {
    FastqReader* reader = new FastqReader(fname);
    if (reader->IsValid()) {
        SeqReader::Item item;

        while (reader->Next(item)) {
            assert(!item.head.empty());
            action(item);
        }
        if (!reader->IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}

template<typename F> 
void LoadReadFileTxt(const std::string &fname, F action) {
    std::ifstream in(fname);
    if (in.is_open()) {
        std::string line;
        while (std::getline(in, line)) {
            auto begin = std::find_if(line.begin(), line.end(), [](char a){return !::isspace(a); });
            if (begin != line.end()) {
                LoadReadFile(line, "", action);
            }
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}


template<typename F> 
void LoadReadFile(const std::string &fname, const std::string &type, F action) {
    std::string t = type != "" ? type : ReadStore::DetectFileType(fname);
    if (t == "fasta" || t == "fasta.gz") {
        LoadReadFileFasta(fname, action);
    } else if (t == "fastq" || t == "fastq.gz") {
        LoadReadFileFastq(fname, action);
    } else if (t == "fofn") {
        LoadReadFileTxt(fname, action);
    } else if (t == "txt") {
        LoadReadFileTxt(fname, action);
    } else {
        LOG(ERROR)("Failed to recognize read files type: %s", t.c_str());
    }
}

} // namespace fsa {

