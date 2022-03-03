#ifndef FSA_READ_STORE_HPP
#define FSA_READ_STORE_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <algorithm>
#include "./utils/logger.hpp"
#include "sequence.hpp"
#include "seq_io.hpp"
#include "utils/string_pool.hpp"

namespace fsa {

class ReadStore {
public:
    ReadStore(StringPool &sp) : string_pool_(sp) {}
    ReadStore() : ReadStore(default_string_pool_) {}
    struct Item {
        Item(){}
        Item(const std::string &s, SeqReader::ItemId i, SeqReader* r) : seq(s), id(i), reader(r) {}
        DnaSeq seq;
        SeqReader::ItemId id {-1};
        SeqReader* reader {nullptr};
    };

    Seq::Id QueryIdByName(const std::string &name) const;
    
    Seq::Id GetIdByNameSafe(const std::string &name);
    Seq::Id GetIdByNameUnsafe(const std::string &name);
    const std::string& QueryNameById(Seq::Id id) const;
    StringPool& GetStringPool() { return string_pool_; }
    const StringPool& GetStringPool() const { return string_pool_; }


    const DnaSeq& GetSeq(Seq::Id id) const { LoadItem(items_[id]); return items_[id].seq;  }
    const DnaSeq& GetSeq(const std::string &name) { return GetSeq(QueryIdByName(name)); }

    std::string GetSeq(const Seq::Tile& sa);
    size_t GetSeqLength(Seq::Id id) const { return GetSeq(id).Size(); }
    const Item& GetItem(Seq::Id id) const { return items_[id]; }
    void SaveIdToName(const std::string& fname) const;
    std::array<Seq::Id, 2> GetIdRange() const { return std::array<Seq::Id, 2>{(Seq::Id)0, (Seq::Id)string_pool_.Size()}; }
    size_t Size() const { return string_pool_.Size(); }

    void Load(const std::string &fname, const std::string &type="", bool all=true, const std::unordered_set<Seq::Id>& seqids=std::unordered_set<Seq::Id>());
    void LoadFasta(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadFastq(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadFofn(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids);
    void LoadTxt(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) { LoadFofn(fname, all, seqids); }
    void LoadItem(Item &item) const;
    
    const std::unordered_set<Seq::Id>& IdsInFile(const std::string &fname) const;

    static std::string DetectFileType(const std::string &fname);
protected:
    Seq::Id Insert(std::string &&name, std::string &&seq);
    Seq::Id Insert(const SeqReader::Item &item, SeqReader *reader, int mode);
    void Insert(Seq::Id id, const SeqReader::Item &item, SeqReader *reader, bool loadseq);

protected:
    std::mutex mutex_;
    StringPool default_string_pool_;
    StringPool &string_pool_;   
    mutable std::vector<Item> items_;
    std::unordered_map<std::string, std::unordered_set<Seq::Id>> ids_in_file_;
    std::vector<SeqReader*> readers_;
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
    
#endif // FSA_READ_STORE_HPP
