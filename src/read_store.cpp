#include "read_store.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>

#include "file_io.hpp"
#include "./utils/logger.hpp"

namespace fsa {

int ReadStore::GetIdByNameSafe(const std::string &name) {
    std::lock_guard<std::mutex> lock(mutex_);
    return GetIdByNameUnsafe(name);
}

int ReadStore::GetIdByNameUnsafe(const std::string &name) {
    int id = string_pool_.GetIdByStringUnsafe(name);
    while (items_.size() <= id) {
        items_.push_back(Item());
    }
    return id;
}

Seq::Id ReadStore::QueryIdByName(const std::string &name) const {
    return string_pool_.QueryIdByString(name);
}

const std::string& ReadStore::QueryNameById(int id) const {
    return string_pool_.QueryStringById(id);
}




void ReadStore::Load(const std::string &fname, const std::string &type, bool all, const std::unordered_set<Seq::Id>& seqids) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "fasta" || t == "fasta.gz") {
        LoadFasta(fname, all, seqids);
    } else if (t == "fastq" || t == "fastq.gz") {
        LoadFastq(fname, all, seqids);
    } else if (t == "fofn") {
        LoadFofn(fname, all, seqids);
    } else if (t == "txt") {
        LoadTxt(fname, all, seqids);
    } else {
        LOG(ERROR)("Failed to recognize read files type: %s", t.c_str());
    }
}


void ReadStore::LoadFasta(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) {
    std::unordered_set<Seq::Id> ids;

    FastaReader* reader = new FastaReader(fname);
    readers_.push_back(reader);

    if (reader->IsValid()) {
        SeqReader::Item item;
        while (reader->Next(item)) {
            assert(!item.head.empty());
            Seq::Id id = GetIdByNameUnsafe(item.head);
            Insert(id, item, reader, all || seqids.find(id) != seqids.end());
            ids.insert(id);
        }
        if (!reader->IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
    ids_in_file_[fname] = ids;
    LOG(INFO)("Load %zd reads from fasta file: %s", ids.size(), fname.c_str());
}


void ReadStore::LoadFastq(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) {
    std::unordered_set<Seq::Id> ids;

    FastqReader* reader = new FastqReader(fname);
    readers_.push_back(reader);
    
    if (reader->IsValid()) {
        SeqReader::Item item;

        while (reader->Next(item)) {
            assert(!item.head.empty());
            Seq::Id id = GetIdByNameSafe(item.head);
            Insert(id, item, reader, all || seqids.find(id) != seqids.end());
            ids.insert(id);
        }
        if (!reader->IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
    ids_in_file_[fname] = ids;
    LOG(INFO)("Load %zd reads from fastq file: %s", ids.size(), fname.c_str());
}

void ReadStore::LoadFofn(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) {
    std::ifstream in(fname);
    if (in.is_open()) {
        std::string line;
        while (std::getline(in, line)) {
            auto begin = std::find_if(line.begin(), line.end(), [](char a){return !::isspace(a); });
            if (begin != line.end()) {
                Load(line, "", all, seqids);
            }
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}

const std::unordered_set<Seq::Id>& ReadStore::IdsInFile(const std::string &fname) const {
    auto iter = ids_in_file_.find(fname);
    assert(iter != ids_in_file_.end());
    return iter->second;
}



std::string ReadStore::DetectFileType(const std::string &fname) {

    if (fname.size() >= 6 && fname.substr(fname.size()-6) == ".fasta") {
        return "fasta";
    } else if (fname.size() >= 9 && fname.substr(fname.size()-9) == ".fasta.gz") {
        return "fasta.gz";
    } else if (fname.size() >= 3 && fname.substr(fname.size()-3) == ".fa") {
        return "fasta";
    } else if (fname.size() >= 6 && fname.substr(fname.size()-6) == ".fa.gz") {
        return "fasta.gz";
    } else if (fname.size() >= 6 && fname.substr(fname.size()-6) == ".fastq") {
        return "fastq";
    } else if (fname.size() >= 9 && fname.substr(fname.size()-9) == ".fastq.gz") {
        return "fastq.gz";
    } else if (fname.size() >= 5 && fname.substr(fname.size()-5) == ".fofn") {
        return "fofn";
    } else if (fname.size() >= 4 && fname.substr(fname.size()-4) == ".txt") {
        return "txt";
    } else {
        LOG(WARNING)("Unrecognized the suffix of %s", fname.c_str());
        return "fasta";
    }
}

void ReadStore::Insert(Seq::Id id, const SeqReader::Item &item, SeqReader *reader, bool loadseq) {
    assert(id >= 0 && id < (int)string_pool_.Size());

    if (loadseq) items_[id].seq = item.seq;
    items_[id].id = item.id;
    items_[id].reader = reader;
}

void ReadStore::LoadItem(Item &item) const {
    if (item.seq.Size() == 0 && item.reader != nullptr) {
        SeqReader::Item i;
        auto r = item.reader->Get(item.id, i);
        if (r) {
            item.seq = i.seq;
        } else {
            LOG(ERROR)("Failed to load a read");
        }
    }
}

void ReadStore::SaveIdToName(const std::string &fname) const {
    string_pool_.Save(fname);
}


std::string ReadStore::GetSeq(const Seq::Tile& sa) {
    auto seq = GetSeq(sa.id).ToString();

    if (sa.strand == 0) {
        return seq->substr(sa.start, sa.end);
    }
    else {
        return Seq::ReverseComplement(seq->substr(sa.start, sa.end));
    }
}

} // namespace fsa {
