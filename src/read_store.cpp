#include "read_store.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>

#include "file_io.hpp"
#include "./utils/logger.hpp"
#include "./utility.hpp"

namespace fsa {

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
    FastaReader reader(fname);
    LoadReader(reader, all, seqids);
}


void ReadStore::LoadFastq(const std::string &fname, bool all, const std::unordered_set<Seq::Id>& seqids) {
    FastqReader reader(fname);
    LoadReader(reader, all, seqids);
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

std::string ReadStore::DetectFileType(const std::string &fname) {
    const std::vector<std::array<std::string,2>> suffix = {
        {".fasta", "fasta"},
        {".fasta.gz", "fasta.gz"},
        {".fa", "fasta"},
        {".fa.gz", "fasta.gz"},
        {".fastq", "fastq"},
        {".fastq.gz", "fastq.gz"},
        {".fq", "fastq"},
        {".fq.gz", "fastq.gz"},
        {".fofn", "txt"},
        {".txt", "txt"},
    };

    for (const auto& sx : suffix) {
        // check the length first
        if (fname.length() >= sx[0].length() && fname.rfind(sx[0]) == fname.length() - sx[0].length()) {
            return sx[1];
        }
    }
    LOG(WARNING)("Unrecognized the suffix of %s", fname.c_str());
    return "fasta";
}

void ReadStore::LoadReader(SeqReader& reader, bool all, const std::unordered_set<Seq::Id>& seqids) {
    const std::string& fname = reader.GetFileName();

    size_t count = 0;
    if (reader.IsValid()) {
        SeqReader::Item item;

        while (reader.Next(item)) {
            AddItem(item, all, seqids);
            count++;
        }
        if (!reader.IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
    LOG(INFO)("Load %zd reads from file: %s", count, fname.c_str());
}

void ReadStore::AddItem(const SeqReader::Item &item, bool all, const std::unordered_set<Seq::Id>& seqids) {

    assert(!item.head.empty());
    // TODO not support muli-threads
    Seq::Id id = string_pool_.GetIdByStringUnsafe(item.head);
    assert(id >= 0);
    assert((size_t)id >= offset_);
    while (items_.size() <= id - offset_) {
        items_.push_back(Item());
    }
    
    if (all || seqids.find(id) != seqids.end()) items_[id - offset_].seq = item.seq;
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
