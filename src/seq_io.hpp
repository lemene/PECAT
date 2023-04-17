/**
 * @brief read/write sequence files (fasta/fastq).
 */

#pragma once

#include "file_io.hpp"

namespace fsa {

class SeqReader {
public:
    using ItemId = int64_t;
    struct Item {
        std::string head;
        std::string sub_head;
        std::string seq;
        std::string quality;
        ItemId id;            // location in file
    };
public:
    SeqReader(const std::string &fname) : fname_(fname) {}
    virtual ~SeqReader() {}

    const std::string& GetFileName() const { return fname_; }

    virtual bool IsValid() const = 0;
    virtual bool IsFileEnd() const = 0;

    virtual bool Next(Item &item) = 0;
    virtual bool Get(ItemId id, Item &item) = 0;

    static std::string NextNonEmptyLine(std::ifstream &in);
protected:
    std::string fname_;
};

class SeqWriter {
public:
    virtual ~SeqWriter() {}
    virtual void Write(const SeqReader::Item& item)=0;
};


class FastaReader : public SeqReader {
public:
    FastaReader(const std::string &fname) : SeqReader(fname), in_(fname) {}
    virtual ~FastaReader() {}

    virtual bool IsValid() const { return in_.Valid(); }
    virtual bool IsFileEnd() const { return in_.IsEnd(); }

    virtual bool Next(Item &item);
    virtual bool Get(ItemId id, Item &item) { Seek(id); return Next(item);}

protected:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq);

    ItemId Tell() { return in_.Tell(); }
    void Seek(ItemId id) { in_.Seek(id); }

protected:
    GzFileReader in_;
};

class FastqReader : public SeqReader {
public:
    FastqReader(const std::string &fname) : SeqReader(fname), in_(fname) { }
    virtual ~FastqReader() {}

    virtual bool IsValid() const { return in_.Valid(); }
    virtual bool IsFileEnd() const  { return in_.IsEnd(); }

    virtual bool Next(Item &item);
    virtual bool Get(ItemId id, Item &item) { Seek(id); return Next(item); }

protected:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq) { seq = in_.GetStrippedLine(); return true; }
    bool GetHead1() { std::string line = in_.GetStrippedLine(); return line[0] == '+'; }
    bool GetQuality(std::string &qua) { qua = in_.GetStrippedLine(); return true;  }

    ItemId Tell() { return in_.Tell(); }
    void Seek(ItemId id) { in_.Seek(id); }
protected:
    GzFileReader in_;
};


class FastaWriter : public SeqWriter {
public:
    FastaWriter(const std::string &fname) : out_(fname) {
    }
    virtual ~FastaWriter() {}
    
    void Write(const SeqReader::Item& item) {
        out_.Write(">");
        out_.Write(item.head);
        out_.Write(" ");
        out_.Write(item.sub_head);
        out_.Write("\n");
        out_.Write(item.seq);
        out_.Write("\n");
    }

protected:
    GzFileWriter out_;
};

class FastqWriter : public SeqWriter {
public:
    FastqWriter(const std::string &fname) : out_(fname) {
    }
    virtual ~FastqWriter() {}
    
    void Write(const SeqReader::Item& item) {
        out_.Write("@");
        out_.Write(item.head);
        out_.Write(" ");
        out_.Write(item.sub_head);
        out_.Write("\n");
        out_.Write(item.seq);
        out_.Write("\n");
        out_.Write("+");
        out_.Write(item.head);
        out_.Write(" ");
        out_.Write(item.sub_head);
        out_.Write("\n");
        out_.Write(item.quality);
        out_.Write("\n");
    }

protected:
    GzFileWriter out_;
};

} // namespace fsa {

