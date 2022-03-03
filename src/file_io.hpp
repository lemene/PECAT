#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <zlib.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <mutex>
#include <memory>

namespace fsa {

class Reader {
public:
    virtual ~Reader() {};
    
    virtual size_t GetLines(std::vector<std::string> &lines) = 0;
    virtual size_t GetBlock0(char* block, size_t max_size) = 0;
    virtual size_t GetBlock(std::vector<char> &block, const std::string& sentinel="", int off=0);

protected:
    std::vector<char> buf;
    size_t bufsize { 0 };
};

class StdioReader : public Reader {
public:
    StdioReader() {std::ios::sync_with_stdio(false);}
    virtual ~StdioReader() {}

    virtual size_t GetLines(std::vector<std::string> &lines);
    virtual size_t GetBlock0(char* block, size_t max_size);
};

class GzFileReader : public Reader {
public:
    GzFileReader(const std::string &fname);
    ~GzFileReader() { if (in_ != nullptr) gzclose(in_); }
    bool Valid() const { return in_ != nullptr; }
    bool IsEnd() const { return gzeof(in_); }
    std::string GetStrippedLine();
    std::string GetNoEmptyLine();
    const std::string &QueryNoEmptyLine() {
        if (!HasNoEmptyLine()) {
            int64_t pos = Tell();   // TODO record pos
            next_line_ = GetNoEmptyLine();
            next_line_pos_ = pos;
        }
        return next_line_;
    }
    void ConsumeStrippedLine() { next_line_pos_ = -1; }     
    bool HasNoEmptyLine() const { return next_line_pos_ >= 0; }

    std::string GetLine();
    bool GetLine(std::string &s);
    size_t GetLines(std::vector<std::string> &lines);

    int64_t Tell() { return !HasNoEmptyLine() ? gztell(in_) : next_line_pos_; }
    void Seek(int64_t pos) { 
        next_line_.clear(); 
        next_line_pos_ = -1;   
        gzseek(in_, pos, SEEK_SET);
    }

    virtual size_t GetBlock0(char* block, size_t max_size);
    
protected:
    gzFile in_ { nullptr };

    std::string next_line_ { "" };
    int64_t next_line_pos_ {-1} ;
};

class MultiFileReader : public Reader {
public:
    MultiFileReader(const std::vector<std::string> &ifnames) : ifnames_(ifnames) {}
    virtual ~MultiFileReader() {};
    
    virtual size_t GetLines(std::vector<std::string> &lines);
    virtual size_t GetBlock0(char* block, size_t max_size);

    std::vector<std::string> ifnames_;
    size_t file_index_ { 0 };
    std::shared_ptr<GzFileReader> ifile_ { nullptr };
};

class Writer {
public:
    virtual ~Writer() {};
    virtual void Write(const std::string &str) = 0;
    virtual void Flush() = 0;
};

inline Writer& operator << (Writer& writer, const std::string &str) {
    writer.Write(str);
    return writer;
}

template<typename T>
inline Writer& operator << (Writer& writer, T s) {
    std::ostringstream oss;
    oss << s;
    writer << oss.str();
    return writer;
}

class StdioWriter : public Writer {
    
public:
    StdioWriter() {}
    virtual ~StdioWriter() {}
    virtual void Write(const std::string &str) { printf("%s", str.c_str()); }//std::cout << str <<std::flush; }
    virtual void Flush() { }
    
};

class GzFileWriter : public Writer {
public:
    GzFileWriter(const std::string &fname, bool plain) : plain_(plain) {
        /*

        */
    }
    bool Valid() const { return (plain_ && out_plain_ != nullptr) || (!plain_ && out_compress_ != nullptr); }

    GzFileWriter(const std::string &fname) {
        plain_ = !(fname.size() >= 3 && fname.substr(fname.size()-3) == ".gz");
        if (plain_) {
            out_plain_ = fopen(fname.c_str(), "wb");
        } else {
            out_compress_ = gzopen(fname.c_str(), "wb");
        }
    }

    void Write(const std::string &str) {
        if (plain_) {
            fputs(str.c_str(), out_plain_);
        } else {
            gzputs(out_compress_, str.c_str());
        }
    }
    virtual void Flush() {/* TOD8 */}
    virtual void Flush(std::ostringstream &oss) {
        Write(oss.str());
        oss.str("");
    }

    ~GzFileWriter() { 
        if (out_compress_ != nullptr) gzclose(out_compress_); 
        if (out_plain_ != nullptr) fclose(out_plain_); 
    }
    
protected:
    gzFile out_compress_ {nullptr};
    FILE* out_plain_ {nullptr};
    bool plain_ ;

};


std::vector<std::string> GetLineFromFile(const std::string& fname);



class LineInBlock {
public:
    LineInBlock(Reader &reader, size_t bsize, std::mutex *m = nullptr) 
     : reader_(reader), block_(bsize), mutex_(m) {

        block_size_ = 0;
        block_index_ = block_.begin();
    }

    bool GetLine(std::string& line);

    Reader& reader_;
    std::vector<char> block_;
    size_t block_size_;
    decltype(block_)::iterator block_index_;
    std::mutex *mutex_ { nullptr };
};

} // namespace fsa {


