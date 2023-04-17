#include "file_io.hpp"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include "./utils/logger.hpp"
#include "utility.hpp"

namespace fsa {


std::string StripString(const std::string &line) {
    
    std::string::size_type s = 0;
    while (s< line.size() && ::isspace(line[s])) s++;
        
    std::string::size_type e = line.size();
    while (e > s && ::isspace(line[e-1])) e--;
        
    return line.substr(s, e-s);
        
}


size_t Reader::GetBlock(std::vector<char> &block, const std::string& sentinel, int off) { 
    size_t block_index = 0;
    if (bufsize > 0) {
        assert(bufsize < block.size());
        std::copy(buf.begin(), buf.begin()+bufsize, block.begin());
        block_index = bufsize;
        bufsize = 0;
    }

    block_index += GetBlock0(&block[0]+block_index, block.size()-block_index);

    if (sentinel.length() > 0 && block_index == block.size()) {
        assert((size_t)off <= sentinel.size());

        size_t end = 0;
        for (size_t i = 0; i < block_index - sentinel.length() + 1; ++i ) {
            end = block_index - sentinel.length() - i;
            if (memcmp(sentinel.c_str(), &block[end], sentinel.size()) == 0) {
                end += off;
                break;
            }
        }
        assert(block_index >= end);

        if (buf.size() > block_index - end) {
            std::copy(block.begin()+end, block.begin()+block_index, buf.begin());
        } else {
            buf.assign(block.begin()+end, block.begin()+block_index);
        }
        bufsize = block_index - end;
        block_index = end;
    }


    return block_index;
}


size_t Reader::GetBlockByMultipleLines(std::vector<char> &block, size_t n) {  
    size_t block_index = 0;
    if (bufsize > 0) {
        assert(bufsize < block.size());
        std::copy(buf.begin(), buf.begin()+bufsize, block.begin());
        block_index = bufsize;
        bufsize = 0;
    }

    block_index += GetBlock0(&block[0]+block_index, block.size()-block_index);

    if (block_index > 0) {

        size_t end = block_index + 1;
        size_t line_no = 0;
        for (size_t i = 0; i < block_index; ++i) {
            if (block[i] == '\n') {
                line_no ++;
                
                if (line_no % n == 0) {
                    end = i+1;  // skip '\n'
                } 
                printf("%zd line\n", line_no);
            }
            
        }
        if (end > block_index) {
            LOG(ERROR)("Buffer is too small and cannot contain %d line of data. %zd %zd", n, end, block_index);
        }
        
        if (buf.size() > block_index - end) {
            std::copy(block.begin()+end, block.begin()+block_index, buf.begin());
        } else {
            buf.assign(block.begin()+end, block.begin()+block_index);
        }
        bufsize = block_index - end;
        block_index = end;
    }
    
    return block_index;
}

size_t StdioReader::GetLines(std::vector<std::string> &lines) {
    size_t sz = 0;
    for (; sz < lines.size(); ++sz) {
        bool r = (bool)std::getline(std::cin, lines[sz]);
        if (!r) break;
    }
    return sz;
}

size_t StdioReader::GetBlock0(char* block, size_t max_size) {
    return fread(block, 1, max_size, stdin);
}



GzFileReader::GzFileReader(const std::string &fname) {
    in_ = gzopen(fname.c_str(), "rb");
}

std::string GzFileReader::GetLine() {
    std::string str;

    const int BUF_SIZE = 1024*2;
    char buffer[BUF_SIZE];

    char* line = gzgets(in_, buffer, BUF_SIZE);
    while (line != nullptr) {
        str += line;
        if (str.back() == '\n') {
            break;
        }
        line = gzgets(in_, buffer, BUF_SIZE);
    }

    return str;
}

bool GzFileReader::GetLine(std::string &str) {
    str.clear();

    const int BUF_SIZE = 1024*2;
    char buffer[BUF_SIZE];

    char* line = gzgets(in_, buffer, BUF_SIZE);
    while (line != nullptr) {
        str += line;
        if (str.back() == '\n') {
            break;
        }
        line = gzgets(in_, buffer, BUF_SIZE);
    }

    return !str.empty();
}

size_t GzFileReader::GetLines(std::vector<std::string> &lines) {

    const int BUF_SIZE = 1024*2;
    char buffer[BUF_SIZE];

    size_t index = 0; 
    while (index < lines.size()) {
        std::string &str = lines[index];
        str.clear();

        char* line = gzgets(in_, buffer, BUF_SIZE);
        while (line != nullptr) {
            str += line;
            if (str.back() == '\n') {
                break;
            }
            line = gzgets(in_, buffer, BUF_SIZE);
        }

        if (str.empty()) break;

        index ++;
    }

    return index;
}

std::string GzFileReader::GetStrippedLine() {
    return StripString(GetLine());
}

std::string GzFileReader::GetNoEmptyLine() {
    if (!HasNoEmptyLine()) {
        std::string str = StripString(GetLine());
        while (str.empty() && !IsEnd()) {
            str = StripString(GetLine());
        } 
        return str;
    } else {
        next_line_pos_ = -1;
        return next_line_;
    }
}

std::vector<std::string> GetLineFromFile(const std::string& fname) {
    GzFileReader in(fname);
    std::vector<std::string> lines;
    if (in.Valid()) {
        std::string s = in.GetNoEmptyLine();
        while (!s.empty()) {
            lines.push_back(s);
            s = in.GetNoEmptyLine();
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
    return lines;
}

size_t GzFileReader::GetBlock0(char* block, size_t max_size) {
    return gzread(in_, block, max_size);
}


size_t MultiFileReader::GetLines(std::vector<std::string> &lines) {
    auto size = ifile_ == nullptr ? 0 : ifile_->GetLines(lines);
    while (size == 0) {
        if (file_index_ < ifnames_.size()) {
            LOG(INFO)("%s", ifnames_[file_index_].c_str());
            ifile_.reset(new GzFileReader(ifnames_[file_index_++]));
            if (ifile_->Valid()) {
                size = ifile_->GetLines(lines);
            } else {
                LOG(ERROR)("Failed to open file: %s", ifnames_[file_index_-1].c_str());
            }
        } else {
            break;
        }
    }
    return size;
}

size_t MultiFileReader::GetBlock0(char* block, size_t max_size) { 
    auto size = ifile_ == nullptr ? 0 : ifile_->GetBlock0(block, max_size);
    while (size == 0) {
        if (file_index_ < ifnames_.size()) {
            //LOG(INFO)("%s", ifnames_[file_index_].c_str());
            ifile_.reset(new GzFileReader(ifnames_[file_index_++]));
            if (ifile_->Valid()) {
                size = ifile_->GetBlock0(block, max_size);
            } else {
                LOG(ERROR)("Failed to open file: %s", ifnames_[file_index_-1].c_str());
            }
        } else {
            break;
        }
    }
    return size;
}

bool LineInBlock::GetLine(std::string& line) {
    if (block_index_ >= block_.begin() + block_size_) {
        if (mutex_ == nullptr) {
            block_size_ = reader_.GetBlock(block_, "\n", 1);
        } else {
            std::lock_guard<std::mutex> lock(*mutex_);
            block_size_ = reader_.GetBlock(block_, "\n", 1);
        }
        block_index_ = block_.begin();
    }

    if (block_index_ < block_.begin() + block_size_) {
        auto next = std::find(block_index_, block_.begin()+block_size_, '\n');
        line = std::string(block_index_, next);

        block_index_ = next;
        if (block_index_ < block_.begin()+block_size_ && *block_index_ == '\n')  block_index_++;
        return true;
    } else {
        return false;
    }
}

} // namespace fsa {