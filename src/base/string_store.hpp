#pragma once

#include <cassert>
#include <cstring>

#include <vector>
#include <string>

namespace fsa {

class StringStore {
public:
    class String {
    public:
        String(const char* s, size_t l) : start_(s), len_(l) {}
    protected:
        const char* start_ { nullptr };
        size_t len_ { 0 };
    };
public:
    StringStore(size_t bs=10*1024*1024) : bsize_(bs) { 
        blocks_.push_back(new char[bsize_]); 
        bindex_ = 0;
    }

    ~StringStore() {
        for (auto b : blocks_) { delete[] b; }
    }

    String Create(const std::string& s) {
        return Create(s.c_str(), s.size());
    }
    String Create(const char* s) {
        return Create(s, strlen(s));
    }
    String Create(const char* s, const char* e) {
        return Create(s, (size_t)(e-s));
    }
    String Create(const char* s, size_t len) {
        if (bindex_ + len + 1> bsize_) {
            blocks_.push_back(new char[bsize_]); 
            bindex_ = 0;
        }
        assert(bindex_ + len + 1< bsize_);
        
        memcpy(blocks_.back()+bindex_, s, len + 1);
        bindex_ += len + 1;
        return String(blocks_.back() - len, len);
    }


protected:
    std::vector<char*> blocks_;
    size_t bsize_;
    size_t bindex_ { 0 };
};
}