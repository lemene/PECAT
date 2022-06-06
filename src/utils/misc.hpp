#pragma once

#include <vector>
#include <sstream>
#include <unordered_map>

namespace fsa {

// 解析形如 "diff:s=100:e=0.9" 的参数
class StringOptions {
public:
    struct Output {
        Output() {
            oss.precision(2);
            oss.setf(std::ios::fixed);
        }

        template<typename T>
        void Add(const std::string& n, const T& v) {
            if (count > 0) {
                oss << ":";
            }
            oss << n << "=" << v;
            count++;
        }
        std::string ToString() const { return oss.str(); }
        std::ostringstream oss;
        int count { 0 };
    };

};


struct PrjFile {
    class ContigTiles {
    public:
        ContigTiles(class StringPool &sp) : string_pool_(sp) {}

        void Load(const std::string &fname);
        const std::vector<int>& GetTiles(int& ctg) {
            auto iter = tiles_.find(ctg);
            return iter != tiles_.end() ? iter->second : empty_;
        }

        const std::unordered_map<int, std::vector<int>> & GetTiles() const { return tiles_; }
    protected:
        std::unordered_map<int, std::vector<int>> tiles_;
        std::vector<int> empty_;
        class StringPool &string_pool_;
    };
};

}