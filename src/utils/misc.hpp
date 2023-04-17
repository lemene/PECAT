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

}