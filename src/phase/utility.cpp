#include "utility.hpp"

#include <sstream>

#include "../utils/logger.hpp"
#include "../utility.hpp"

namespace fsa {
void CoverageOptions::From(const std::string& str) {
    auto items = SplitStringByChar(str, ':');

    for (auto &i : items) {
        auto kv = SplitStringByChar(i, '=');
        if (kv[0] == "lc") {
            range[0] = std::stoi(kv[1]);
        } else if (kv[0] == "hc") {
            range[1] = std::stoi(kv[1]);
        } else if (kv[0] == "lvc") {
            valid_range[0] = std::stoi(kv[1]);
        } else if (kv[0] == "hvc") {
            valid_range[1] = std::stoi(kv[1]);
        } else if (kv[0] == "lr") {
        } else if (kv[0] == "hr") {
            rate[1] = std::stod(kv[1]);
        } else {
            LOG(ERROR)("Unrecoginze coverage option %s", kv[0].c_str());
        }
    }
}

std::string CoverageOptions::ToString() const  {
    std::ostringstream oss;
    oss.precision(2);
    oss.setf(std::ios::fixed);
    oss << "lc=" << range[0]
        << ":hc=" << range[1]
        << ":lr=" << rate[0]
        << ":hr=" << rate[1]
        << ":lvc=" << valid_range[0]
        << ":hvc=" << valid_range[1];;
    return oss.str();
}


void Variant::Comfirm(const CoverageOptions& covopts) {
    
    int cov = counts[9]; // std::accumulate(counts, counts+9, 0);

    if (covopts.Effective(cov)) {
        int threshold = covopts.Threshold(cov);
        std::vector<int> ok;
        for (size_t j=0; j<9; ++j) {
            if (counts[j] >= threshold) {
                ok.push_back(j);
            }
        }
        
        if (ok.size() >=2) {
            std::sort(ok.begin(), ok.end(),[this](int a, int b) { return counts[a] > counts[b]; });
            var[0] = ok[0];
            var[1] = ok[1];
        }
    }

}

} // namespace fsa