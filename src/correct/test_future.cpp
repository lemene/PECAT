#include "test_future.hpp"

#include "../file_io.hpp"
#include "../utility.hpp"

namespace fsa {

void SnpFile::Load(const std::string &fname) {
    GzFileReader ofile(fname);

    if (ofile.Valid()) {
        while (!ofile.IsEnd()) {
            std::string line = ofile.GetNoEmptyLine();
            auto its = SplitStringBySpace(line);
            if (its.size() > 2) {
                auto id = sp_.QueryIdByString(its[0]);
                if (id != sp_.NID) {
                    auto &snpset = GetSnpSet(id);
                    for (size_t i = 1; i < its.size(); ++i) {
                        snpset.insert(std::stoi(its[i]));
                    }
                }
            }
        }
    }

}

const std::unordered_set<int>& SnpFile::QuerySnps(StringPool::ID id) const {
    auto iter = snps_.find(id);
    if (iter != snps_.end()) {
        return iter->second;
    } else {
        return empty_;
    }
}

std::unordered_set<int>& SnpFile::GetSnpSet(StringPool::ID id) {
    return snps_[id];
}


} // namespace fsa {
