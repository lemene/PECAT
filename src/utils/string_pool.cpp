#include "string_pool.hpp"

#include <cassert>
#include "file_io.hpp"
#include "logger.hpp"
namespace fsa {

StringPool::ID StringPool::GetIdByString(const std::string &str) {
    std::lock_guard<std::mutex> lock(mutex_);
    return GetIdByStringUnsafe(str);
}

StringPool::ID StringPool::GetIdByStringUnsafe(const std::string &str) {
    auto iter = ids_.find(str);
    if (iter != ids_.end()) {
        return iter->second;
    } else {
        StringPool::ID id = (StringPool::ID)strings_.size();
        ids_.insert(std::make_pair(str, id));
        strings_.push_back(str);
        return id;
    }
}

const std::string& StringPool::QueryStringById(ID i) const {
    assert(i >= 0 && i < (int)strings_.size());
    return strings_[i];
}

StringPool::ID  StringPool::QueryIdByString(const std::string &str) const {
    auto iter = ids_.find(str);
    return iter != ids_.end() ? iter->second : -1;
}

void StringPool::Clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    strings_.clear();
    ids_.clear();
}

std::unordered_map<StringPool::ID, StringPool::ID> StringPool::Merge(const StringPool &sp) {
    std::lock_guard<std::mutex> lock(mutex_);

    std::unordered_map<StringPool::ID, StringPool::ID> mapper;
    for (size_t i=0; i<sp.strings_.size(); ++i) {
        mapper[i] = GetIdByStringUnsafe(sp.strings_[i]);
    }

    return mapper;
}

void StringPool::Save(const std::string &fname) const {
    GzFileWriter writer(fname);
    if (writer.Valid()) {
        for (size_t i=0; i<strings_.size(); ++i) {
            writer << i << " " << strings_[i] << "\n";
        }
    } else {
         LOG(ERROR)("Failed to open outfile: %s", fname.c_str());
    }
}


std::unordered_map<Seq::Id, Seq::Id> StringPool::MergeNameId(const TempNameId &ni) {
    std::lock_guard<std::mutex> lock(mutex_);

    std::unordered_map<Seq::Id, Seq::Id> id2id;
    for (auto &i : ni.names_to_ids) {
        id2id[i.second] = GetIdByStringUnsafe(i.first);
    }
    return id2id;
}


}
