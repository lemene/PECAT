#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <mutex>

#include "logger.hpp"
#include "../sequence.hpp"
namespace fsa {


class StringPool {
public:
    typedef int ID;
    const static ID NID = -1;

    struct NameId {
        virtual Seq::Id GetIdByName(const std::string &name) = 0;
        virtual std::string QueryNameById(Seq::Id id) const = 0;
    };

    struct TempNameId : public NameId {
        Seq::Id GetIdByName(const std::string &name) {
            auto it = names_to_ids.find(name);
            if (it != names_to_ids.end()) {
                return it->second;
            } else {
                size_t i = names_to_ids.size();
                names_to_ids[name] = i;
                names.push_back(name);
                return i;
            }
        }
        std::string QueryNameById(Seq::Id id) const {
            return id < (int)names.size() ? names[id] : "";
        }

        std::unordered_map<std::string, Seq::Id> names_to_ids;
        std::vector<std::string> names;
    };

    struct TempNameId2 : public NameId {
        TempNameId2(const StringPool& sp) : string_pool_(sp) {}
        Seq::Id GetIdByName(const std::string &name) {
            Seq::Id id = string_pool_.QueryIdByString(name);
            if (id == StringPool::NID) {
                id = name_id_.GetIdByName(name);
                id += string_pool_.Size();
            }
            return id;
        }
        std::string QueryNameById(Seq::Id id) const {
            assert(id >= 0);
            if (id < (int)string_pool_.Size()) {
                return string_pool_.QueryStringById(id);
            } else {
                return name_id_.QueryNameById(id - string_pool_.Size());
            }
        }
        TempNameId name_id_;
        const StringPool& string_pool_;
    };

    struct SafeNameId : public NameId {
        SafeNameId(StringPool& rs) : rd_store(rs) {}
        virtual Seq::Id GetIdByName(const std::string &name) { 
            return rd_store.GetIdByString(name);
        }

        std::string QueryNameById(Seq::Id id) const  { return rd_store.QueryStringById(id); }

        StringPool& rd_store;

    };

    struct UnsafeNameId : public NameId {
        UnsafeNameId(const StringPool& rs) : rd_store(rs) {}
        virtual Seq::Id GetIdByName(const std::string &name) { 
            Seq::Id id = rd_store.QueryIdByString(name);
            return id;
        }

        std::string QueryNameById(Seq::Id id) const { return rd_store.QueryStringById(id); }
        const StringPool& rd_store;

    };

    ID GetIdByString(const std::string &str);
    ID GetIdByStringUnsafe(const std::string &str);
    const std::string& QueryStringById(int ID) const;
    ID QueryIdByString(const std::string &str) const;
    size_t Size() const { return strings_.size(); }
    bool Valid(ID id) { return id >= 0 && id < (int)strings_.size();}
    void Clear();
    void Save(const std::string &fname) const;
    std::unordered_map<StringPool::ID, StringPool::ID> Merge(const StringPool &sp);
    std::unordered_map<Seq::Id, Seq::Id> MergeNameId(const TempNameId &ni);
protected:
    std::mutex mutex_;
    std::vector<std::string> strings_;
    std::unordered_map<std::string, ID> ids_;
    std::string empty_string_;
};

}