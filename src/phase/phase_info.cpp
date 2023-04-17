#include "phase_info.hpp"

#include <cassert>
#include <array>

#include "../file_io.hpp"
#include "../utils/logger.hpp"
#include "../utility.hpp"

namespace fsa {

void PhaseInfoFile::Load(const std::string &fname) {
    
    GzFileReader in(fname);
    if (in.Valid()) {
        std::string line;
        while (in.GetLine(line)) {
            auto sets = SplitStringByChar(line, ':');
            assert(sets.size() == 2);
            auto set1 = SplitStringByChar(sets[1], ',');
            auto tid = string_pool_.GetIdByString(sets[0]);

            auto &www = phased_[tid];

            for (auto &i : set1) {
                auto item = SplitStringByChar(i, '|');
                assert(item.size() == 3 || item.size() == 1);
                if (item.size() == 1) break; //TODO 
                if (item[0] != "-1") {
                    auto qid = string_pool_.GetIdByString(item[0]);
                    www[qid].push_back({std::stoi(item[1]), std::stoi(item[2])});
                } else {
                    removed_.insert(tid);
                }
            }
        }
    }

    LOG(INFO)("Load pairs: %zd", phased_.size());
}

void PhaseInfoFile::Load(const std::string &fname, size_t thread_size) {
    std::mutex mutex_combine;
    std::mutex mutex_generate;

    GzFileReader in(fname);
    std::array<long long, 2> done = {0,0};

    int block_size = 2000;
    auto generate = [&mutex_generate, &in, &done](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_generate);
        auto size = in.GetLines(lines);
        done[0] += size;

        return size;
    };

    auto combine = [this, &mutex_combine](std::unordered_map<StringPool::ID, StringPool::ID> &mapper, 
                                          std::unordered_map<int, std::unordered_map<int, std::vector<Item>>> &phased) {
        std::lock_guard<std::mutex> lock(mutex_combine);

        for (auto i : phased) {
            auto tid = mapper[i.first];
            auto &www = phased_[tid];
            for (auto j : i.second) {
                auto qid = mapper[j.first];
                www[qid].insert(www[qid].end(), j.second.begin(), j.second.end());
            }
        }
    };

    auto work = [&](size_t threadid) {
        StringPool subpool;
        std::vector<std::string> lines(block_size);

        std::unordered_map<int, std::unordered_map<int, std::vector<Item>>> phased;
        
        size_t size = generate(lines);
        while (size > 0) {
            for (size_t i=0; i<size; ++i) {
                const std::string& line = lines[i];
                auto sets = SplitStringByChar(line, ':');
                assert(sets.size() == 2);
                auto set1 = SplitStringByChar(sets[1], ',');
                auto tid = subpool.GetIdByString(sets[0]);

                auto &www = phased[tid];

                for (auto &i : set1) {
                    auto item = SplitStringByChar(i, '|');
                    assert(item.size() == 3 || item.size() == 1);
                    if (item.size() == 1) break; //TODO 
                    auto qid = subpool.GetIdByString(item[0]);
                    www[qid].push_back({std::stoi(item[1]), std::stoi(item[2])});
                }
            }
            

            auto mapper = string_pool_.Merge(subpool);
            subpool.Clear();

            combine(mapper, phased);
            phased.clear();
            size = generate(lines);
        }
    };


    if (in.Valid()) {
        MultiThreadRun(thread_size, work);
    } else {
        LOG(ERROR)("Failed to load file: %s", fname.c_str());
    }
    LOG(INFO)("Load ignored pairs: %zd", phased_.size());

}

bool PhaseInfoFile::Contain(const Overlap &o, int threshold, bool both) const {

    int off_1to0 = o.MappingToSource<1>({0})[0];
    int off_0to1 = o.MappingToTarget<1>({0})[0];
    return Contain(o.a_.id, o.b_.id, o.SameDirect(), off_0to1, off_1to0, threshold, both);
}


size_t PhaseInfoFile::MinDistance(const Overlap &o, bool both) const {
    size_t min_distance0 = MAX_DISTANCE;
    auto ia = phased_.find(o.a_.id);
    if (ia != phased_.end()) {
        auto iab = ia->second.find(o.b_.id);
        if (iab != ia->second.end()) {
            int offset = o.MappingToSource<1>({0})[0];
            for (auto &l : iab->second) {
                if (((o.SameDirect() && l.strand == 0) || (!o.SameDirect() && l.strand == 1) ) && std::abs(offset - l.offset) < min_distance0) {
                    min_distance0 = std::abs(offset - l.offset);
                }
            }
        }
    }

    size_t min_distance1 = MAX_DISTANCE;
    auto ib = phased_.find(o.b_.id);
    if (ib != phased_.end()) {
        auto iba = ib->second.find(o.a_.id);
        if (iba != ib->second.end()) {

            int offset = o.MappingToTarget<1>({0})[0];
            for (auto &l : iba->second) {
                if (((o.SameDirect() && l.strand == 0) || (!o.SameDirect() && l.strand == 1) ) && std::abs(offset - l.offset) < min_distance1) {
                    min_distance1 = std::abs(offset - l.offset);
                }
            }
            
        }
        
    }
    return both ? std::max(min_distance0, min_distance1) : std::min(min_distance0, min_distance1);
}

bool PhaseInfoFile::Contain(const std::string &name0, const std::string &name1) {
    int id0 = string_pool_.QueryIdByString(name0);
    int id1 = string_pool_.QueryIdByString(name1);
    return Contain(id0, id1);
}

bool PhaseInfoFile::Contain(int id0, int id1) const {
    auto ia = phased_.find(id0);
    if (ia != phased_.end()) {
        auto iab = ia->second.find(id1);
        if (iab != ia->second.end()) {
            return true;
        }
    }

    auto ib = phased_.find(id1);
    if (ib != phased_.end()) {
        auto iba = ib->second.find(id0);
        if (iba != ia->second.end()) {
            return true;
        }
    }

    return false;
}

bool PhaseInfoFile::Contain(int id0, int id1, bool d, int off_0to1, int off_1to0, int threshold, bool both) const {
    std::array<bool,2> test = {false, false};
    auto ia = phased_.find(id0);
    if (ia != phased_.end()) {
        auto iab = ia->second.find(id1);
        if (iab != ia->second.end()) {
            if (threshold >= 0) {
                for (auto &l : iab->second) {
                    if (((l.strand == 0 && d) || (l.strand == 1 && !d)) && std::abs(off_1to0 - l.offset) < threshold) {
                        test[0] = true;
                        break;
                    }
                }
            } else {
                test[0] = true;
            }
        }
    }

    auto ib = phased_.find(id1);
    if (ib != phased_.end()) {
        auto iba = ib->second.find(id0);
        if (iba != ia->second.end()) {
            if (threshold >= 0) {
                for (auto &l : iba->second) {
                    if (((l.strand == 0 && d) || (l.strand == 1 && !d)) && std::abs(off_0to1 - l.offset) < threshold) {
                        test[1] = true;
                        break;
                    }
                }
            } else {
                test[1] = true;
            }
        }
    }

    return both ? (test[0] && test[1]) : (test[0] || test[1]);
}

std::unordered_set<int> PhaseInfoFile::Get(int id) const {
    std::unordered_set<int> result;
    auto it = phased_.find(id);
    if (it != phased_.end()) {
        for (auto i : it->second) {
            result.insert(i.first);
        }
    }

    return result;

}


}