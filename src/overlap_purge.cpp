#include "overlap_purge.hpp"
#include "./utils/logger.hpp"

#include <numeric> 
#include <algorithm>
#include <list>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>

namespace fsa {

ArgumentParser OverlapPurge::GetArgumentParser() {
    ArgumentParser ap("fsa_ol_purge", "tools about reads", "1.0");

    ap.AddPositionOption(ifname_, "ifname", "overlap file name");
    ap.AddPositionOption(ofname_, "ofname", "overlap file name");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    ap.AddNamedOption(read_file_, "read_file", "read file");
    return ap;
}


void OverlapPurge::Running() {
    
    LoadFiles();

    LOG(INFO)("Group overlaps");
    ol_store_.Group(groups_, thread_size_);

    std::array<size_t, 2> range = rd_store_.GetIdRange();
    for (size_t i=range[0]; i<range[1]; ++i) {
        read_ids_.push_back(i);
    }

    GroupReadIds();
    printf("tick: %zd\n", group_ticks.size());
    SaveName(ofname_, 0);
}
 
void OverlapPurge::LoadFiles() {
   if (!read_file_.empty()) {
        rd_store_.Load(read_file_, "", false);
    }

    if (rd_store_.GetIdRange()[1] > 0) {    // read names have been loaded
        ol_store_.LoadFast(ifname_, "", (size_t)thread_size_);
    } else {
        ol_store_.Load(ifname_, "", (size_t)thread_size_);
    }
}

void OverlapPurge::GroupReadIds() {
    std::unordered_map<Seq::Id, bool> done;
    std::unordered_map<Seq::Id, int> lens;

    for (auto i : read_ids_) {
        done[i] = false;
        auto iter = groups_.find(i);
        if (iter == groups_.end()) {
            lens[i] = 0;
        } else {
            lens[i] = iter->second.begin()->second->GetRead(i).len;
        }
    }

    std::sort(read_ids_.begin(), read_ids_.end(), [&lens](int a, int b) {return lens[a] > lens[b]; });

    group_ticks.push_back(0);

    for (auto i : read_ids_) {
        if (!done[i]) {
            size_t s = grouped_ids_.size();
            grouped_ids_.push_back(i);
            done[i] = true;

            while (s < grouped_ids_.size()) {
                auto &ols = groups_[grouped_ids_[s]];
                for (auto &o : ols) {
                    auto d = done.find(o.first);
                    if (d != done.end() && !d->second) {
                        grouped_ids_.push_back(o.first);
                        d->second = true;
                        if (grouped_ids_.size() - group_ticks.back() >= group_size[1]) {
                            break;
                        }
                    }
                }
                if (grouped_ids_.size() - group_ticks.back() >= group_size[0]) {
                    break;
                } else {
                    s++;
                }

            }
            if (grouped_ids_.size() - group_ticks.back() >= group_size[0]) {
                group_ticks.push_back(grouped_ids_.size());
            }
        }
    }

    if (grouped_ids_.size() > group_ticks.back()) {
        group_ticks.push_back(grouped_ids_.size());
    }


    // check whether grouping result is correct
    // assert(grouped_ids_.size() == read_ids_.size());
    //for (auto i: grouped_ids_) {
    //    assert(done.find(i) != done.end());
    //}

}

void OverlapPurge::SaveName(const std::string &opattern, size_t group_size) {

    auto format = [](const std::string& pattern, int d) {
        std::string result = pattern;
        std::string s = std::to_string(d);
        result.replace(pattern.find("{}"), 2, s);
        return result;
    };

    
    for (size_t i=1; i<group_ticks.size(); ++i) {
        std::ofstream of(format(opattern, i-1));
        for (size_t j = group_ticks[i-1]; j < group_ticks[i]; ++j) {
            
            of << rd_store_.QueryNameById(grouped_ids_[j]) << "\n";
        }
    }

}

} // namespace fsa {
    