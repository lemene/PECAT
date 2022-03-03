#include "read_variants.hpp"

#include "../file_io.hpp"
#include "../utils/logger.hpp"
#include "../utility.hpp"
#include "../overlap.hpp"

namespace fsa {

void ReadVariants::Load(const std::string& fname) {
    // format: ctg read d off [items(ctgoff-b-rdoff-b)]

    const int item_start = 4;
    GzFileReader reader(fname);

    if (reader.Valid()) {
        size_t count = 0;
        std::string line = reader.GetLine();
        while (!line.empty()) {
            auto items = SplitStringBySpace(line);
            assert(items.size() >= item_start);

            if (items.size() > item_start) {
                Variants v;
                v.contig = string_pool_.QueryIdByString(items[0]);
                auto rid = string_pool_.QueryIdByString(items[1]);
                v.d = std::stoi(items[2]);
                v.offset = std::stoi(items[3]);

                for (size_t i=item_start; i < items.size(); ++i) {
                    auto ss = SplitStringByChar(items[i], '|');
                    assert(ss.size() >= 3);
                    v.vars.insert(std::make_pair<int, int>(std::stoi(ss[0]), std::stoi(ss[1])));
                }

                auto r = reads.find(rid);
                if (r != reads.end()) {
                    r->second.push_back(v);
                } else {
                    reads.insert({rid, {v}});
                }
            }
            line = reader.GetLine();
        } 
    } else {
        LOG(WARNING)("Failed to open file: %s", fname.c_str());
    }
}


bool ReadVariants::IsCompatible(const Variants &a, const Variants &b, const Overlap &ol) const {
    if (a.contig != b.contig) return false;
    
    if (a.d == b.d && !ol.SameDirect()) return false;
    if (a.d != b.d && ol.SameDirect()) return false;

    auto atob = b.d == 0 ? (a.offset - b.offset) : -(a.offset - b.offset) ;
    if (std::abs(ol.MappingToTarget<1>({0})[0] - atob) > std::max<int>(1000,atob*0.10)) return false;

    return true;
}

std::array<int, 2> ReadVariants::Test(int a, int b) const {
    auto ra = reads.find(a);
    auto rb = reads.find(b);

    std::array<int, 2> result = {0, 0};
    if (ra != reads.end() && rb != reads.end()) {

        for (auto & ira : ra->second) {
            for (auto & irb : rb->second) {
                if (ira.contig == irb.contig) {
                    std::array<int, 2> rs {0, 0};
                    for (auto & irav : ira.vars) {
                        auto irbv = irb.vars.find(irav.first);
                        if (irbv != irb.vars.end()) {
                            if (irav.second != -1 && irbv->second != -1) {
                                if (irav.second == irbv->second) {
                                    rs[0] ++;
                                } else {
                                    rs[1] ++;
                                }
                            }
                        }
                    }
                    if (result[0] + result[1] < rs[0] + rs[1]) {
                        result = rs;
                    }
                }
            }
        }
    }
    return result;
}

std::array<int, 2> ReadVariants::Test(const Overlap& ol, bool bad) const  {
    auto ra = reads.find(ol.a_.id);
    auto rb = reads.find(ol.b_.id);

    std::array<int, 2> result = {0, 0};
    if (ra != reads.end() && rb != reads.end()) {

        for (auto & ira : ra->second) {
            for (auto & irb : rb->second) {
                if (IsCompatible(ira, irb, ol)) {
                    std::array<int, 2> rs {0, 0};
                    for (auto & irav : ira.vars) {
                        auto irbv = irb.vars.find(irav.first);
                        if (irbv != irb.vars.end()) {
                            if (irav.second != -1 && irbv->second != -1) {
                                if (irav.second == irbv->second) {
                                    rs[0] ++;
                                } else {
                                    rs[1] ++;
                                }
                            }
                        }
                    }
                    if (!bad) {
                        if (result[0] + result[1] < rs[0] + rs[1]) {
                            result = rs;
                        }
                    } else {
                        if ( result[1] < rs[1]) {
                            result = rs;
                        }
                    }
                }
            }
        }
    }
    return result;

}

}