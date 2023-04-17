#include "project_file.hpp"

#include <regex>

#include "../file_io.hpp"
#include "string_pool.hpp"
#include "./utility.hpp"

namespace fsa {

// std::unordered_map<Seq::Id, std::vector<Seq::Id>> LoadTiles(const std::string& fname, StringPool &sp) {
//     GzFileReader ftile(fname);

//     std::unordered_map<Seq::Id, std::vector<Seq::Id>> tiles;
//     if (ftile.Valid()) {
//         std::vector<Seq::Id> *tls = nullptr;
//         Seq::Id curr_ctgid = Seq::NID;

//         std::string line = ftile.GetNoEmptyLine();
//         while (!line.empty()) {
//             auto its = SplitStringBySpace(line);
//             auto ctgid = sp.GetIdByString(its[0]);
//             bool isnew = false;
//             if (curr_ctgid == Seq::NID || curr_ctgid != ctgid) {
//                 tls = &(tiles[ctgid]);
//                 curr_ctgid = ctgid;
//                 isnew = true;
//             }

//             auto nodes = SplitStringByChar(SplitStringByChar(its[1], '=')[1], '~');
//             if (isnew) {
//                 tls->push_back(sp.GetIdByString(SplitStringByChar(nodes[0], ':')[0]));
//             }
//             tls->push_back(sp.GetIdByString(SplitStringByChar(nodes[1], ':')[0]));
//             line = ftile.GetNoEmptyLine();
//         }


//     } else {
//         LOG(INFO)("Failed to load file: %s", fname.c_str());
//     }

//     return tiles;
// }

void PrjContigTiles::Load(const std::string& fname) {

    GzFileReader reader(fname);

    if (reader.Valid()) {
        std::string line = reader.GetLine();
        int curr = StringPool::NID;
        std::vector<int> *ts;
        while (!line.empty()) {

            std::regex pattern("(ctg\\d+) edge=(\\d+):[BE]\\~(\\d+):[BE]");
            
            std::smatch m;
            bool r = std::regex_search(line, m, pattern);
            assert(r);
            auto ctg = string_pool_.GetIdByString(m.str(1));
            auto r0 = string_pool_.GetIdByString(m.str(2));
            auto r1 = string_pool_.GetIdByString(m.str(3));
            
            if (ctg != curr) {
                ts = &tiles_[ctg];
                curr = ctg; 
                ts->push_back(r0);
                ts->push_back(r1);
            } else {
                ts->push_back(r1);
            }

            line = reader.GetLine();
        } 
    } else {
        LOG(WARNING)("Failed to open file: %s", fname.c_str());
    }
}


void PrjVariants::Load(const std::string& fname) {

    GzFileReader reader(fname);

    if (reader.Valid()) {
        std::string line = reader.GetLine();
        while (!line.empty()) {
            auto its = SplitStringBySpace(line);
            if (its.size() == 13 && its.back() == "1") {
                auto id = string_pool_.GetIdByString(its[0]);
                size_t pos = std::stoul(its[1]);
                std::array<size_t, 4> count = { 
                    std::stoul(its[2]), 
                    std::stoul(its[3]), 
                    std::stoul(its[4]), 
                    std::stoul(its[5])};
                std::array<uint8_t,2> snp = {0, 1};
                if (count[1] > count[0]) snp = {1, 0};
                for (size_t i = 2; i < 4; ++i) {
                    if (count[i] > count[snp[0]]) {
                        snp[1] = snp[0];
                        snp[0] = i;
                    } else if (count[i] > count[snp[1]]) {
                        snp[1] = i;
                    }
                }

                variants_[id][pos] = snp;
            }

            line = reader.GetLine();
        } 
    } else {
        LOG(WARNING)("Failed to open file: %s", fname.c_str());
    }
}
}