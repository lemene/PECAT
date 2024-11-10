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


void SnpStore::Load(const std::string& fname) {

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

void SnpStore::LoadFromVcf(const std::string &fname) {
    GzFileReader vcf(fname);

    size_t count = 0;
    auto line = vcf.GetNoEmptyLine();
    for (auto line = vcf.GetNoEmptyLine(); !line.empty(); line = vcf.GetNoEmptyLine()) {
        if (line[0] == '#') continue;
        const int CHROM = 0;
        const int POS = 1;
        const int REF = 3;
        const int ALT = 4;
        const int FILTER = 6;
        const int SAMPLE = 9;
        auto its = SplitStringBySpace(line);
        assert(its.size() >= 10);
        std::array<uint8_t,2> bases = {(uint8_t)-1, (uint8_t)-1};
        if (its[FILTER] == "PASS" && its[SAMPLE].size() >= 3) {

            if (its[SAMPLE].find("0/1") == 0 ) {
                if (its[REF].size() == 1 && its[ALT].size() == 1) {
                    bases[0] = DnaSeq::Serial(its[REF][0]);
                    bases[1] = DnaSeq::Serial(its[ALT][0]);
                    count ++;
                }

            } else if (its[SAMPLE].find("1/2") == 0) {
                if (its[ALT].size() == 3) { //A,C
                    bases[0] = DnaSeq::Serial(its[ALT][0]);
                    bases[1] = DnaSeq::Serial(its[ALT][2]);
                    count ++;
                }
            } else {
                // pass
            }
        }

        if (bases[0] != (uint8_t)-1) {
            auto id = string_pool_.QueryIdByString(its[CHROM]);
            auto pos = std::stoi(its[POS]) - 1;
            if (bases[0] > bases[1]) std::swap(bases[0], bases[1]);

            variants_[id][pos] = bases;
        }

        
    }
    LOG(INFO)("Load SNPs: count = %zd", count);
}

}