#include "misc.hpp"

#include <regex>

#include "../file_io.hpp"
#include "string_pool.hpp"


namespace fsa {
void PrjFile::ContigTiles::Load(const std::string& fname) {

    const int item_start = 4;
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

} // namespace fsa {