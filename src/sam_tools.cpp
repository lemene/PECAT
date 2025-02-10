#include "sam_tools.hpp"

#include <htslib/sam.h>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "utils/project_file.hpp"
#include "assemble/read_variants.hpp"

#include "overlap/mapping.hpp"

namespace fsa {



void Program_GroupRead::Running() {


    samFile *in = hts_open(ifname_.c_str(), "r");
    assert(in != nullptr);
    
    bam_hdr_t *header = sam_hdr_read(in);
    assert(header != nullptr);
    auto idx = sam_index_load2(in, ifname_.c_str(), (ifname_+".bai").c_str());
    assert(idx != nullptr);


    std::unordered_set<std::string> done;
    std::vector<std::string> reads;
    size_t accu_len = 0;
    size_t sn = 0;
    for (size_t tgtid = 0; tgtid < header->n_targets; tgtid++) {

        hts_itr_t* itr = sam_itr_queryi(idx, tgtid, 0, header->target_len[tgtid]);
        bam1_t *b = bam_init1();
            
        //while (sam_itr_next(in, itr, b) >= 0) {
        while (sam_read1(in, header, b) >=0 ) {
            //LOG(INFO)("%d: %d", b->core.tid, b->core.pos);
            if (done.find(bam_get_qname(b)) == done.end()) {
                reads.push_back(bam_get_qname(b));
                done.insert(bam_get_qname(b));
                accu_len += b->core.l_qseq;

                if (accu_len >= block_size_) {
                    Save(sn, reads);
                    reads.clear();
                    sn++;
                    accu_len = 0;
                }
            }
            //
        }

        if (accu_len > 0) {
            Save(sn, reads);
            reads.clear();
            sn++;
            accu_len = 0;
        }
        bam_destroy1(b);
    }
    bam_hdr_destroy(header);
    hts_close(in);
    LOG(INFO)("EEEE");
    
}

inline std::string Format(const std::string& pattern, int d) {
    std::string result = pattern;
    std::string s = std::to_string(d);
    result.replace(pattern.find("{}"), 2, s);
    return result;
};
void Program_GroupRead::Save(int id, const std::vector<std::string>& names) {


    std::string ofname = Format(opattern_, id);

    
    std::ofstream of(ofname);
    if (of.is_open()) {
        for (const auto &n : names) {
            of << n << "\n";
        }
    }
    
}

void Program_Test::Running() {
    OverlapStore ol_store;
    ol_store.Load(ifname_);
    
}

} // namespace fsa
