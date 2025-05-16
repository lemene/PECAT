#include "read_call.hpp"

#include <unordered_set>
#include <iostream>
#include <atomic>
#include "../utils/logger.hpp"
#include "overlap_store.hpp"
#include "utility.hpp"
#include "contig_phaser.hpp"
//#include "local_phaser.hpp"

namespace fsa {
    
bool ReadCall::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}


void ReadCall::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser ReadCall::GetArgumentParser() {
    ArgumentParser ap;
    ap.AddNamedOption(vcf_fname_, "vcf_fname", "VCF filename that contains snp information");
    return ap;
}


void ReadCall::CheckArguments() {
}


void ReadCall::Running() {
    snp_store_.LoadFromVcf(vcf_fname_);

    // CallVariantInReads

    // Save

}


void ReadCall::CallVariantsInReads(std::unordered_map<ReadOffset, ReadInfo>& read_infos, const std::vector<Variant> &vars) {
    const int C = 3;
    const int variant_type = 1;
    
    for (auto &iter : read_infos) {
        const Overlap& o = *(iter.second.o);

        std::unordered_map<int, std::array<int,5>> cand_vars;

        const auto &rd = rd_store_.GetSeq(o.a_.id); 
        size_t ctg_off = 0;
        size_t rd_off = 0;
        for (const auto &d : o.detail_) {

            switch (d.type) {
            case 'M':
                if (variant_type & VARIANT_TYPE_M) {
                    if (d.len >= C) {
                        for (int i=C/2; i<d.len-C/2; ++i) {
                            size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+i : o.b_.end-ctg_off-i-1;
                            size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+i : o.a_.end-rd_off-i-1;

                            uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : 3 - rd[rd_i];

                            if (vars[ctg_i].Valid()) {
                                if (vars[ctg_i].AtM(rd_b)) {
                                    cand_vars[ctg_i] = {(int)ctg_i, (int)rd_i, rd_b, vars[ctg_i].Offset(rd_b), -1 };
                                    //iter.second.vars.push_back({(int)ctg_i, (int)rd_i, rd_b, vars[ctg_i].Offset(rd_b), -1 });
                                } 
                            }
                        }
                    }
                }
                ctg_off += d.len;
                rd_off += d.len;
                break;

            case 'D':
                if (variant_type & VARIANT_TYPE_D) {
                    for (int i=0; i<d.len; ++i) {
                        size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+i : o.b_.end-ctg_off-i-1;
                        size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+i : o.a_.end-rd_off-i-1;
                        if (vars[ctg_i].AtD()) {
                            //iter.second.vars.push_back({(int)ctg_i, (int)rd_i, (uint8_t)8});
                            cand_vars[ctg_i] = {(int)ctg_i, (int)rd_i, (uint8_t)8, vars[ctg_i].Offset(8), -1};
                        }
                    }
                }
                ctg_off += d.len;
                break;

            case 'I':
                if (variant_type & VARIANT_TYPE_I) {
                    size_t ctg_i = o.b_.strand == 0 ? o.b_.start+ctg_off+0 : o.b_.end-ctg_off-0-1;
                    size_t rd_i = o.a_.strand == 0 ? o.a_.start+rd_off+0 : o.a_.end-rd_off-0-1;
                    uint8_t rd_b = o.a_.strand == o.b_.strand ? rd[rd_i] : 3 - rd[rd_i];
                    if (vars[ctg_i].AtI(rd_b)) {
                        //iter.second.vars.push_back({(int)ctg_i, (int)rd_i, rd_b+4});
                        cand_vars[ctg_i] = {(int)ctg_i, (int)rd_i, rd_b+4, vars[ctg_i].Offset(rd_b+4), -1};
                    }
                }
                rd_off += d.len;
                break;

            case 'S':
            case 'H':
                break;
            case '=':
            default:
                LOG(ERROR)("Not support cigar type '%c'.", d.type);
            }
            
        }

        for (size_t i = o.b_.start; i < (size_t)o.b_.end; ++i) {
            if (vars[i].Valid()) {
                auto it = cand_vars.find(i);
                if (it != cand_vars.end()) {
                    iter.second.vars.push_back(it->second);
                } else {
                    iter.second.vars.push_back({int(i), -1, -1, -1, -1});
                }
            }
        }
 
    }


}

} // namespace fsa {
