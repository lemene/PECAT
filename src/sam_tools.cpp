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

void Program_N50::Running() {
    std::vector<int> ctg_lens;
    
    samFile *in = hts_open(ifname_.c_str(), "r");
    assert(in != nullptr);
    
    bam_hdr_t *header = sam_hdr_read(in);
    assert(header != nullptr);

    for (size_t tgtid = 0; tgtid < header->n_targets; tgtid++) {
        ctg_lens.push_back(header->target_len[tgtid]);
    }
    bam_hdr_destroy(header);


    std::sort(ctg_lens.begin(), ctg_lens.end(), [](int a, int b) { return a > b; });

    long long total_length = std::accumulate(ctg_lens.begin(), ctg_lens.end(), (long long)0);
    long long genome_size = genome_size_ == 0 ? total_length : genome_size_;

    std::cout << "Genome: " << genome_size << "\n";
    std::cout << "Count: " << ctg_lens.size() << "\n";
    std::cout << "Total: " << total_length << "\n";
    std::cout << "Max: " << ctg_lens.front() << "\n";
    std::cout << "Min: " << ctg_lens.back() << "\n";
    
    long long accu = 0;
    int ns[] = { 25, 50, 75};
    size_t ins = 0;
    for (size_t i=0; i<ctg_lens.size(); ++i) {
        accu += ctg_lens[i];
        for (; ins < sizeof(ns) / sizeof(ns[0]) && accu > (long long)ns[ins]*1.0/100 * genome_size; ++ins) {
            std::cout << "N" << ns[ins] << ": " << ctg_lens[i] << "\n";
            std::cout << "L" << ns[ins] << ": " << i+1 << "\n";
        }
        if (ins >= sizeof(ns) / sizeof(ns[0])) {
            break;
        }
    }

}

void Program_Stat::Running() {
    ReadStore ref;
    ref.Load(ref_fname_);
    InitializeTables();

    
    samFile *in = hts_open(ifname_.c_str(), "r");
    assert(in != nullptr);
    
    bam_hdr_t *header = sam_hdr_read(in);
    assert(header != nullptr);
    auto idx = sam_index_load2(in, ifname_.c_str(), (ifname_+".bai").c_str());
    assert(idx != nullptr);

    size_t accu_len = 0;
    size_t sn = 0;
    for (size_t it = 0; it < header->n_targets; it++) {
        auto rseq = ref.GetSeq(header->target_name[it]);
        hts_itr_t* itr = sam_itr_queryi(idx, it, 0, header->target_len[it]);
        bam1_t *b = bam_init1();
        auto& table = tables_[header->target_name[it]];
            
    char nt16_2_nt4[16] = {-1};
    nt16_2_nt4[1] = 0;
    nt16_2_nt4[2] = 1;
    nt16_2_nt4[4] = 2;
    nt16_2_nt4[8] = 3;
        while (sam_itr_next(in, itr, b) >= 0) {
        //while (sam_read1(in, header, b) >=0 ) {
            //LOG(INFO)("%d: %d", b->core.tid, b->core.pos);

            //if ((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FSUPPLEMENTARY)) continue; 
            if ((b->core.flag & BAM_FSECONDARY)) continue; 
            
            auto r_start = b->core.pos;
            uint8_t *qseq = bam_get_seq(b);
            
            uint32_t *cigar = bam_get_cigar(b);
            int rpos = 0;
            int qpos = 0;
            int match = 0;
            int clip = 0;
            for(int ic = 0; ic < b->core.n_cigar; ++ic){
                int icigar = cigar[ic];
                int n = bam_cigar_oplen(icigar);
                char t = bam_cigar_opchr(icigar);
                switch (t) {
                    
                case '=':
                    for (size_t i = 0; i < n; ++i) {
                        table[r_start + rpos + i].match ++;
                    }
                    match += n;
                    rpos += n;
                    qpos += n;
                    break;
                case 'M':
                    for (size_t i = 0; i < n; ++i) {
                        if (rseq[r_start + rpos + i] == nt16_2_nt4[bam_seqi(qseq,qpos+i)]) {
                            table[r_start + rpos + i].match ++;
                            if (r_start + rpos + i == 554322) {
                                printf("m= %zd, %zd\n", rseq[r_start + rpos + i], nt16_2_nt4[bam_seqi(qseq,qpos+i)]);
                                printf("m= %s\n", bam_get_qname(b));

                            } 
                        } else {
                            table[r_start + rpos + i].mismatch[nt16_2_nt4[bam_seqi(qseq,qpos+i)]] ++;
                            if (r_start + rpos + i == 554322) {
                                printf("mx %zd, %zd, %d\n", rseq[r_start + rpos + i], nt16_2_nt4[bam_seqi(qseq,qpos+i)],bam_seqi(qseq,qpos+i));
                                printf("mx %s\n", bam_get_qname(b));
                            } 
                        }
                    }
                    match += n;
                    rpos += n;
                    qpos += n;
                    break;
                case 'X':
                    for (size_t i = 0; i < n; ++i) {
                        table[r_start + rpos + i].mismatch[qseq[qpos+i]] ++;
                    }
                    rpos += n;
                    qpos += n;
                    break;
                case 'I':
                    qpos += n;
                    table[r_start + rpos].insertion ++;
                    table[r_start + rpos].inssize += n;
                    break;
                case 'D':
                    for (size_t i = 0; i < n; ++i) {
                        table[r_start + rpos + i].deletion ++;
                    }
                    rpos += n;
                    break;
                case 'S':
                    qpos += n;
                    clip += n;
                    break;
                case 'H':
                    clip += n;
                    break;
                default:
                    //LOG(ERROR)("Not support %c", t);
                    break;
                }
            }
            //ol.identity_ = match * 2.0 / (qpos + rpos);

            //ol.a_.len = clip+qpos;

            // assert (b->core.n_cigar >= 1);
            // char cigar0 = bam_cigar_opchr(cigar[0]);
            // char cigar1 = bam_cigar_opchr(cigar[b->core.n_cigar-1]);

            // if ((cigar0 == 'S' || cigar0 == 'H') && ol.a_.strand == 0) {
            //     ol.a_.start = bam_cigar_oplen(cigar[0]);
            //     ol.a_.end = ol.a_.start + qpos;
            // } else if ((cigar1 == 'S' || cigar1 == 'H')  && ol.a_.strand == 1) {
            //     ol.a_.start = bam_cigar_oplen(cigar[b->core.n_cigar-1]);
            //     ol.a_.end = ol.a_.start + qpos;
            // } else {
            //     ol.a_.start = 0;
            //     ol.a_.end = ol.a_.start + qpos;
            // }
            
        }

        bam_destroy1(b);
    }
    PrintAccuracy();
    bam_hdr_destroy(header);
    hts_close(in);
    LOG(INFO)("EEEE");
}


void Program_Stat::InitializeTables() {
    
    samFile *in = hts_open(ifname_.c_str(), "r");
    assert(in != nullptr);
    
    bam_hdr_t *header = sam_hdr_read(in);
    assert(header != nullptr);

    for (size_t i = 0; i < header->n_targets; i++) {
        tables_[header->target_name[i]] = std::vector<Position>(header->target_len[i]);
    }
    bam_hdr_destroy(header);
}

void Program_Stat::PrintAccuracy() {
    size_t ref_ok = 0;
    size_t ref_err = 0;
    size_t ref_total = 0;
    size_t ref_ins = 0;
    
    size_t rd_ok = 0;
    size_t rd_err = 0;
    size_t rd_total = 0;
    size_t rd_ins = 0;
    
    size_t i = 0;
    for (auto& it : tables_) {
        for (auto &t : it.second) {
            // if ( i >= 0) {
            //     printf("T: %zd %zd, (%zd, %zd, %zd, %zd), %zd, %zd\n", i, t.match, t.mismatch[0], t.mismatch[1], t.mismatch[2], t.mismatch[3], 
            //         t.deletion, t.insertion);
            //     //break;
            // }
            assert(it.second[i].match == t.match);
            
            std::array<size_t, 7> count = {t.match, 
                t.mismatch[0], t.mismatch[1], t.mismatch[2], t.mismatch[3], 
                t.deletion, t.insertion
            };

            auto mx = std::max_element(count.begin(), count.end()) - count.begin();
            if (mx == 0) {
                ref_ok ++;
                rd_err += std::accumulate(count.begin()+1, count.begin()+5, 0);
                rd_ok += count[0];
            } else if (mx >= 1 && mx <= 4) {   
                printf("ref_err: %zd, %zd, %zd\n", ref_err, i)     ;        
                printf("%zd, (%zd, %zd, %zd, %zd), %zd, %zd\n", t.match, t.mismatch[0], t.mismatch[1], t.mismatch[2], t.mismatch[3], 
                     t.deletion, t.insertion); 
                //assert(0);
                ref_err ++;
                rd_ok += count[mx];
                rd_err = std::accumulate(count.begin(), count.begin()+5, 0) - count[mx];
            } else if (mx == 5) {
                rd_ins += std::accumulate(count.begin(), count.begin()+5, 0);
            } else {
                assert(mx == 6);
                ref_ins += t.inssize / t.insertion;
            }
            rd_total += std::accumulate(count.begin(), count.begin()+5, 0) + t.inssize;
            rd_ins += t.insertion;
            
            i++;
        }
        ref_total += it.second.size();

    }
    printf("%.02f\n", ref_err*1.0/ref_total);
    printf("%zd, %zd, %zd, %zd, %zd\n",ref_total, ref_err, ref_ok, rd_ok, rd_err);
}

void Program_Test::Running() {
    OverlapStore ol_store;
    ol_store.Load(ifname_);
    
}

} // namespace fsa
