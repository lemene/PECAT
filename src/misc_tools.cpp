#include "misc_tools.hpp"

#include <sstream>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "sequence_store.hpp"
#include "file_io.hpp"
#include "utils/project_file.hpp"

#include "phase/hic_read_infos.hpp"

namespace fsa {


void Program_SplitOverlaps::Running() {
    // 只支持 paf
    
    std::string itype = OverlapStore::DetectFileType(ifname_);
    assert(itype == "txt" || itype == "txt.gz" || itype == "paf" || itype == "paf.gz" );

    std::vector<std::string> ifnames = itype == "txt" ? GetLineFromFile(ifname_) : std::vector<std::string>({ifname_});
    std::shared_ptr<GzFileReader> ifile;
    size_t ifile_index = 0;    
 

    LOG(INFO)("Start loading reads");
    ReadStore rd_store(string_pool_);
    rd_store.Load(rd_fname_, "", false);

    LOG(INFO)("Start loading overlaps");
    OverlapStore ol_store(rd_store.GetStringPool());
    ol_store.LoadFast(ifname_, "", thread_size_);
    LOG(INFO)("End loading overlaps");
    
    return;
    std::mutex mutex_combine;
    std::mutex mutex_generate;

    auto generate = [&ifile, &ifnames, &ifile_index,&mutex_generate](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_generate);
        auto size = ifile == nullptr ? 0 : ifile->GetLines(lines);
        while (size == 0) {
            if (ifile_index < ifnames.size()) {
                LOG(INFO)("%s", ifnames[ifile_index].c_str());
                ifile.reset(new GzFileReader(ifnames[ifile_index++]));
                if (ifile->Valid()) {
                    size = ifile->GetLines(lines);
                } else {
                    LOG(ERROR)("Failed to open file: %s", ifnames[ifile_index-1].c_str());
                }
            } else {
                break;
            }
        }
        return size;
    };

    auto work = [&](size_t threadid) {

        int line_size = 10000;
        std::vector<std::string> lines(line_size);
        
        size_t size = generate(lines);
        int count = 0;
        while (size > 0) {
            count += size;
            for (size_t i=0; i<size; ++i) {
                  auto its = SplitStringBySpace(lines[i]);
            //    infos.CheckFull(block_size_);
            //     char aname[1024];
            //     char bname[1024];
            //     int alen;
            //     int blen;
            //     get_name_len_from_paf_line(lines[i], aname, alen, bname, blen);
            //     auto aid = string_pool_.GetIdByStringUnsafe(aname);
            //     auto bid = string_pool_.GetIdByStringUnsafe(bname);

            //     std::array<int,2> dist = { 0, 0};
            //     // auto its = SplitStringBySpace(lines[i]);
            //     // auto aid = string_pool_.GetIdByStringUnsafe(its[names[0]]);
            //     // auto bid = string_pool_.GetIdByStringUnsafe(its[names[1]]);
            //     // auto alen = std::stoi(its[len_off[0]]);
            //     // auto blen = std::stoi(its[len_off[1]]);
            //     auto gid = infos.IsAssigned(aid);
            //     if (gid >= 0) {
            //         dist[0] = gid;
            //     } else {
            //         infos.groups.back().Assign(aid, alen);
            //         dist[0] = infos.groups.size() - 1;
            //         infos.Assign(aid, dist[0]);
            //     }

            //     gid = infos.IsAssigned(bid);
            //     if (gid >= 0) {
            //         dist[1] = gid;
            //     } else {
            //         infos.groups.back().Assign(bid, blen);
            //         dist[1] = infos.groups.size() - 1;
            //         infos.Assign(bid, dist[1]);
            //     }

            //     infos.groups[dist[0]].Add(lines[i], bid);
            //     if (dist[0] == dist[1]) {
            //         infos.groups[dist[1]].Add(bid);
            //     } else {
            //         infos.groups[dist[1]].Add(lines[i], aid);
            //     }
            } 
            size = generate(lines);
        }

    };   
 

    //SaveReadName(infos.groups);

    // auto combine = [&mutex_combine](const std::vector<const std::string*> &lines) {
    //     std::lock_guard<std::mutex> lock(mutex_combine);
       
    // };
 

    MultiThreadRun((size_t)thread_size_, work);

}

void Program_SplitOverlaps::SaveReadName(const std::vector<Group>& groups) {
    auto format = [](const std::string& pattern, int d) {
        std::string result = pattern;
        std::string s = std::to_string(d);
        result.replace(pattern.find("{}"), 2, s);
        return result;
    };
    for (auto &g : groups) {
        if (!rdfname0_.empty()) {
            SaveReadName(g.reads, format(rdfname0_, g.seqid));
        }
        if (!rdfname1_.empty()) {
            SaveReadName(g.reads2, format(rdfname1_, g.seqid));
        }
    }
}

void Program_SplitOverlaps::SaveReadName(const std::unordered_set<Key>& reads, const std::string &fname) {
    GzFileWriter writer(fname);
    
    if (writer.Valid()) {
        for (const auto &n : reads) {
            writer.Write(string_pool_.QueryStringById(n));
            writer.Write("\n");
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}



void Program_SplitName::Running() {
    
    LoadReadnames();

    if (fn_ols_.empty()) {
        GroupReadsRandomly();
    } else {
        GroupReadsByOverlaps();
    }

    SaveReadnames();

    if (!fn_part_ols_.empty()) {
        SaveOverlaps(fn_part_ols_);
    }
}

void Program_SplitName::LoadReadnames() {
    LOG(INFO)("load read names");
    LoadReadFile(fn_rds_, "", [&](const SeqReader::Item& item) {
        auto id = string_pool_.GetIdByStringUnsafe(item.head);
        if ((int)lengths_.size() <= id) {
            lengths_.push_back(0);
        }
        lengths_[id] = (int)item.seq.size();
    });

    size_t minlen = 0;
    if (base_size_ > 0) {
        std::vector<int> lens = lengths_;
        FindLongestXHeap(lens, base_size_);
        minlen = lens[0];
    }
    LOG(INFO)("minimum length = %zd", minlen);

    for (size_t i = 0; i < lengths_.size(); ++i) {
        if (lengths_[i] >= (int)minlen) {
            read_ids_index_[i] = read_ids_.size();
            read_ids_.push_back(i);
        }
    }
    LOG(INFO)("read size = %zd", read_ids_.size());
}


void Program_SplitName::GroupReadsRandomly() {
    groups_.assign(read_ids_.size(), -1);

    auto ig = 0;
    long long int gsize = 0;
    for (size_t i = 0; i < groups_.size(); ++i) {
        groups_[i] = ig;
        gsize += lengths_[read_ids_[i]];

        if (gsize >= block_size_) {
            ig += 1;
            gsize = 0;
        }
    }
}

void Program_SplitName::GroupReadsByOverlaps() {
    std::string itype = OverlapStore::DetectFileType(fn_ols_);
    assert(itype == "txt" || itype == "txt.gz" || itype == "paf" || itype == "paf.gz" );

    std::vector<std::string> ifnames = itype == "txt" ? GetLineFromFile(fn_ols_) : std::vector<std::string>({fn_ols_});
    MultiFileReader reader(ifnames);

    std::vector<std::vector<Seq::Id>> groups(read_ids_.size());
    
    std::mutex mutex_combine;
    std::mutex mutex_generate;


    // auto generate = [&reader, &mutex_generate](std::vector<char> &block) {
    //     std::lock_guard<std::mutex> lock(mutex_generate);
    //     return reader.GetBlock(block, "\n", 1);
    // };

    
    auto combine = [&mutex_combine, &groups](const std::vector<std::vector<int>> &gs) {
        std::lock_guard<std::mutex> lock(mutex_combine);

        for (size_t i = 0; i < gs.size(); ++i) {
            groups[i].insert(groups[i].end(), gs[i].begin(), gs[i].end());
        }
    };
 

    auto work = [&](size_t threadid) {

        LineInBlock line_in_block(reader, 100000000, &mutex_generate);

        std::vector<std::vector<int>> gs(read_ids_.size());
        std::string line;
        while (line_in_block.GetLine(line)) {
            auto its = SplitStringBySpace(line);
            // TODO 目前只支持paf
            auto aid = string_pool_.QueryIdByString(its[0]);
            auto bid = string_pool_.QueryIdByString(its[5]);

            auto ait = read_ids_index_.find(aid);
            auto bit = read_ids_index_.find(bid);
            if (ait != read_ids_index_.end() && bit != read_ids_index_.end()) {
                gs[ait->second].push_back(bit->second);
                gs[bit->second].push_back(ait->second);
            }
        }
  
        combine(gs);
        for (auto &ig : gs) {
            ig.clear();
        }
    };   
 
    MultiThreadRun((size_t)thread_size_, work);

    groups_.assign(read_ids_.size(), -1);

    int curr_group = 0;
    long long int curr_group_size = 0;

    for (size_t i = 0; i < groups_.size(); ++i) {
        if (groups_[i] == -1) {
            groups_[i] = curr_group;
            curr_group_size += lengths_[read_ids_[i]];

            for (auto j : groups[i]) {
                if (groups_[j] == -1) {
                    groups_[j] = curr_group;
                    curr_group_size += lengths_[read_ids_[j]];
                }
            }
            if (curr_group_size > block_size_) {
            printf("%lld\n", curr_group_size);
                curr_group += 1;
                curr_group_size = 0;
            }
        }
    }
}

void Program_SplitName::SaveReadnames() {
    auto maxitem = std::max_element(groups_.begin(), groups_.end());
    for (int ig = 0; ig <= *maxitem; ++ig) {
        auto fn = Format(fn_part_rdnames_, ig);
        std::ofstream of(fn);
        if (of.is_open()) {
            for (size_t i = 0; i < groups_.size(); ++i) {
                if (groups_[i] == ig) {
                    of << string_pool_.QueryStringById(read_ids_[i]) << "\n";
                }
            }
        }
    }
}


void Program_SplitName::SaveOverlaps(const std::string &opattern) {
    auto maxitem = std::max_element(groups_.begin(), groups_.end());

    std::vector<std::ofstream> ofs(*maxitem + 1);
    for (size_t i = 0; i < ofs.size(); ++i) {
        auto fn = Format(opattern, i);
        ofs[i].open(fn);
        if (!ofs[i]) {
            LOG(ERROR)("Failed to open file: %s", fn.c_str());
        }
    }

    std::string itype = OverlapStore::DetectFileType(fn_ols_);
    assert(itype == "txt" || itype == "txt.gz" || itype == "paf" || itype == "paf.gz" );

    std::vector<std::string> ifnames = itype == "txt" ? GetLineFromFile(fn_ols_) : std::vector<std::string>({fn_ols_});
    MultiFileReader reader(ifnames);
    
    std::mutex mutex_combine;
    std::mutex mutex_generate;

    auto combine = [&ofs, &mutex_combine](std::vector<std::ostringstream>& osss, size_t flushsize=0) {
        std::lock_guard<std::mutex> lock(mutex_combine);

        assert(osss.size() == ofs.size());
        for (size_t i = 0; i < osss.size(); ++i) {
            if (osss[i].tellp() >= (int)flushsize) {
                ofs[i] << osss[i].str();
                osss[i].str("");
            }
        }
    };
 

    auto work = [&](size_t tid) {

        const size_t BSIZE = 100000000;
        LineInBlock line_in_block(reader, BSIZE, &mutex_generate);

        std::vector<std::ostringstream> osss(ofs.size());

        std::string line;
        while (line_in_block.GetLine(line)) {
            auto its = SplitStringBySpace(line);
            // TODO 目前只支持paf
            auto aid = string_pool_.QueryIdByString(its[0]);
            auto bid = string_pool_.QueryIdByString(its[5]);

            int writed = -1;
            auto ait = read_ids_index_.find(aid);
            auto bit = read_ids_index_.find(bid);
            
            if (ait != read_ids_index_.end()) {
                osss[groups_[ait->second]] << line << "\n";
                writed = groups_[ait->second];
            }

            if (bit != read_ids_index_.end()) {
                if (writed != groups_[bit->second]) {
                    osss[groups_[bit->second]] << line << "\n";
                }
            }

            combine(osss, BSIZE);

        }
        combine(osss, 0);

    };   
 
    MultiThreadRun((size_t)thread_size_, work);
}

void Program_SplitName::Group(const std::vector<std::vector<int>>& groups, const ReadStore &rd_store) {
    std::vector<int> work(groups.size(), -1);

    int curr_group = 0;
    int curr_group_size = 0;

    for (size_t i = 0; i < groups.size(); ++i) {
        if (work[i] == -1) {
            work[i] = curr_group;
            curr_group_size += rd_store.GetSeqLength(i);

            for (auto j : groups[i]) {
                if (j < (int)work.size()) { // max_read
                    if (work[j] == -1) {
                        work[j] = curr_group;
                        curr_group_size += rd_store.GetSeqLength(j);
                    }
                }
            }

            if (curr_group_size > block_size_) {
                curr_group += 1;
                curr_group_size = 0;
            }
        }
    }
}

void Program_SplitOverlaps2::Running() {
    // 只支持 paf
    
    std::string itype = OverlapStore::DetectFileType(ifname_);
    assert(itype == "txt" || itype == "txt.gz" || itype == "paf" || itype == "paf.gz" );

    std::array<int,2> names {0, 5};
    std::array<int, 2> len_off { 1, 6};

    std::vector<std::string> ifnames = itype == "txt" ? GetLineFromFile(ifname_) : std::vector<std::string>({ifname_});
    std::shared_ptr<GzFileReader> ifile;
    size_t ifile_index = 0;    
 

    ReadStore rd_store(string_pool_);
    rd_store.Load(rd_fname_, "", false);
    string_pool_.Save("id2name");
    std::mutex mutex_combine;
    std::mutex mutex_generate;

    Infos infos(ofname_, sub_size_, string_pool_.Size());
    LoadSubNames(rdfname0_, sub_size_, infos);

    auto generate = [&ifile, &ifnames, &ifile_index,&mutex_generate](std::vector<std::string> &lines) {
        std::lock_guard<std::mutex> lock(mutex_generate);
        auto size = ifile == nullptr ? 0 : ifile->GetLines(lines);
        while (size == 0) {
            if (ifile_index < ifnames.size()) {
                ifile.reset(new GzFileReader(ifnames[ifile_index++]));
                if (ifile->Valid()) {
                    size = ifile->GetLines(lines);
                } else {
                    LOG(ERROR)("Failed to open file: %s", ifnames[ifile_index-1].c_str());
                }
            } else {
                break;
            }
        }
        return size;
    };

    auto combine = [&mutex_combine, &infos](size_t size, WorkArea &wa, const std::vector<std::string>& lines) {
        std::lock_guard<std::mutex> lock(mutex_combine);
        infos.Merge(size, wa, lines);
    };
 
    // auto get_name_len_from_paf_line = [](const std::string& line, char* a, int& alen, char* b, int &blen) {
    //     int astart, aend, bstart, bend;
    //     char buf0[1000], buf1[1000]; 
    //     sscanf(line.c_str(), "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s", a, &alen, &astart, &aend, buf0, b, &blen, &bstart, &bend, buf1);
    // };
    

    const int line_size = 100000;
    std::atomic<int> finished { 0 };

    auto work = [&](size_t threadid) {
        std::vector<std::string> lines(line_size);
        WorkArea workarea(line_size);

        size_t size = generate(lines);
        while (size > 0) {
            finished += size;
            if (finished % 100000000 == 0) {
                LOG(INFO)("Fininshed: %d ", finished.load());
            }

            for (size_t i=0; i<size; ++i) {
                // char aname[1024];
                // char bname[1024];
                // int alen;
                // int blen;
                // get_name_len_from_paf_line(lines[i], aname, alen, bname, blen);
                // auto aid = string_pool_.GetIdByStringUnsafe(aname);
                // auto bid = string_pool_.GetIdByStringUnsafe(bname);

                auto its = SplitStringBySpace(lines[i]);
                auto aid = string_pool_.GetIdByStringUnsafe(its[names[0]]);
                auto bid = string_pool_.GetIdByStringUnsafe(its[names[1]]);
                auto alen = std::stoi(its[len_off[0]]);
                auto blen = std::stoi(its[len_off[1]]);
                std::array<int,2> dist = { infos.GetGroup(aid), infos.GetGroup(bid)};

                workarea.Set(i, aid, alen, dist[0], bid, blen, dist[1]);
            } 
            combine(size, workarea, lines);
            size = generate(lines);
        }

    };   
    

    MultiThreadRun((size_t)thread_size_, work);
    SaveReadName(infos.groups);

}

void Program_SplitOverlaps2::LoadSubNames(const std::string& pattern, size_t subsize, Infos &infos) {
    for (size_t i = 0; i < subsize; ++i) {
        auto fname = Format(pattern, i);

        std::ifstream file(fname);
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                infos.SetGroup(string_pool_.GetIdByStringUnsafe(line), i);

            }
        } else {
            LOG(ERROR)("Failed to open file: %s", fname.c_str());
        }
            
    }
}

void Program_SplitOverlaps2::SaveReadName(const std::vector<Group>& groups) {
    for (auto &g : groups) {
        if (!rdfname1_.empty()) {  
            SaveReadName(g.reads2, Format(rdfname1_, g.seqid));
        }
    }
}

void Program_SplitOverlaps2::SaveReadName(const std::unordered_set<int>& reads, const std::string &fname) {
    GzFileWriter writer(fname);
    
    if (writer.Valid()) {
        for (const auto &n : reads) {
            writer.Write(string_pool_.QueryStringById(n));
            writer.Write("\n");
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}

void Program_Test::Running() {

    ReadStore rd_store;
    rd_store.Load(ifname_, "", false);
    LOG(INFO)("Load data");
    rd_store.Load(ifname_, "", true);

}

void Program_Hic::Running() {
    StringPool string_pool;
    HicReadInfos hic_infos(string_pool);
    hic_infos.Build(fn_hic1_, fn_paf1_, fn_hic2_, fn_paf2_, fn_vars_);
    hic_infos.Save(fn_snp_in_hic_);

 

    // LOG(INFO)("Start grouping ols");
    // std::unordered_map<Seq::Id, std::unordered_map<Seq::Id, std::vector<const Overlap*>>> groups;
    // ol_store.GroupQuery(groups, 1);

    // return;
    // LOG(INFO)("Start load tiles");
    // auto tiles = LoadTiles(ifname1_, string_pool);

    // std::unordered_map<Seq::Id, std::vector<const Overlap*>> seg_ols_begin;
    // std::unordered_map<Seq::Id, std::vector<const Overlap*>> seg_ols_end;
    // for (const auto & t : tiles) {
    //     const auto &reads = t.second;
    //     for (size_t i = 0; i < std::min<size_t>(3, reads.size()); ++i) {
    //         auto &r = reads[i];
    //         auto iter = groups.find(r);
    //         if (iter != groups.end()) {
    //             for (auto &ols : iter->second) {
    //                 auto &one_seg = seg_ols_begin[t.first];
    //                 one_seg.insert(one_seg.end(), ols.second.begin(), ols.second.end());
    //             }
    //         }
    //     }

    //     for (size_t i = reads.size() - std::min<size_t>(3, reads.size()); i < reads.size(); ++i) {
    //         auto &r = reads[i];
    //         auto iter = groups.find(r);
    //         if (iter != groups.end()) {
    //             for (auto &ols : iter->second) {
    //                 auto &one_seg = seg_ols_end[t.first];
    //                 one_seg.insert(one_seg.end(), ols.second.begin(), ols.second.end());
    //             }
    //         }
    //     }
    // }

    // LOG(INFO)("Start load hic mappings");
    // auto hic1 = IndexHiCMapping(ifname2_, string_pool);
    // auto hic2 = IndexHiCMapping(ifname3_, string_pool);
    
    // std::unordered_map<Seq::Id, std::unordered_set<Seq::Id>> seg_hic1_begin;
    // std::unordered_map<Seq::Id, std::unordered_set<Seq::Id>> seg_hic1_end;
    // std::unordered_map<Seq::Id, std::unordered_set<Seq::Id>> seg_hic2_begin;
    // std::unordered_map<Seq::Id, std::unordered_set<Seq::Id>> seg_hic2_end;
    // for (auto seg : seg_ols_begin) {
    //     for (auto ol : seg.second) {
    //         // 
    //         auto iter1 = hic1.find(ol->b_.id);
    //         if (iter1 != hic1.end()) {
    //             auto& shic = seg_hic1_begin[seg.first];
    //             const size_t N = 1000;
    //             for (size_t i = ol->b_.start / N; i < (ol->b_.end - 1) / N + 1; ++i) {
    //                 shic.insert(iter1->second[i].begin(), iter1->second[i].end());
    //             }
    //         }

    //         auto iter2 = hic2.find(ol->b_.id);
    //         if (iter2 != hic2.end()) {
    //             auto& shic = seg_hic2_begin[seg.first];
    //             const size_t N = 1000;
    //             for (size_t i = ol->b_.start / N; i < (ol->b_.end - 1) / N + 1; ++i) {
    //                 shic.insert(iter2->second[i].begin(), iter2->second[i].end());
    //             }
    //         }
    //     }
    // }

    // for (auto seg : seg_ols_end) {
    //     for (auto ol : seg.second) {
    //         // 
    //         auto iter1 = hic1.find(ol->b_.id);
    //         if (iter1 != hic1.end()) {
    //             auto& shic = seg_hic1_end[seg.first];
    //             const size_t N = 1000;
    //             for (size_t i = ol->b_.start / N; i < (ol->b_.end - 1) / N + 1; ++i) {
    //                 shic.insert(iter1->second[i].begin(), iter1->second[i].end());
    //             }
    //         }

    //         auto iter2 = hic2.find(ol->b_.id);
    //         if (iter2 != hic2.end()) {
    //             auto& shic = seg_hic2_end[seg.first];
    //             const size_t N = 1000;
    //             for (size_t i = ol->b_.start / N; i < (ol->b_.end - 1) / N + 1; ++i) {
    //                 shic.insert(iter2->second[i].begin(), iter2->second[i].end());
    //             }
    //         }
    //     }
    // }

    // LOG(INFO)("Scoring ctgs");
    // for (auto seg0 : seg_hic1_begin) {
    //     auto s0 = seg0.first;
    //     for (auto seg1 : seg_hic1_begin) {
    //         auto s1 = seg1.first;
    //         if (s0 >= s1) continue;
    //         size_t count = CountIntersect(seg_hic1_begin[s0], seg_hic2_begin[s1]) + CountIntersect(seg_hic2_begin[s0], seg_hic1_begin[s1]);

    //         printf("%s:B <-> %s:B -> %zd\n", string_pool.QueryStringById(s0).c_str(), string_pool.QueryStringById(s1).c_str(), count);
    //     }

    //     for (auto seg1 : seg_hic1_end) {
    //         auto s1 = seg1.first;
    //         if (s0 >= s1) continue;
    //         size_t count = CountIntersect(seg_hic1_begin[s0], seg_hic2_end[s1]) + CountIntersect(seg_hic2_begin[s0], seg_hic1_end[s1]);

    //         printf("%s:B <-> %s:E -> %zd\n", string_pool.QueryStringById(s0).c_str(), string_pool.QueryStringById(s1).c_str(), count);
    //     }
    // }

    // for (auto seg0 : seg_hic1_end) {
    //     auto s0 = seg0.first;
    //     for (auto seg1 : seg_hic1_begin) {
    //         auto s1 = seg1.first;
    //         if (s0 >= s1) continue;
    //         size_t count = CountIntersect(seg_hic1_end[s0], seg_hic2_begin[s1]) + CountIntersect(seg_hic2_end[s0], seg_hic1_begin[s1]);

    //         printf("%s:E <-> %s:B -> %zd\n", string_pool.QueryStringById(s0).c_str(), string_pool.QueryStringById(s1).c_str(), count);
    //     }

    //     for (auto seg1 : seg_hic1_end) {
    //         auto s1 = seg1.first;
    //         if (s0 >= s1) continue;
    //         size_t count = CountIntersect(seg_hic1_end[s0], seg_hic2_end[s1]) + CountIntersect(seg_hic2_end[s0], seg_hic1_end[s1]);

    //         printf("%s:E <-> %s:E -> %zd\n", string_pool.QueryStringById(s0).c_str(), string_pool.QueryStringById(s1).c_str(), count);
    //     }
    // }
}


std::unordered_map<Seq::Id, std::vector<std::unordered_set<Seq::Id>>> IndexHiCMapping(const std::string &fname, StringPool& sp) {

    std::unordered_map<Seq::Id, std::vector<std::unordered_set<Seq::Id>>> index;
    size_t i = 0;

    GzFileReader ols(fname);
    if (ols.Valid()) {
        for (std::string line = ols.GetNoEmptyLine(); !line.empty(); line = ols.GetNoEmptyLine()) {
            i++;
            auto its = SplitStringBySpace(line);    // paf format

            auto ctgid = sp.GetIdByString(its[5]);
            auto rdid = sp.GetIdByString(its[0]);
            auto aligned = std::stoul(its[3]) - std::stoul(its[2]);
            if (aligned < 100) continue;

            const size_t N = 1000;
            if (index.find(ctgid) == index.end()) {
                size_t ctglen = std::stoul(its[6]);
                index[ctgid].assign(ctglen + N-1, std::unordered_set<Seq::Id>());
            }

            auto start = std::stoul(its[7]) / N;
            auto end = (std::stoul(its[8])  - 1) / N;
            index[ctgid][start].insert(rdid);
            if (end != start) {
                index[ctgid][end].insert(rdid);
            }


            if (i % 1000000 == 0) {
                LOG(INFO)("%zd", i);
            }
        }
    } else {
        LOG(INFO)("Failed to load file: %s", fname.c_str());
    }

    return index;

}

size_t CountIntersect(const std::unordered_set<Seq::Id>& s0, const std::unordered_set<Seq::Id> &s1) {
    size_t count = 0;
    for (auto i : s0) {
        if (s1.find(i) != s1.end()) {
            count++;
        }
    }
    return count;
}

} // namespace fsa
