#include "misc_tools.hpp"

#include "overlap_store.hpp"
#include "read_store.hpp"

#include <sstream>

namespace fsa {


void Program_SplitOverlaps::Running() {
    // 只支持 paf
    
    std::string itype = OverlapStore::DetectFileType(ifname_);
    assert(itype == "txt" || itype == "txt.gz" || itype == "paf" || itype == "paf.gz" );

    std::array<int,2> names {0, 5};
    std::array<int, 2> len_off { 1, 6};

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

    const int line_size = 100000;
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

    auto combine = [&mutex_combine](const std::vector<const std::string*> &lines) {
        std::lock_guard<std::mutex> lock(mutex_combine);
       
    };
 

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
        if (lengths_.size() <= id) {
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
        if (lengths_[i] >= minlen) {
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


    auto generate = [&reader, &mutex_generate](std::vector<char> &block) {
        std::lock_guard<std::mutex> lock(mutex_generate);
        return reader.GetBlock(block, "\n", 1);
    };

    
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
            printf("%zd\n", curr_group_size);
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
            if (osss[i].tellp() >= flushsize) {
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
                if (j < work.size()) { // max_read
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
 
    auto get_name_len_from_paf_line = [](const std::string& line, char* a, int& alen, char* b, int &blen) {
        int astart, aend, bstart, bend;
        char buf0[1000], buf1[1000]; 
        sscanf(line.c_str(), "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s", a, &alen, &astart, &aend, buf0, b, &blen, &bstart, &bend, buf1);
    };
    

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
} // namespace fsa
