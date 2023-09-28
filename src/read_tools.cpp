#include "read_tools.hpp"

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "sequence_store.hpp"
#include "file_io.hpp"
#include "utils/project_file.hpp"

#include "phase/hic_read_infos.hpp"

namespace fsa {

inline std::string Format(const std::string& pattern, int d) {
            std::string result = pattern;
            std::string s = std::to_string(d);
            result.replace(pattern.find("{}"), 2, s);
            return result;
        };

void Program_Check::Running() {
    DnaSerialTable table;
    LoadReadFile(ifname_, "", [&table](const SeqReader::Item& item) {
        for (auto c : item.seq) {
            if (table[c] == (decltype(table[c]))-1) {
                LOG(ERROR)("Found bad base %c(%d) in %s", c, (int)c, item.head.c_str());
            }
        }
    });
}


void Program_N50::Running() {
    LOG(INFO)("Load read file %s", ifname_.c_str());

    std::vector<int> read_lens;
    LoadReadFile(ifname_, "", [&read_lens, this](const SeqReader::Item& item) {
        if (DnaSeq::Check(item.seq)) {
            read_lens.push_back((int)item.seq.size());
        } else {
            LOG(WARNING)("Found bad base in %s", item.head.c_str());
        }
    });

    std::sort(read_lens.begin(), read_lens.end(), [](int a, int b) { return a > b; });

    long long total_length = std::accumulate(read_lens.begin(), read_lens.end(), (long long)0);
    long long genome_size = genome_size_ == 0 ? total_length : genome_size_;

    std::cout << "Genome: " << genome_size << "\n";
    std::cout << "Count: " << read_lens.size() << "\n";
    std::cout << "Total: " << total_length << "\n";
    std::cout << "Max: " << read_lens.front() << "\n";
    std::cout << "Min: " << read_lens.back() << "\n";
    
    long long accu = 0;
    int ns[] = { 25, 50, 75};
    size_t ins = 0;
    for (size_t i=0; i<read_lens.size(); ++i) {
        accu += read_lens[i];
        for (; ins < sizeof(ns) / sizeof(ns[0]) && accu > (long long)ns[ins]*1.0/100 * genome_size; ++ins) {
            std::cout << "N" << ns[ins] << ": " << read_lens[i] << "\n";
            std::cout << "L" << ns[ins] << ": " << i+1 << "\n";
        }
        if (ins >= sizeof(ns) / sizeof(ns[0])) {
            break;
        }
    }
}

long long FileLength(const std::string &fname) {
    std::ifstream in(fname);
	in.seekg(0, std::ios::end);
	return in.tellg(); 
}

void Program_Split::Running() {
    
    if (block_size_ == 0) LOG(ERROR)("block_size should be greater than 0");

    auto format = [](const std::string& pattern, int d) {
        std::string result = pattern;
        std::string s = std::to_string(d);
        result.replace(pattern.find("{}"), 2, s);
        return result;
    };

    size_t min_length = 0;
    if (base_size_ > 0) {

        std::vector<int> lengths;
        LoadReadFile(ifname_, "", [&lengths](const SeqReader::Item& item) {
            lengths.push_back((int)item.seq.size());
        });

        LOG(INFO)("length size = %zd", lengths.size());
    
        FindLongestXHeap(lengths, base_size_);
        min_length = lengths[0];
    }

    int index = 0;
    FastaWriter *writer = nullptr;
    long long accu = 0;
    
    std::vector<int> lengths;
    LoadReadFile(ifname_, "", [&](const SeqReader::Item& item) {
        if (writer == nullptr) {
            writer = new FastaWriter(format(opattern_, index++));
        }

        assert(accu < block_size_);
        
        if (item.seq.size() >= min_length)  {
            writer->Write(item);
            accu += item.seq.size();
        } 
        
        if (accu >= block_size_) {
            delete writer;
            writer = nullptr;
            accu = 0;
        }
    });

    if (writer != nullptr) {
        delete writer;
    }

}

// void Program_SplitName::Running() {
    
//     if (block_size_ == 0) LOG(ERROR)("block_size should be greater than 0");

//     auto format = [](const std::string& pattern, int d) {
//         std::string result = pattern;
//         std::string s = std::to_string(d);
//         result.replace(pattern.find("{}"), 2, s);
//         return result;
//     };

//     size_t min_length = 0;
//     if (base_size_ > 0) {

//         std::vector<int> lengths;
//         LoadReadFile(ifname_, "", [&lengths](const SeqReader::Item& item) {
//             lengths.push_back((int)item.seq.size());
//         });

//         LOG(INFO)("length size = %zd", lengths.size());
    
//         FindLongestXHeap(lengths, base_size_);
//         min_length = lengths[0];
//     }

//     int index = 0;
//     GzFileWriter *writer = nullptr;
//     long long accu = 0;
    
//     std::vector<int> lengths;
//     LoadReadFile(ifname_, "", [&](const SeqReader::Item& item) {
//         if (writer == nullptr) {
//             writer = new GzFileWriter(format(opattern_, index++));
//         }

//         assert(accu < block_size_);

//         if (item.seq.size() >= min_length)  {
//             writer->Write(item.head);
//             writer->Write("\n");
//             accu += item.seq.size();
//         } 

//         if (accu >= block_size_) {
//             delete writer;
//             writer = nullptr;
//             accu = 0;

//         }
//     });

//     if (writer != nullptr) {
//         delete writer;
//     }   
// }


void Program_SplitName::Running() {
    
    LoadReadnames(fn_reads_);

    if (method_ == "plain") {
        GroupReadsPlainly();
    } else {
        GroupReadsByOverlaps();
    }

    LOG(INFO)("Save read names");
    SaveReadnames(fn_sub_reads_);

    if (!fn_sub_ols_.empty() && !fn_ols_.empty()) {
        LOG(INFO)("Save overlaps");
        SaveOverlaps(fn_ols_, fn_sub_ols_);
    }
}

void Program_SplitName::LoadReadnames(const std::string& fname) {
    LOG(INFO)("load read names");
    LoadReadFile(fname, "", [&](const SeqReader::Item& item) {
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

void Program_SplitName::GroupReadsPlainly() {
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

    auto combine = [&mutex_combine, &groups](const std::vector<std::vector<int>> &gs) {
        std::lock_guard<std::mutex> lock(mutex_combine);

        for (size_t i = 0; i < gs.size(); ++i) {
            groups[i].insert(groups[i].end(), gs[i].begin(), gs[i].end());
        }
    };
 

    auto work = [&](size_t threadid) {

        LineInBlock line_in_block(reader, 10000000, &mutex_generate);

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

void Program_SplitName::SaveReadnames(const std::string &opattern) {
    auto maxitem = std::max_element(groups_.begin(), groups_.end());
    for (int ig = 0; ig <= *maxitem; ++ig) {
        auto fn = Format(opattern, ig);
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


void Program_SplitName::SaveOverlaps(const std::string &fn_ols, const std::string &fn_sub_ols) {
    auto maxitem = std::max_element(groups_.begin(), groups_.end());

    std::vector<std::ofstream> ofs(*maxitem + 1);
    for (size_t i = 0; i < ofs.size(); ++i) {
        auto fn = Format(fn_sub_ols, i);
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
            if ((size_t)osss[i].tellp() >= flushsize) {
                ofs[i] << osss[i].str();
                osss[i].str("");
            }
        }
    };
 

    auto work = [&](size_t tid) {

        const size_t BSIZE = 1000000;
        LineInBlock line_in_block(reader, BSIZE, &mutex_generate);

        std::vector<std::ostringstream> osss(ofs.size());

        size_t lineno = 0;
        
        std::string line;
        for (bool is_valid = line_in_block.GetLine(line); is_valid; is_valid = line_in_block.GetLine(line)) {
            lineno ++;

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

            if (lineno % 100000 == 0) {
                combine(osss, BSIZE);
            } 


        }
        combine(osss, 0);

    };   
 
    MultiThreadRun(std::min<size_t>(8, thread_size_), work);
}

void Program_Longest::Running() {
    int min_length = min_length_;
    if (base_size_ > 0) {
        std::vector<int> lengths;
        LoadReadFile(ifname_, "", [&lengths, this](const SeqReader::Item& item) {
            if (DnaSeq::Check(item.seq)) {
                lengths.push_back((int)item.seq.size());
            } else {
                LOG(WARNING)("Found bad base in %s", item.head.c_str());
            }
        });

        LOG(INFO)("length size = %zd", lengths.size());
        
        FindLongestXHeap(lengths, base_size_);
        if (lengths[0] > min_length)  min_length = lengths[0];
    }

    LOG(INFO)("min_length = %zd", min_length);

    FilterReadFile(ifname_, ofname_, id2name_, [min_length, this](SeqReader::Item& item) {
        return (int)item.seq.size() >= min_length && DnaSeq::Check(item.seq);
    });
}

void Program_Random::Running() {
    assert(base_size_ > 0);
    long long total = 0;
    LoadReadFile(ifname_, "", [&total, this](const SeqReader::Item& item) {
        if ((int)item.seq.size() >= min_length_) {
            if (DnaSeq::Check(item.seq)) {
                total += item.seq.size();
            } else {
                LOG(WARNING)("Found bad base in %s", item.head.c_str());
            }
        }
    });

    double rate = base_size_ * 1.01 / total;    
    LOG(INFO)("size = %lld, rate = %f", total, rate);

    srand(clock());
    long long accu = 0;
    FilterReadFile(ifname_, ofname_, id2name_, [&accu, this, rate](SeqReader::Item& item) {
        auto random = []() -> double {
            const int N = 10000;
            return rand() % (N+1) * 1.0 / N;
        };

        bool r = accu < base_size_ && (int)item.seq.size() >= min_length_ && (DnaSeq::Check(item.seq)) && rate >= random();
        if (r) accu += item.seq.size();
        return r;
    });

}


void Program_Sub::Running() {
    std::unordered_set<std::string> names;
    if (!names_.empty()) {
        auto ns = SplitStringByChar(names_, ',');
        names.insert(ns.begin(), ns.end());
    }
    if (!names_fname_.empty()) {
        GzFileReader names_file(names_fname_);
        if (names_file.Valid()) {
            std::string line = names_file.GetNoEmptyLine();
            while (!line.empty()) {
                names.insert(line);
                line = names_file.GetNoEmptyLine();
            }
        } else {
            LOG(ERROR)("Failed to open file: %s", names_.c_str());
        }
    }

    LOG(INFO)("subreads size: %zd", names.size());
    FilterReadFile(ifname_, ofname_, "", [&names](SeqReader::Item& item) {
        return names.find(item.head) != names.end();
    });
}


void Program_Test::Running() {
    std::unordered_set<Seq::Id> ids;
    StringPool sp;
    ReadStore rd(sp);

    std::ifstream file(ifname1_);
    std::string line;
    while (std::getline(file, line)) {
        ids.insert(sp.GetIdByStringUnsafe(line));
    }

    rd.Load(ifname0_, "", true, ids);

    printf("Load %zd, %zd\n", ids.size(), rd.Size());
}

template<typename C>
void FilterReadFile(const std::string &ifname, const std::string &ofname, const std::string &id2name, C check) {
    std::shared_ptr<GzFileWriter> i2n(id2name.empty() ? nullptr : new GzFileWriter(id2name));

    
    std::string t = ReadStore::DetectFileType(ofname);
    SeqWriter* writer = nullptr;
    if (t == "fastq" || t == "fastq.gz") {
        writer = new FastqWriter(ofname);
    } else if (t == "fasta" || t == "fasta.gz") {
        writer = new FastaWriter(ofname);
    } else {
        LOG(ERROR)("Failed to recognize filetype: %s", ofname.c_str());
    }

    size_t index = 0;
    LoadReadFile(ifname, "", [&writer, i2n, check, &index](SeqReader::Item& item) { 
        if (i2n != nullptr) {
            char buff[24];
            sprintf(buff, "%zd ", index);
            i2n->Write(buff);
            i2n->Write(item.head);
            i2n->Write("\n");
            
            item.head = buff;   
        }

        if (check(item)) {
            writer->Write(item);
        }

        index++;
    });

    delete writer;
}

} // namespace fsa
