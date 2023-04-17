#include "overlap_tools.hpp"

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "utils/project_file.hpp"
#include "assemble/read_variants.hpp"

namespace fsa {

void Program_Filter::Running() {
    OverlapStore ol_store;
    
    std::string itype = OverlapStore::DetectFileType(ifname_);
    std::string otype = OverlapStore::DetectFileType(ofname_);

    std::vector<std::string> ifnames = itype == "txt" ? GetLineFromFile(ifname_) : std::vector<std::string>({ifname_});
    itype = OverlapStore::DetectFileType(ifnames[0]);

    std::shared_ptr<GzFileReader> ifile;
    size_t ifile_index = 0;
    GzFileWriter ofile(ofname_);
    std::shared_ptr<GzFileWriter> flt_file(filtered_.empty() ? (GzFileWriter*)nullptr : new GzFileWriter(filtered_));

    StringPool string_pool;
    LOG(INFO)("Load inconsistent pairs");
    PhaseInfoFile phased_info(string_pool);
    phased_info.Load(inconsistent_);

    LOG(INFO)("Load consistent pairs");
    PhaseInfoFile consistent(string_pool);
    consistent.Load(consistent_);

    if (!range_fn_.empty()) {
        LOG(INFO)("Load range of corrected reads");
        LoadRanges(range_fn_, string_pool);
    }

    std::mutex mutex_combine;
    std::mutex mutex_generate;
    std::array<long long, 2> done = {0,0};
    int block_size = 1000;

    auto generate = [&mutex_generate, &ifile, &ifnames, &ifile_index, &done](std::vector<std::string> &lines) {
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

        done[0] += size;
        if (done[0] >= done[1]) {
             LOG(INFO)("Done %lld", done[1]);
             done[1] += 500000 ;
        }
        return size;
    };

    auto combine = [&ofile, &flt_file, &mutex_combine](const std::vector<const std::string*> &keeped, const std::vector<const std::string*> &filtered) {
        std::lock_guard<std::mutex> lock(mutex_combine);
        for (auto l : keeped) {
            ofile.Write(*l);
        }
        if (flt_file) {
            for (auto l : filtered) {
                flt_file->Write(*l);
            }
        }
    };

    auto work = [&](size_t threadid) {
        std::vector<std::string> lines(block_size);
        decltype(&OverlapStore::FromPafLine) from_line = itype == "paf" ? &OverlapStore::FromPafLine : &OverlapStore::FromM4aLine;

        size_t size = generate(lines);

        StringPool::UnsafeNameId nameid(phased_info.GetStringPool());

        while (size > 0) {
            std::vector<const std::string*> keeped;
            std::vector<const std::string*> filtered;
            for (size_t i=0; i<size; ++i) {
                Overlap o;
                if (from_line(lines[i], o, nameid)) {
                    // if (consistent.Contain(o, threshold_) || (!phased_info.IsRemoved(o.a_.id) && !phased_info.IsRemoved(o.b_.id) &&
                    //     !phased_info.Contain(o, threshold_, strict_))) {
                            
                    //     keeped.push_back(&lines[i]);
                    // } else {
                    //     filtered.push_back(&lines[i]);
                    // }
                    // auto distance_consistent = consistent.MinDistance(o);
                    // auto distance_inconsistent = phased_info.MinDistance(o);
                    // if (distance_consistent < PhaseInfoFile::MAX_DISTANCE && distance_consistent <= distance_inconsistent + 1000) {
                    //     keeped.push_back(&lines[i]);
                    // } else {
                    //     if (distance_inconsistent < PhaseInfoFile::MAX_DISTANCE && (threshold_ == -1 || distance_inconsistent <= threshold_)) {
                    //         filtered.push_back(&lines[i]);
                    //     } else {
                    //         keeped.push_back(&lines[i]);
                    //     }
                    // }

                    int off_1to0 = o.MappingToSource<1>({0})[0];
                    int off_0to1 = o.MappingToTarget<1>({0})[0];
                    auto r0 = ranges_.find(o.a_.id);
                    auto r1 = ranges_.find(o.b_.id);
                    if (r0 != ranges_.end()) {
                        if (o.SameDirect()) {
                            off_1to0 += r0->second[0];
                            off_0to1 -= r0->second[0];
                        } else {
                            off_1to0 += r0->second[0];
                            off_0to1 += r0->second[0];
                        }
                    }
                    
                    if (r1 != ranges_.end()) {
                        if (o.SameDirect()) {
                            off_1to0 -= r1->second[0];
                            off_0to1 += r1->second[0];
                        } else {
                            off_1to0 += r1->second[0];
                            off_0to1 += r1->second[0];
                        }
                    }
                    int th = std::max<int>(threshold_, rate_*std::max(off_1to0, off_0to1));
                    if (!phased_info.Contain(o.a_.id, o.b_.id, o.SameDirect(), off_0to1, off_1to0, th, strict_)) {
                        keeped.push_back(&lines[i]);
                    } else {
                        filtered.push_back(&lines[i]);
                    }
                } else {
                    keeped.push_back(&lines[i]);
                }
            }
            combine(keeped, filtered);
            size = generate(lines);
        }

    };    

    if (ofile.Valid()) {
        MultiThreadRun((size_t)thread_size_, work);
    } else {
        if (!ofile.Valid()) LOG(ERROR)("Failed to save file: %s", ofname_.c_str());
    }
}

void Program_Filter::LoadRanges(const std::string& fname, StringPool& sp) {
    GzFileReader reader(fname);
    if (reader.Valid()) {
        std::string line = reader.GetNoEmptyLine();
        while (!line.empty()) {
            if (line[0] == '>') {
                // >120819 range=116-116311
                auto its = SplitStringBySpace(line);
                if (its.size() >= 2) {
                    auto id = sp.GetIdByString(std::string(its[0].begin()+1, its[0].end()));
                    auto range = SplitStringByChar(its[1], '=');
                    if (range.size() == 2 && range[0] == "range") {
                        auto pos = SplitStringByChar(range[1], '-');
                        if (pos.size() == 2) {
                            size_t start = std::stoul(pos[0]);
                            size_t end = std::stoul(pos[1]);
                            ranges_[id] = { start, end };
                        }
                    }
                }

            }
            line = reader.GetNoEmptyLine();
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
}

void Program_Split::Running() {
    auto nameset_fnames = SplitStringByChar(namesets_, ',');
    auto output_fnames = SplitStringByChar(ofnames_, ',');
    if (nameset_fnames.size() != output_fnames.size() && nameset_fnames.size()+1 != output_fnames.size()) {
        LOG(ERROR)("Parameters 'ofnames' and 'namesets' do not math.");
    }

    std::vector<std::unordered_set<std::string>> namesets;
    for (auto& fn : nameset_fnames) {
        namesets.push_back(LoadNameset(fn));
    }

    std::vector<std::shared_ptr<GzFileWriter>> ofiles;
    for (auto& fn : output_fnames) {
        std::shared_ptr<GzFileWriter> f(new GzFileWriter(fn));
        if (f->Valid()) {
            ofiles.push_back(f);
        } else {
            LOG(ERROR)("Failed to save file: %s", fn.c_str());
        }
    }

    std::string itype = OverlapStore::DetectFileType(ifname_);

    std::vector<std::string> ifnames = itype == "txt" ? GetLineFromFile(ifname_) : std::vector<std::string>({ifname_});
    itype = OverlapStore::DetectFileType(ifnames[0]);
    if (itype != "paf") LOG(ERROR)("Only 'paf' format is supported");

    std::shared_ptr<GzFileReader> ifile;
    size_t ifile_index = 0;

    std::mutex mutex_combine;
    std::mutex mutex_generate;
    std::array<long long, 2> done = {0,0};
    int block_size = 1000;

    auto generate = [&mutex_generate, &ifile, &ifnames, &ifile_index, &done](std::vector<std::string> &lines) {
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

        done[0] += size;
        if (done[0] >= done[1]) {
             LOG(INFO)("Done %lld", done[1]);
             done[1] += 500000 ;
        }
        return size;
    };

    auto combine = [&ofiles, &mutex_combine](const std::vector<std::string> &lines, const std::vector<size_t> &bins, size_t sz) {
        std::lock_guard<std::mutex> lock(mutex_combine);
        for (size_t i = 0; i < sz; ++i) {
            if (bins[i] < ofiles.size()) {
                ofiles[bins[i]]->Write(lines[i]);
            }
        }
    };

    auto work = [&](size_t threadid) {
        std::vector<std::string> lines(block_size);
        std::vector<size_t>      bins(block_size);

        size_t size = generate(lines);
        
        while (size > 0) {
            for (size_t i=0; i<size; ++i) {
                auto its = SplitStringBySpace(lines[i]);

                size_t ins = 0;
                for (; ins < namesets.size(); ++ins) {
                    if (namesets[ins].find(its[0]) != namesets[ins].end() || namesets[ins].find(its[5]) != namesets[ins].end()) {
                        break;
                    }
                }
                bins[i] = ins;

            }
            combine(lines, bins, size);
            size = generate(lines);
        }

    };    

    MultiThreadRun((size_t)thread_size_, work);

}


std::unordered_set<std::string> Program_Split::LoadNameset(const std::string &fname) {
    ReadStore rdstore;
    rdstore.Load(fname, "", false);

    std::unordered_set<std::string> names;
    const auto &sp = rdstore.GetStringPool();
    for (size_t i = 0; i < sp.Size(); ++i) {
        names.insert(sp.QueryStringById(i));
    }
    return names;
}

void Program_Purge::Running() {
    StringPool sp;
    PrjContigTiles tiles(sp);

    tiles.Load(tiles_);

    ReadVariants rvs(sp);
    rvs.Load(readinfos_);
    // TODO
}


template<typename C>
void FilterFileByLine(GzFileReader &reader, GzFileWriter &writer, size_t thread_size, C check) {
    assert(reader.Valid() && writer.Valid());

    std::mutex mutex_reader;
    std::mutex mutex_writer;

    auto write_func = [&mutex_writer, &writer](std::ostringstream  &oss) {
        std::lock_guard<std::mutex> lock(mutex_writer);
        size_t size = oss.tellp();
        writer.Write(oss.str());
        writer.Flush();
        oss.str("");
        return size;
    };

    auto work_func = [&mutex_reader, &reader, write_func, check](size_t id) {
        size_t total = 0;
        std::ostringstream oss;
        LineInBlock line_in_block(reader, 10000000, &mutex_reader);
        std::string line;
        while (line_in_block.GetLine(line)) {
            if (check(line)) {
                oss << line << "\n";
            }

            if (oss.tellp() > 10000000) {
                total += write_func(oss);
            }
        }
        total += write_func(oss);
    };
    MultiThreadRun(thread_size, work_func);
}

std::vector<std::string> LoadLinesFromFile(const std::string &fname) {
    std::vector<std::string> lines;
    GzFileReader reader(fname);
    if (reader.Valid()) {
        std::string line = reader.GetNoEmptyLine();
        while (!line.empty()) {
            lines.push_back(line);
            line = reader.GetNoEmptyLine();
        }
    } else {
        LOG(ERROR)("Failed to open file: %s", fname.c_str());
    }
    return lines;
}

void Program_Sub::Running() {
    OverlapStore ol_store;
    
    std::unordered_set<std::string> names = LoadNames();
    auto is_valid_line = [&names](const std::string& line) {
        auto items = SplitStringBySpace(line); // PAF Format
        return names.find(items[0]) != names.end() || names.find(items[5]) != names.end();
    };

    auto type = OverlapStore::DetectFileType(ifname_);

    std::vector<std::string> ifnames;
    if (type == "txt") {
        ifnames = LoadLinesFromFile(ifname_);
    } else {
        ifnames.push_back(ifname_);
    }

    GzFileWriter writer(ofname_);
    for (const auto& fn : ifnames) {
        LOG(INFO)("Get overlaps from %s", fn.c_str());
        GzFileReader reader(fn);
        if (reader.Valid() && writer.Valid()) {
            FilterFileByLine(reader, writer, thread_size_, is_valid_line);
        } else {
            if (!reader.Valid()) LOG(ERROR)("Failed to open file: %s", fn.c_str());
            if (!writer.Valid()) LOG(ERROR)("Failed to open file: %s", ofname_.c_str());
        }
    }
}


std::unordered_set<std::string> Program_Sub::LoadNames() const {
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
    return names;
}


void Program_Accuracy::Running() {
    OverlapStore ol_store;
    ol_store.Load(ifname_, "", thread_size_);

    std::mutex mutex;
    
    struct Info {
        size_t len { 0 };
        size_t match { 0 };
        size_t mismatch { 0 };
        size_t clip { 0 };
        size_t insert { 0 };
        size_t dels { 0 };
        size_t hclip { 0 };
        void Print() const {
            printf("%zd, %zd, %zd, %zd, %zd, %zd, %zd\n", len, match, mismatch, clip, insert, dels, hclip);
        }
    };

    std::unordered_map<Seq::Id, Info> all_infos;

    auto merge_info = [](std::unordered_map<Seq::Id, Info>& infos, Seq::Id id, const Info& inf) {
        auto iter = infos.find(id);
        if (iter != infos.end()) {
            if (iter->second.match < inf.match) {
                iter->second = inf;
            } 
        } else {
            infos[id] = inf;
        }
    };
    
    auto combine_func = [&all_infos, &mutex, merge_info](const std::unordered_map<Seq::Id, Info>& infos) {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto &i : infos) {
            merge_info(all_infos, i.first, i.second);
        }
    };

    std::atomic<size_t> index { 0 };

    auto work_func = [&index, &ol_store, combine_func, merge_info](int tid) {
        std::unordered_map<Seq::Id, Info> infos;

        for (size_t curr = index.fetch_add(1); curr < ol_store.Size(); curr = index.fetch_add(1)) {
            auto& ol = ol_store.Get(curr);
            Info inf;

            inf.len = ol.a_.len;
            if (ol.SameDirect()) {
                inf.clip += std::min(ol.a_.start, ol.b_.start);
                inf.clip += std::min(ol.a_.len - ol.a_.end, ol.b_.len - ol.b_.end);
            } else {
                inf.clip += std::min(ol.a_.start, ol.b_.len - ol.b_.end);
                inf.clip += std::min(ol.a_.len - ol.a_.end, ol.b_.start);
            }
            inf.hclip += ol.a_.start + ol.a_.len - ol.a_.end;

            for (auto &d : ol.detail_) {
                switch (d.type) {
                case '=':
                    inf.match += d.len;
                    break;
                case 'X':
                    inf.mismatch += d.len;
                    break;
                case 'I':
                    inf.insert += d.len;
                    break;
                case 'D':
                    inf.dels += d.len;
                    break;
                default:
                    LOG(ERROR)("Not support %c", d.type);

                }
                
            }

            merge_info(infos, ol.a_.id, inf);
        }

        if (infos.size() > 0) {
            combine_func(infos);
            infos.clear();
        }

    };
    
    MultiThreadRun(thread_size_, work_func);

    LOG(INFO)("%zd", all_infos.size());
    Info accu;
    for (auto &i : all_infos) {
        accu.len += i.second.len;
        accu.match += i.second.match;
        accu.mismatch += i.second.mismatch;
        accu.clip += i.second.clip;
        accu.insert += i.second.insert;
        accu.dels += i.second.dels;
        accu.hclip += i.second.hclip;
    }

    
    size_t len_all = accu.len + accu.dels;
    size_t len_mat = len_all - accu.hclip;
    printf("size:      %zd\n", accu.len);
    printf("match:     %2.06f\t%2.06f\n", accu.match*1.0/len_mat,    accu.match*1.0/len_all);
    printf("mismatch:  %2.06f\t%2.06f\n", accu.mismatch*1.0/len_mat, accu.mismatch*1.0/len_all);
    printf("insertion: %2.06f\t%2.06f\n", accu.insert*1.0/len_mat,   accu.insert*1.0/len_all);
    printf("deletion:  %2.06f\t%2.06f\n", accu.dels*1.0/len_mat,     accu.dels*1.0/len_all);
    printf("clip:      %2.06f\t%2.06f\n", accu.hclip*1.0/len_mat,    accu.hclip*1.0/len_all);
}


void Program_Accuracy2::Running() {
    StringPool string_pool;
    OverlapStore ol_store(string_pool);
    ol_store.Load(ifname_, "", thread_size_);
    LOG(INFO)("Load overlaps: %zd", ol_store.Size());

    struct Item {
        double accu;
        const Overlap* ol; 
    };
    std::unordered_map<int, Item> best;

    for (size_t i = 0; i < ol_store.Size(); ++i) {
        const auto& ol = ol_store.Get(i);
        size_t match = 0;
        size_t del = 0;
        for (auto &d : ol.detail_) {
            switch (d.type) {
            case '=':
                match += d.len;
                break;
            case 'X':
                break;
            case 'I':
                break;
            case 'D':
                del += d.len;
                break;
            default:
                LOG(ERROR)("Not support %c", d.type);
            }
        }
        double accu = match*1.0 / (ol.a_.len + del);

        auto iter = best.find(ol.a_.id);
        if (iter != best.end()) {
            if (accu > iter->second.accu) {
                iter->second.accu = accu;
                iter->second.ol = &ol;
            }

        } else {
            best[ol.a_.id] = {accu, &ol};
        }
    }

    LOG(INFO)("Map items: %zd\n", best.size());

    struct BedItem {
        int start;
        int end;
    };
    auto load_bed = [&string_pool](const std::string& fname) {    
        GzFileReader reader(fname);

        std::unordered_map<int, std::vector<BedItem>> bedinfos;
        auto &sp = string_pool;
        if (reader.Valid()) {
            std::string line = reader.GetLine();
            while (!line.empty()) {
                auto items = SplitStringBySpace(line);
                assert(items.size() == 3);
                auto& pos = bedinfos[sp.GetIdByString(items[0])];
                int start = std::stoi(items[1]);
                int end = std::stoi(items[2]);
                pos.push_back({start, end});
                line = reader.GetLine();
            }
        }
        for (auto& i : bedinfos) {
            std::sort(i.second.begin(), i.second.end(), [](const BedItem& a, const BedItem& b) {
                return a.start < b.start;
            });
        }
        return bedinfos;
    };
 
    auto bedinfos = load_bed(bed_fname_);
    LOG(INFO)("load bed ctg size: %zd\n", bedinfos.size());

    auto find_start_bed_item = [](const std::vector<BedItem>& beds, int pos) {
        size_t s = 0;
        size_t e = beds.size();

        while (s < e) {
            size_t m = (s+e) / 2;
            if (beds[m].start > pos) {
                e = m;
            } else if (beds[m].start < pos) {
                s = m;
            } else {
                return m;
            }
            if (s == e || s + 1 == e) return s;
        }
        return e;
    };


    struct Diff {
        size_t mat {0};
        size_t mis {0};
        size_t del {0};
        size_t ins {0};
    };
        
    Diff alldiff;
    for (auto& b : best) {
        const auto &ol = *b.second.ol;
        auto &bed = bedinfos[ol.b_.id];

        std::vector<Diff> diffs(ol.b_.end - ol.b_.start + 1);
        
        size_t apos = 0;
        size_t bpos = 0;
        for (auto &d : ol.detail_) {
            switch (d.type) {
            case '=':
                for (int i = 0; i < d.len; ++i) {
                    diffs[bpos].mat++;
                    bpos++;
                    apos++;
                }
                break;
            case 'X':
                for (int i = 0; i < d.len; ++i) {
                    diffs[bpos].mis++;
                    bpos++;
                    apos++;
                }
                break;
            case 'I':
                for (int i = 0; i < d.len; ++i) {
                    diffs[bpos].ins++;
                    apos++;
                }
                break;
            case 'D':
                for (int i = 0; i < d.len; ++i) {
                    diffs[bpos].del++;
                    bpos++;
                }
                break;
            default:
                LOG(ERROR)("Not support %c", d.type);
            }
        }

        for (size_t i=1; i < diffs.size(); ++i) {
            diffs[i].mat += diffs[i-1].mat;
            diffs[i].mis += diffs[i-1].mis;
            diffs[i].del += diffs[i-1].del;
            diffs[i].ins += diffs[i-1].ins;
        }

        for (auto i = find_start_bed_item(bed, ol.b_.start); i < bed.size(); ++i) {
            //LOG(INFO)("BED: %zd (%d %d) - (%d %d), ", i, bed[i].start, bed[i].end, ol.b_.start, ol.b_.end);
            if (bed[i].end > ol.b_.start && bed[i].start < ol.b_.end ) {
                size_t s = std::max(bed[i].start, ol.b_.start) - ol.b_.start;
                size_t e = std::min(bed[i].end, ol.b_.end) - ol.b_.start;
                assert(s < diffs.size() && e < diffs.size());
                
                //LOG(INFO)("REG: (%d %d) - (%d %d)", s, e, diffs[s].mat, diffs[e].mat);
                alldiff.del += diffs[e].del - diffs[s].del;
                alldiff.ins += diffs[e].ins - diffs[s].ins;
                alldiff.mat += diffs[e].mat - diffs[s].mat;
                alldiff.mis += diffs[e].mis - diffs[s].mis;
            } else if (bed[i].start >= ol.b_.end) {
                break;
            }
        }
    }
    LOG(INFO)("Diff: %zd, %zd, %zd, %zd",  alldiff.mat, alldiff.mis, alldiff.ins, alldiff.del);

}

} // namespace fsa
