#include "overlap_tools.hpp"

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "utils/misc.hpp"
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
                    if (consistent.Contain(o, threshold_) || (!phased_info.IsRemoved(o.a_.id) && !phased_info.IsRemoved(o.b_.id) &&
                        !phased_info.Contain(o, threshold_))) {
                            
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
    PrjFile::ContigTiles tiles(sp);

    tiles.Load(tiles_);

    ReadVariants rvs(sp);
    rvs.Load(readinfos_);

    std::unordered_map<int, std::unordered_map<uint64_t, int>> ctginfos;

    for (auto &ctg : tiles.GetTiles()) {
        auto &cinfo = ctginfos[ctg.first];

        for (auto r : ctg.second) {
            const auto& vs = rvs.GetVariants(r);
            
        }
    }
}

} // namespace fsa
