#include "overlap_refine.hpp"

#include <numeric> 
#include <algorithm>
#include <list>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>

#include "correct/align/tool_aligner.hpp"
#include "./utils/logger.hpp"

namespace fsa {


ArgumentParser OverlapRefine::GetArgumentParser() {
    ArgumentParser ap("fsa_ol_refine", "refine overlaps", "1.0");

    ap.AddPositionOption(ifname_, "ifname", "overlap file name");
    ap.AddPositionOption(ofname_, "ofname", "overlap file name");
    ap.AddNamedOption(itype_, "itype", "input overlap files type");
    ap.AddNamedOption(otype_, "otype", "output overlap files type");
    ap.AddNamedOption(read_fname_, "read_fname", "read file name");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    ap.AddNamedOption(task_, "task", "task performed by program.");
    ap.AddNamedOption(filter0_opts_, "filter0", "overlap filtering options", "");
    ap.AddNamedOption(filter1_opts_, "filter1", "overlap filtering options", "");
    ap.AddNamedOption(aligner_opts_, "aligner", "set an aligner", "");

    return ap;
}


void OverlapRefine::CheckArguments() {
    filter0_.From(filter0_opts_);
    filter1_.From(filter1_opts_);
    
    filter0_opts_ = filter0_.ToString();
    filter1_opts_ = filter1_.ToString();
}

void OverlapRefine::Running() {
    if (task_ == "extend" || task_ == "identity") {
        if (!read_fname_.empty()) rd_store_.Load(read_fname_);
    }

    bool is_same_type = DetectFileType();

    std::unique_ptr<Reader> reader(GetReader(ifname_));
    std::unique_ptr<Writer> writer(GetWriter(ofname_));

    if (reader != nullptr && writer != nullptr) {
        std::mutex mutex_read;
        std::mutex mutex_write;


        auto write_func = [&mutex_write, &writer](std::ostringstream  &filtered) {
            std::lock_guard<std::mutex> lock(mutex_write);
            writer->Write(filtered.str());
            writer->Flush();
            filtered.str("");
        };

        auto work_func = [this, &reader, &mutex_read, write_func, is_same_type](size_t id) {
            std::ostringstream filtered;
            StringPool::TempNameId ni;
            LineInBlock line_in_block(*reader.get(), 10000000, &mutex_read);

            size_t total = 0;
            auto task = GetTask(task_);
            std::string line;
            while (line_in_block.GetLine(line)) {
                fsa::Overlap ol;
                int replen = 0;
                auto r = fromLine(line, ol, ni, replen); 
                if (r > 0) {
                    if (filter0_.Valid(ol, replen)) {
                        (this->*task)(ol, ni);

                        if ( task == &OverlapRefine::TaskIdle || filter1_.Valid(ol)) {
                            if (is_same_type && task == &OverlapRefine::TaskIdle) {
                                filtered << line << "\n";
                            } else {
                                filtered << toLine(ol, ni) << "\n";
                            }
                        }
                    }
                }
                if (filtered.tellp() > 10000000) {
                    total += filtered.tellp();
                    write_func(filtered);
                }
            }
            total += filtered.tellp();
            write_func(filtered);
        };

        MultiThreadRun(thread_size_, work_func);
        
    } else {
        if (reader == nullptr) LOG(ERROR)("Failed to open file: %s", ifname_.c_str());
        if (writer == nullptr) LOG(ERROR)("Failed to open file: %s", ofname_.c_str());
    }

}

 
void OverlapRefine::TaskExtend(Overlap &o, const StringPool::NameId &ni) {
    thread_local auto aligner = ToolAligner::Create(aligner_opts_);


    const auto & query = rd_store_.GetSeq(ni.QueryNameById(o.a_.id));
    const auto & target = rd_store_.GetSeq(ni.QueryNameById(o.b_.id));

    if (o.SameDirect()) {
    
        std::vector<uint8_t> qseq(query.Size(), 0);
        std::vector<uint8_t> tseq(target.Size(), 0);
        for (size_t i=0; i<qseq.size(); ++i) {
            qseq[i] = query[i];
        }
        for (size_t i=0; i<tseq.size(); ++i) {
            tseq[i] = target[i];
        }
        if (o.a_.start > 0 && o.b_.start > 0) {    
            size_t qlen = o.a_.start;
            size_t tlen = o.b_.start;
            size_t qpos = qlen - 1;
            size_t tpos = tlen - 1;
            Alignment al;
            auto r = aligner->Align((const char*)&qseq[0], qlen, (const char*)&tseq[0], tlen, {qpos, qpos}, {tpos, tpos}, al);
            if (r) {
                o.a_.start = al.query_start;
                o.b_.start = al.target_start;
            }
        }

        if (o.a_.len - o.a_.end > 0 && o.b_.len - o.b_.end > 0) {
            size_t qlen = o.a_.len - o.a_.end;
            size_t tlen = o.b_.len - o.b_.end;
            size_t qpos = 0;
            size_t tpos = 0;
            Alignment al;
            auto r = aligner->Align((const char*)&qseq[0]+o.a_.end, qlen, 
                                    (const char*)&tseq[0]+o.b_.end, tlen, {qpos, qpos}, {tpos, tpos}, al);
            if (r) {
                o.a_.end += al.query_end;
                o.b_.end += al.target_end;
            }
        }
    } else {
        std::vector<uint8_t> qseq(query.Size(), 0);
        std::vector<uint8_t> tseq(target.Size(), 0);
        for (size_t i=0; i<qseq.size(); ++i) {
            qseq[i] = 3 - query[qseq.size()-1-i];
            assert(qseq[i] <4 && qseq[i] >= 0);
        }
        for (size_t i=0; i<tseq.size(); ++i) {
            tseq[i] = target[i];
        }
        
        if (o.a_.len - o.a_.end > 0 && o.b_.start > 0) {    
            size_t qlen = o.a_.len - o.a_.end;
            size_t tlen = o.b_.start;
            size_t qpos = qlen - 1;
            size_t tpos = tlen - 1;
            Alignment al;
            auto r = aligner->Align((const char*)&qseq[0], qlen, (const char*)&tseq[0], tlen, {qpos, qpos}, {tpos, tpos}, al);
            if (r) {
                o.a_.end += qlen - al.query_start;
                o.b_.start = al.target_start;
            }
        }

        if (o.a_.start > 0 && o.b_.len - o.b_.end > 0) {
            size_t qlen = o.a_.start;
            size_t tlen = o.b_.len - o.b_.end;
            size_t qpos = 0;
            size_t tpos = 0;
            Alignment al;
            auto r = aligner->Align((const char*)&qseq[0]+(o.a_.len-qlen), qlen, 
                                    (const char*)&tseq[0]+o.b_.end, tlen, {qpos, qpos}, {tpos, tpos}, al);
            if (r) {
                o.a_.start = qlen - al.query_end;
                o.b_.end += al.target_end;
            }
        }
    }
}

void OverlapRefine::TaskIdentity(Overlap &o, const StringPool::NameId &ni) {
    thread_local auto aligner = ToolAligner::Create(aligner_opts_);


    const auto & query = rd_store_.GetSeq(ni.QueryNameById(o.a_.id));
    const auto & target = rd_store_.GetSeq(ni.QueryNameById(o.b_.id));

    if (o.SameDirect()) {
    
        std::vector<uint8_t> qseq(query.Size(), 0);
        std::vector<uint8_t> tseq(target.Size(), 0);
        for (size_t i=0; i<qseq.size(); ++i) {
            qseq[i] = query[i];
        }
        for (size_t i=0; i<tseq.size(); ++i) {
            tseq[i] = target[i];
        }
        size_t qlen = o.a_.end - o.a_.start;
        size_t tlen = o.b_.end - o.b_.start;
        size_t qpos = 0;
        size_t tpos = 0;
        Alignment al;
        auto r = aligner->Align((const char*)&qseq[0]+o.a_.start, qlen, 
                                (const char*)&tseq[0]+o.b_.start, tlen, {qpos, qpos}, {tpos, tpos}, al);
        if (r) {
            printf("%0.02f %s\n", al.Identity(), toLine(o, ni).c_str());
        }
    } else {
        std::vector<uint8_t> qseq(query.Size(), 0);
        std::vector<uint8_t> tseq(target.Size(), 0);
        for (size_t i=0; i<qseq.size(); ++i) {
            qseq[i] = 3 - query[qseq.size()-1-i];
            assert(qseq[i] <4 && qseq[i] >= 0);
        }
        for (size_t i=0; i<tseq.size(); ++i) {
            tseq[i] = target[i];
        }
        
        size_t qlen = o.a_.end - o.a_.start;
        size_t tlen = o.b_.end - o.b_.start;
        size_t qpos = 0;
        size_t tpos = 0;
        Alignment al;
        auto r = aligner->Align((const char*)&qseq[0]+(o.a_.len-o.a_.end), qlen, 
                                (const char*)&tseq[0]+o.b_.start, tlen, {qpos, qpos}, {tpos, tpos}, al);
        if (r) {
            printf("%0.02f %s\n", al.Identity(), toLine(o, ni).c_str());
        }
    }
}

Reader* OverlapRefine::GetReader(const std::string &fname) {
    // ;
    if (fname == "-") {
        return new StdioReader;
    } else {
        GzFileReader* reader = new GzFileReader(fname);
        if (!reader->Valid()) {
            delete reader;
            reader = nullptr;
        }
        return reader;
    }
}

Writer* OverlapRefine::GetWriter(const std::string &fname) {
    if (fname == "-") {
        return new StdioWriter;
    } else {
        GzFileWriter* writer = new GzFileWriter(fname);
        if (!writer->Valid()) {
            delete writer;
            writer = nullptr;
        }
        return writer;
    }
}

bool OverlapRefine::DetectFileType() {
    if (itype_ == "") itype_ = OverlapStore::DetectFileType(ifname_);
    if (otype_ == "") otype_ = OverlapStore::DetectFileType(ifname_);

    std::unordered_map<std::string, int> type_mapper = {{"m4", 0}, {"m4.gz", 0},{"m4a", 1}, {"m4a.gz", 1},{"paf", 2}, {"paf.gz", 2}};

    LOG(INFO)("itype = %s, otype = %s", itype_.c_str(), otype_.c_str());

    if (itype_ == "m4" || itype_ == "m4.gz") {
        fromLine = &OverlapStore::FromM4LineEx;
    } else if (itype_ == "m4a" || itype_ == "m4a.gz") {
        fromLine = &OverlapStore::FromM4aLineEx;
    } else if (itype_ == "paf" || itype_ == "paf.gz") {
        fromLine = &OverlapStore::FromPafLineEx;
    } else {
        LOG(ERROR)("Failed to recognize input overlap files type: %s", itype_.c_str());
    }

    if (otype_ == "m4" || otype_ == "m4.gz") {
        toLine = &OverlapStore::ToM4Line;
    } else if (otype_ == "m4a" || otype_ == "m4a.gz") {
        toLine = &OverlapStore::ToM4aLine;
    } else if (otype_ == "paf" || otype_ == "paf.gz") {
        toLine = &OverlapStore::ToPafLine;
    } else {
        LOG(ERROR)("Failed to recognize output overlap files type: %s", itype_.c_str());
    }

    return type_mapper[itype_] == type_mapper[otype_];
}


void (OverlapRefine::* OverlapRefine::GetTask(const std::string& task))(Overlap& o, const StringPool::NameId &ni) {
    if (task == "extend") {
        return &OverlapRefine::TaskExtend;
    } else if (task == "identity") {
        return &OverlapRefine::TaskIdentity;
    } else {
        return &OverlapRefine::TaskIdle;
    }
}



} // namespace fsa {
    