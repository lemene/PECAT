#pragma once

#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "utils/program.hpp"
#include "utils/string_pool.hpp"
#include "utility.hpp"
#include "file_io.hpp"

namespace fsa {

inline std::string Format(const std::string& pattern, int d) {
            std::string result = pattern;
            std::string s = std::to_string(d);
            result.replace(pattern.find("{}"), 2, s);
            return result;
        };

class Program_SplitOverlaps : public Program {
public:
    Program_SplitOverlaps() {
        name_ = "split_ols";
        desc_ = "split overlap files for correction";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "overlap file name");
        ap.AddPositionOption(ofname_, "ofname", "output file name pattern '{}' is replaced by seqnum");
        ap.AddPositionOption(rd_fname_, "reads", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(rdfname0_, "rdfname0", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(rdfname1_, "rdfname1", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(block_size_, "block_size", "size of split files");
        ap.AddNamedOption(thread_size_, "thread_size", "");

        return ap;
    }

    virtual void Running();
        using Key = int;
    struct Group {
        
        Group(int id, std::string fname) {
            seqid = id;
            writer.reset(new GzFileWriter(fname));
        }
        ~Group() {
            Flush();
        }
        size_t Size() const { return size; }

        bool IsAssigned(const Key &n) const {
            return reads.find(n) != reads.end();
        }

        void Add(const std::string& line, const Key &b) {
            writer->Write(line);
            reads2.insert(b);
        }
        void Add( const Key &b) {
            reads2.insert(b);
        }

        void Assign(const Key &a, int len) {
            reads.insert(a);
            reads2.insert(a);
            size += len;
        }

        void AddLine(const std::string &line) {
            if (line_buff_index >= line_buff.size()) {
                Flush();
            } 
            
            line_buff[line_buff_index] = line;
            line_buff_index++;
        }
    
        void Flush() {
            for (size_t i=0; i<line_buff_index; ++i) {
                writer->Write(line_buff[i]);
            }
            line_buff_index = 0;
        }
        int seqid;
        std::unordered_set<Key> reads;
        size_t size { 0 };
        std::unordered_set<Key> reads2;
        std::shared_ptr<GzFileWriter> writer;
        std::vector<std::string> line_buff{10000};
        size_t line_buff_index {0};
    };
    
    struct WorkArea {
        WorkArea(int block_size) : assigned(block_size) {

        }

        void Assign(size_t i, int aid, int alen, int adist, int bid, int blen, int bdist) {
            assigned[i][0] = aid;
            assigned[i][1] = alen;
            assigned[i][2] = adist;
            assigned[i][3] = bid;
            assigned[i][4] = blen;
            assigned[i][5] = bdist;
        }
        std::vector<std::array<int, 6>> assigned;
    };

    struct Infos {
        Infos(const std::string &ofn): ofname(ofn) {
            groups.push_back(Group(0, Format(ofname, 0)));
        }
        struct Item {
            std::vector<std::string> lines;
            int count { 0 };
        };

   


        void CheckFull(int block_size) {
            
            if (groups.back().Size() > (size_t)block_size) {
                LOG(INFO)("size: %zd %zd", groups.back().Size(), groups.size());
                groups.push_back(Group(groups.size(), Format(ofname, groups.size())));
            }
        }

        int IsAssigned(const Key &n) const {
            auto iter = assigned.find(n);
            return iter != assigned.end() ? iter->second : -1;
        }

        void Assign(const Key &n, int g) { assigned[n] = g; }

        void Merge(size_t size, WorkArea& wa, const std::vector<std::string> &lines) {
            
            for (size_t i = 0; i < size; ++i) {
                auto &assd = wa.assigned[i];
                if (assd[2] > 0) {
                    groups[assd[2]].Add(lines[i], assd[3]);
                } else {

                }
            }
        }

        std::unordered_map<Key, int> assigned;
        std::vector<Group> groups;
        std::string ofname;
    };


    void GetNameLenFromPafLine(std::string& aname, int alen, std::string &bname, int &blen);
    void SaveReadName(const std::vector<Group> &group);
    void SaveReadName(const std::unordered_set<Key>& read, const std::string &fname);


    std::string ifname_ ;
    std::string ofname_ ;
    std::string rd_fname_;
    std::string rdfname0_ {""};
    std::string rdfname1_ {""};
    long long int block_size_ { 4000000000L };

    StringPool string_pool_;
    int thread_size_ { 1 };

};


class Program_SplitName : public Program {
public:
    Program_SplitName() {
        name_ = "split_name";
        desc_ = "split overlap files for correction";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(fn_rds_, "reads", "output file name pattern '{}' is replaced by seqnum");
        ap.AddPositionOption(fn_part_rdnames_, "partitoned_readnames", "output file name pattern '{}' is replaced by seqnum");
        
        ap.AddNamedOption(fn_part_ols_, "partitioned_overlaps", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(fn_ols_, "overlaps", "overlap file name");
        ap.AddNamedOption(fn_rdnames_, "readname", "which reads are grouped");
        ap.AddNamedOption(block_size_, "block_size", "size of split files");
        ap.AddNamedOption(base_size_, "base_size", "");
        ap.AddNamedOption(thread_size_, "thread_size", "");
        return ap;
    }

    virtual void Running();

    void LoadReadnames();
    void GroupReadsRandomly();
    void GroupReadsByOverlaps();
    void Group(const std::vector<std::vector<int>>& groups, const class ReadStore &rd_store);
    void SaveReadnames();
    void SaveOverlaps(const std::string &opattern);
    std::string fn_ols_ ;
    std::string fn_part_ols_ ;
    std::string fn_rds_ ;
    std::string fn_part_rdnames_ ;
    std::string fn_rdnames_ { "" };
    long long int block_size_ { 4000000000L };
    long long int base_size_ { 0 };

    StringPool string_pool_;
    std::vector<int> lengths_;

    int thread_size_ { 1 };

    std::unordered_map<Seq::Id, int> read_ids_index_;
    std::vector<Seq::Id> read_ids_;
    std::vector<int> groups_;
};

class Program_SplitOverlaps2 : public Program {
public:
    Program_SplitOverlaps2() {
        name_ = "split_ols2";
        desc_ = "split overlap files for correction";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "overlap file name");
        ap.AddPositionOption(ofname_, "ofname", "output file name pattern '{}' is replaced by seqnum");
        ap.AddPositionOption(rd_fname_, "reads", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(rdfname0_, "rdfname0", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(rdfname1_, "rdfname1", "output file name pattern '{}' is replaced by seqnum");
        ap.AddNamedOption(sub_size_, "sub_size", "");
        ap.AddNamedOption(thread_size_, "thread_size", "");
        return ap;
    }

    virtual void Running();

    
    struct Group {
        
        Group(int id, std::string fname) {
            seqid = id;
            writer.reset(new GzFileWriter(fname));
        }
        ~Group() {
            Flush();
        }
        
        void Add(const std::string& line, const int &b) {
            writer->Write(line);
            reads2.insert(b);
        }
        void Add( const int &b) {
            reads2.insert(b);
        }


        void AddLine(const std::string &line) {
            if (line_buff_index >= line_buff.size()) {
                Flush();
            } 
            
            line_buff[line_buff_index] = line;
            line_buff_index++;
        }
    
        void Flush() {
            for (size_t i=0; i<line_buff_index; ++i) {
                writer->Write(line_buff[i]);
            }
            line_buff_index = 0;
        }
        int seqid;
        std::unordered_set<int> reads2;

        std::shared_ptr<GzFileWriter> writer;
        std::vector<std::string> line_buff{10000};
        size_t line_buff_index {0};
    };
    
    struct WorkArea {
        WorkArea(int block_size) : assigned(block_size) {

        }

        void Set(size_t i, int aid, int alen, int adist, int bid, int blen, int bdist) {
            assigned[i][0] = aid;
            assigned[i][1] = alen;
            assigned[i][2] = adist;
            assigned[i][3] = bid;
            assigned[i][4] = blen;
            assigned[i][5] = bdist;
        }
        std::vector<std::array<int, 6>> assigned;
    };

    struct Infos {
        Infos(const std::string &ofn, size_t subsize, size_t readsize): ofname(ofn), assigned(readsize, -1) {
            for (size_t i = 0; i < subsize; ++i) {
                groups.push_back(Group(i, Format(ofname, i)));
            }
        }

        int GetGroup(int i) { return assigned[i]; }
        void SetGroup(int i, int g) { assigned[i] = g; groups[g].Add(i); }

        void Merge(size_t size, WorkArea& wa, const std::vector<std::string> &lines) {
            for (size_t i = 0; i < size; ++i) {
                auto &assd = wa.assigned[i];
                //assert(assd[2] >= 0 && assd[5] >= 0);
                if (assd[2] >= 0) {
                    groups[assd[2]].Add(lines[i], assd[3]);
                }

                if (assd[5] >= 0) {
                    if (assd[2] != assd[5]) { // 包括assd[2] < 0 的情况
                        groups[assd[5]].Add(lines[i], assd[0]);
                    } else {
                        groups[assd[5]].Add(assd[0]);
                    }

                }

            }
        }

        std::string ofname;
        std::vector<int> assigned;
        std::vector<Group> groups;
    };


    void LoadSubNames(const std::string& fname, size_t subsize, Infos &infos);
    void SaveReadName(const std::vector<Group> &group);
    void SaveReadName(const std::unordered_set<int>& reads, const std::string &fname);

    std::string ifname_ ;
    std::string ofname_ ;
    std::string rd_fname_;
    std::string rdfname0_ {""};
    std::string rdfname1_ {""};
    int sub_size_;

    StringPool string_pool_;
    int thread_size_ { 1 };
};

class Program_Test : public Program {
public:
    Program_Test() {
        name_ = "test";
        desc_ = "test some fucntion";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(ifname_, "ifname", "input file");
        ap.AddPositionOption(ofname_, "ofname", "output file");
        ap.AddNamedOption(thread_size_, "thread_size", "");
        return ap;
    }

    virtual void Running();

    
    std::string ifname_ ;
    std::string ofname_ ;
    int thread_size_ { 1 };
};

class Program_Hic : public Program {
public:
    Program_Hic() {
        name_ = "snp_in_hic";
        desc_ = "identify snps in hic reads";
    }
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        ap.AddPositionOption(fn_hic1_, "hic1", "hic1 reads");
        ap.AddPositionOption(fn_paf1_, "hic1_2_ctg", "alignment(PAF format) between hic1 reads and contigs");
        ap.AddPositionOption(fn_hic2_, "hic2", "hic2 reads");
        ap.AddPositionOption(fn_paf2_, "hic2_2_ctg", "alignment(PAF format) between hic2 reads and contigs");
        ap.AddPositionOption(fn_vars_, "variants", "SNP information in contigs");
        ap.AddPositionOption(fn_snp_in_hic_, "snp_in_hic", "result, SNP information in hic reads");

        return ap;
    }
    virtual void Running();
protected:
    std::string fn_hic1_;
    std::string fn_paf1_;
    std::string fn_hic2_;
    std::string fn_paf2_;
    std::string fn_vars_;
    std::string fn_snp_in_hic_;
};

class MiscTools : public MultiProgram {
public:
    MiscTools() {
        Add(new Program_SplitOverlaps());
        Add(new Program_SplitOverlaps2());
        Add(new Program_SplitName());
        Add(new Program_Test());
        Add(new Program_Hic());
    }
};

} // namespace fsa

