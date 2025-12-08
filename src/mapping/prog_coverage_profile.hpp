#pragma once

#include "utils/program.hpp"
#include "overlap_store.hpp"
#include "read_store.hpp"
#include "overlap/mapping.hpp"
#include "align/match_info.hpp"
#include "polish/multi_coverage.hpp"

namespace fsa {
class Program_CoverageProfile : public Program {
public:
    Program_CoverageProfile() {
        name_ = "covprof";
        desc_ = "generate coverage profile of contigs or references from mapping files";
    }

    class Options {
    public: 
        void SetArguments(ArgumentParser &ap) {
            
            ap.AddPositionOption(mapping_fname_, "mapping", "mapping file of reads to reference/contigs");
            ap.AddPositionOption(ctg_fname_, "contig", "contig file");
            ap.AddPositionOption(read_fname_, "reads", "raw read file");

            ap.AddNamedOption(min_identity_, "min_identity", "");
            ap.AddNamedOption(min_local_identity_, "min_local_identity", "");
            
            ap.AddNamedOption(thread_size, "thread_size", "thread size");
            ap.AddNamedOption(min_coverage, "min_coverage", "");
        
            ap.AddNamedOption(output_directory, "output_directory", "The directory for temporary files");
        }
    
            std::string OutputPath(const std::string &fname) { return output_directory+"/"+fname; }
    
            std::string output_directory {"."};
    
            int min_coverage_ { 4 };
            double min_identity_ { 60 };
            double min_local_identity_ { 50 };
            int local_window_size_ { 1000 };
            bool check_local_identity_ { false };
    
            int coverage_ { 50 };
            int window_size_ { 50000 };
            int overlap_size_ { 500 };
    
            std::string ctg_fname_;
            std::string mapping_fname_;
            std::string read_fname_;
    
            int thread_size { 4 };
            int min_coverage { 4 };
            bool debug { false };
            double secondary_to_primary_ratio { 0.80 };
        };
    class Dataset {
    public:
        Dataset(class Options &opt) : opts_(opt) {}
        void Load();
        std::vector<Seq::Id> LoadContigIds(const ReadStore& seq_store);
        void LoadOverlaps(const std::string &fname);

            
        const std::string& QueryStringById(Seq::Id id) const { return string_pool_.QueryStringById(id); }    
        const StringPool& GetStringPool() const { return string_pool_; }
        size_t CountReadMap(Seq::Id id) const ;
        size_t MaxReadLength() const { return max_read_length_; }
    
        void Stat();
        /** reads in  */
        std::vector<Bed> CollectBedFromBam(const std::vector<Seq::Id>& read_ids);
        
        void SelectBestMapping();
            Options& opts_;

            StringPool string_pool_;
            ReadStore seq_store_ {string_pool_};
            std::array<size_t,2> rd_ids_;
            
            OverlapStore rd_2_ctg_ {string_pool_ };
            OverlapGrouper grouper_ { rd_2_ctg_ };
            
            
            std::vector<Seq::Id> ctg_ids_;
            double overlap_quality_median_ {0.0};
            double overlap_quality_mad_ {0.0};
            size_t max_read_length_ {0};
            size_t ave_read_length_ {0};
    
        };

    class Worker {
  
    public:
        Worker(Seq::Id tid, const Dataset& ds) 
         : tid_(tid), dataset_(ds), multi_cov_(ds.seq_store_.GetSeq(tid), ds.overlap_quality_median_, ds.overlap_quality_mad_) {
        }
    
        void ComputeCoverage(size_t thread_size=1);
    
    
        /** Save infomations */
        void DumpMultiCoverage(std::ofstream &of);
    
    protected:
        Seq::Id tid_;
        const Dataset& dataset_;
    
        MultiCoverage multi_cov_;
    };
    
    virtual ArgumentParser GetArgumentParser() {
        ArgumentParser ap(Name(), Description(), "");
        opts_.SetArguments(ap);
        return ap;
    }
    virtual void Running() ;
protected:
    void AnalyzeContigs();

protected:
    Options opts_;
    Dataset dataset_ {opts_};

};
} // namespace fsa