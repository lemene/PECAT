#pragma once

#include "utility.hpp"

namespace fsa {

class StringPool;


struct HicReadMapping {
    uint32_t ctg;
    uint32_t offset;
    //double quality;
    std::vector<SnpAllele> alleles;
};

struct HicReadInfo {
    std::array<std::vector<HicReadMapping>, 2> hic;
};

class HicReadInfos {
public:
    HicReadInfos(StringPool &sp) : string_pool_(sp) {}

    void Build(const std::string& fn_hic1, const std::string& fn_paf1, 
               const std::string& fn_hic2, const std::string& fn_paf2, 
               const std::string& fn_vars, size_t thread_size);

    void Save(const std::string& fname) const;
    void Load(const std::string& fname);
    const std::unordered_map<size_t, HicReadInfo>& GetInfos() const { return infos_; }
protected:
    void BuildOne(const std::string& fn_hic, const std::string& fn_paf, size_t i, size_t sp_offset, const class SnpStore &vars, size_t thread_size);
    
    void LoadLine(const std::string &line);
    HicReadMapping GetMappingInfo(const Overlap& ol, const class ReadStore& hic_store, const class SnpStore &vars);
protected:
    StringPool& string_pool_;
    std::unordered_map<size_t, HicReadInfo> infos_;
};
    
} // namespace fsa 
