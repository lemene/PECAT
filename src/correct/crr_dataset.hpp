#pragma once

#include <atomic>
#include <memory>

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "corrector.hpp"
#include "overlap/mapping.hpp"

namespace fsa {

class CrrOptions;
class Dispatcher;

struct CrrDataset {
public:
    CrrDataset(CrrOptions &opt) : opts_(opt) {}

    void Load();
    const std::string& QueryStringById(Seq::Id id) const { return string_pool_.QueryStringById(id); }    

    
    std::unique_ptr<Dispatcher> GetDispatcher();
    std::vector<std::vector<Seq::Id>> GroupReadIds() const;

    class OlGroup {
    public:
        OlGroup(Seq::Id i) : id(i) {}
        OlGroup(const OverlapGrouper::Group& gp) : id (gp.id), ols(gp.ols), index(gp.index) {}

        bool Empty() const { return ols.size() == 0; }
        size_t Size() const { return index.size(); }
        size_t Size(size_t i) const { return index[i][1] - index[i][0]; }
        const Overlap* Get(size_t i, size_t j) { return ols[index[i][0]+j]; }

        void Sort(double opt_ohwt);
        std::vector<double> GetWeight(double opt_ohwt);

        void BuildIndex();

        Seq::Id id;
        std::vector<const Overlap*> ols;
        std::vector<std::array<size_t, 2>> index;
        std::shared_ptr<std::vector<Overlap>> from_mapping;
    };
    OlGroup Get(int id) const { return OlGroup(grouper_.Get(id)); }
    OlGroup Get2(int id) const;
protected:
    void LoadReadIds();
    void LoadOverlaps();
    void LoadReads();
    void EstimateParameters();
    
public:
    CrrOptions& opts_;
    
    StringPool string_pool_;
    ReadStore read_store_ {string_pool_};
    OverlapStore ol_store_{string_pool_ };
    
    OverlapGrouper grouper_ { ol_store_ };
    OverlapStore rd_2_ref_ {string_pool_ };
    Mapping mapping_ { rd_2_ref_ };

    
    std::vector<Seq::Id> read_ids_;
};


struct Dispatcher { 
    virtual std::vector<Seq::Id> Get() = 0;
};
struct SimpleDispatcher : public Dispatcher {
    SimpleDispatcher(const CrrDataset& d, size_t n=1)
        : read_ids_(d.read_ids_), N(n) {}
    virtual std::vector<Seq::Id> Get() {
        std::vector<Seq::Id> ids;
        const int N = 1;
        auto curr = index.fetch_add(N);
        for (size_t i = 0; i < N; ++i) {
            if (curr + i < read_ids_.size()) {
                ids.push_back(read_ids_[curr+i]);
            }
        }
        return ids;
    }
    
    size_t N { 1 };
    const std::vector<Seq::Id>& read_ids_;
    std::atomic<size_t> index { 0 };
};


struct GroupDispatcher : public Dispatcher {
    GroupDispatcher(const CrrDataset& d)
        : clu_ids_(d.GroupReadIds()) { }

    virtual std::vector<Seq::Id> Get() {
        auto curr = index.fetch_add(1);
        if (curr < clu_ids_.size()) {
            return clu_ids_[curr];
            
        } else {
            return std::vector<Seq::Id>();
        }
    }


    std::atomic<size_t> index { 0 };
    std::vector<std::vector<Seq::Id>> clu_ids_;

};

} // namespace fsa
