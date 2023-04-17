#pragma once

#include <array>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <numeric> 

#include "../sequence.hpp"
#include "../overlap.hpp"
#include "../utils/logger.hpp"


#include "phs_options.hpp"

namespace fsa {

class ReadStore;

class LocalPhaser {


public:
    LocalPhaser(const PhsOptions::PhaserOptions& opts, ReadOffset tid, const std::unordered_set<ReadOffset>& qids, 
        const std::unordered_map<ReadOffset, ReadInfo>& rdinfos,
        const ReadStore &read_store, bool correct=true);

    void Run();
    void Run2();

    const std::unordered_map<ReadOffset, std::vector<PhaseItem>>& GetInconsistent() const { return inconsistent_; }
    const std::unordered_map<ReadOffset, std::vector<PhaseItem>>& GetConsistent() const { return consistent_; }

    std::vector<int> MapVariants(const std::vector<std::array<int, 5>>& vars);

    struct CompareResult {
        size_t asize { 0 };
        size_t bsize { 0 };
        size_t common_size { 0 };
        double Coverage() const { return common_size*1.0 / bsize; }
        double Similary() const { return common_size > 0 ? same*1.0 / common_size : 0; }
        int diff { 0 };
        int same { 0 };
    };

    struct Query {
        const ReadInfo* query;
        const std::vector<int> vars;  // -1, 0, 1
        size_t Size() const { 
            return std::accumulate(vars.begin(), vars.end(), 0, [](size_t a, int v) {
                return a += (v != 0 ? 1 : 0);
            });
        }

        size_t PSize() const {
            return std::accumulate(vars.begin(), vars.end(), 0, [](size_t a, int v) {
                return a += (v == 1 ? 1 : 0);
            });
        }

        size_t NSize() const {
            return std::accumulate(vars.begin(), vars.end(), 0, [](size_t a, int v) {
                return a += (v == -1 ? 1 : 0);
            });
        }

        size_t Distance(const Query* b) const {
            std::array<size_t, 3> stat = {0, 0, 0};
            for (size_t i= 0; i < vars.size(); ++i) {
                stat[1 + vars[i]*b->vars[i]] ++;
            }
            return stat[0];
        }

        std::array<size_t, 3> Stat(const Query* b) const {
            return Stat(b->vars);
        }

        std::array<size_t, 3> Stat() const {
            std::array<size_t, 3> stat = {0, 0, 0};
            for (size_t i= 0; i < vars.size(); ++i) {
                stat[1 + vars[i]] ++;
            }
            return stat;
        }

        double Similariry(const Query* b) const {
            auto stat = Stat(b);
            return stat[0] * 1.0 / (stat[0] + stat[2]);
        }

        std::array<size_t, 3> Stat(const std::vector<int> &c) const {
            assert(c.size() == vars.size());
            std::array<size_t, 3> stat = {0, 0, 0};
            for (size_t i = 0; i < vars.size(); ++i) {
                stat[1 + vars[i]*c[i]] ++;
            }
            return stat;
        }

        static double Coverage(const std::array<size_t,3>& stat, size_t size) {
            return (stat[0] + stat[2]) * 1.0 / size;
        }

        static double Similary(const std::array<size_t,3>& stat)  {
            return stat[0] + stat[2] != 0 ? stat[2] * 1.0 / (stat[0] + stat[2]) : 0.0;
        }

    };



    bool IsClosedToCentroid(const std::array<size_t,3>& stat, size_t qsize) const;

    struct Group {
        struct Loc {
            Loc(int a = 0, int c = 0) : accu(a), cov(c) {}
            int accu { 0 };
            int cov { 0 };
            bool Notable(double support_rate) const {
                double msr = support_rate*2 - 1;
                return cov > 0 && std::abs(accu) >= cov * msr; 
            }

            bool Valid(const PhsOptions::PhaserOptions &opts) const {
                return std::abs(accu) >= opts.min_support_count && Notable(opts.min_support_rate);
            }
            bool Valid(double support_rate, int support) const {
                return std::abs(accu) >= support && Notable(support_rate);
            }

            int GetValidType(const PhsOptions::PhaserOptions &opts) const  {
                return Valid(opts) ? GetType() : 0;
            }

            int GetValidType(double support_rate, int support) const  {
                return Valid(support_rate, support) ? GetType() : 0;
            }

            int GetNotableType(double support_rate) const  {
                return Notable(support_rate) ? GetType() : 0;
            }

            int GetType() const {
                return accu < 0 ? -1 : (accu > 0 ? 1 : 0);
            }

            void Merge(const Loc& b) {
                accu += b.accu;
                cov  += b.cov;
            }
            void Remove(const Loc& b) {
                accu -= b.accu;
                cov  -= b.cov;
            }
        };

        struct Centroid {
            Centroid(const std::vector<Loc>& ls) : locs(ls), range({0, ls.size()}) {}

            double Distance(const Query& q) const;

            void SetRange2(const Group& g, const Centroid& alt, const Group& alt_group, const PhsOptions::PhaserOptions& opts);
            std::vector<std::array<size_t,3>> SplitRange(const std::vector<int> &cs, const Group& group, const std::vector<int>& alt_cs, const Group& alt_group, double rate, int count);
            std::vector<std::array<size_t,3>> CombineRange(const std::vector<std::array<size_t,3>>& ranges, const std::vector<int> &cs, const Group& group, const std::vector<int>& alt_cs, const Group& alt_group, const PhsOptions::PhaserOptions& opts);
            std::vector<std::array<size_t,2>> GetRangeLink(const std::vector<std::array<size_t,3>>& ranges, const std::vector<int> &cs, const Group& group, double rate, int count);
            CompareResult Compare(const Query &q) const;
            std::vector<int> Values(double rate, int count) const;
            std::vector<int> Values(double rate, int count, const Centroid& alt) const;
            std::vector<Loc> locs;
            std::array<size_t, 2> range;
        };

        Group() {};
        Group(const Query* q);
        Group(const std::unordered_set<const Query*>& qs);
        Group(const std::vector<Query>& qs);



        size_t Size() const { return queries.size(); }
        std::array<int, 3> StatDiffs() const ;
        std::array<int, 3> StatDiffs(double support_rate) const ;
        std::array<int, 3> StatDiffs(double support_rate, int support) const ;
        void Merge(const Group& b);
        void Merge(const Query* q);
        void Remove(const Query* q);

        std::vector<Centroid> SelectCentroids() const;
        Centroid GetCentroid() const { return Centroid(vars); };
        std::array<double,2> CalcGroupSSE(const PhsOptions::PhaserOptions& opts) const;

        double DistanceFromTarget() const {
            auto stat = StatDiffs();
            return stat[0] + stat[2] == 0 ? 0.0 : stat[0] * 1.0 / (stat[0] + stat[2]);
        }
        std::array<int, 3> Distance(const Group& b, double support_rate) const;
        std::array<int, 3> Distance(const Group& b, double support_rate, int support) const;
        int MaxDistance(const Group &b, double support_rate, int support) const;

    
        std::vector<std::array<int,2>> StatVars(const std::unordered_map<int, int>& mapper) const;

        bool IsConsistentWithGroup0(const Group& b, const PhsOptions::PhaserOptions &opts) const;
        bool IsConsistentWithGroup0(const Query& q, const PhsOptions::PhaserOptions &opts) const;
        bool IsInconsistentWithGroup0(const Group& b, const PhsOptions::PhaserOptions &opts, const std::vector<int>& coverage) const;
        bool IsInconsistentWithGroup0(const Query& q, const Group &b, const PhsOptions::PhaserOptions &opts) const;
        bool IsConsistentWith(const Group& b, const PhsOptions::PhaserOptions &opts) const;

        size_t FirstVaildVar(double support_rate, int support) const {
            for (size_t i = 0; i < vars.size(); ++i) {
                if (vars[i].Valid(support_rate, support)) {
                    return i;
                }
            }
            return vars.size();
        }
        
        bool Connected() const ;
        std::vector<Group> Split() const;

        std::unordered_set<const Query*> queries;
        std::vector<Loc> vars;  // ((-n, 0, n), cov
    };

    std::vector<Group> Divide(const Group& g);
    std::vector<Group> Combine(const std::vector<Group>& groups);
    std::vector<Group> Combine2(const std::vector<Group>& groups);
    std::vector<Group> Combine3(const std::vector<Group>& groups);
    void FindConsistent(const Group& g0);
    void FindInconsistent(const Group& g0, const Group& gx);

    std::vector<Group> BiKMeans(const Group& g);
    std::vector<Group> SplitGroups(const Group& g, const std::vector<Group::Centroid> &centroids);
    size_t AdjustGroups(std::vector<LocalPhaser::Group> & groups);
    size_t DistributeQuery(const std::vector<Group::Centroid>& centroids, const Query* query) const;

    bool IsConsistent(const Group &g0, const Query& query) const;
    void CorrectTarget(const std::vector<Group>& groups, const std::vector<size_t>& inconsist);
    void CorrectTarget(const std::vector<Group>& groups);
    void ClearTarget();
    PhaseItem GetMapItem(const ReadInfo& target, const ReadInfo& query) const;
    void PrintGroups(const std::string& msg, const std::vector<Group>& groups) const;
    void PrintSet(const std::string &name, std::unordered_set<size_t> &set);


    // 通过Merge
    std::vector<Group> MergeGroups(std::vector<Group> &groups, double similary, double coverage, int diff);
protected:
    PhsOptions::PhaserOptions opts_;

    const ReadStore &rd_store_;
    bool correct_ { true };
    
    const ReadInfo* target_ { nullptr };
    std::vector<Query> queries_;
    const Query* target_in_queires_;

    std::unordered_map<int, int> mapper; 
    std::unordered_map<ReadOffset, std::vector<PhaseItem>> inconsistent_;
    std::unordered_map<ReadOffset, std::vector<PhaseItem>> consistent_;
    std::vector<int> coverages_;
};

} // namespace fsa {

