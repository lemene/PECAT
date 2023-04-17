#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <deque>
#include <numeric>
#include <functional>

#include "../sequence.hpp"
#include "../utility.hpp"
#include "../overlap.hpp"
#include "asm_dataset.hpp"
#include "string_node.hpp"
#include "graph/cross_node.hpp"
#include "string_edge.hpp"

namespace fsa {
class ReadVariants;

class SgGraph {
public:
    SgGraph(class AsmDataset &asmdata) : asmdata_(asmdata) {}
    virtual ~SgGraph();
    
    AsmDataset &GetAsmData() { return asmdata_; }
    const class AsmOptions& Options() const { return asmdata_.opts_; }
    const class ReadStore& GetReadStore() const { return asmdata_.rd_store_; }
    
    template<typename F>
    std::vector<SgNode*> CollectNodes(F f) {
        std::vector<SgNode*> ns;
        for (auto& i : org_nodes_) {
            if (f(i.second)) ns.push_back(i.second);
        }
        return ns;
    }

    template<typename F> 
    std::vector<SgEdge*> CollectEdges(F f) {
        std::vector<SgEdge*> es;
        for (auto& i : org_edges_) {
            if (f(i.second)) es.push_back(i.second);
        }
        return es;
    }

    template<typename T, typename F, typename C>
    std::vector<T> CollectObjectsByEdge(F f, C check) {
        std::vector<T> es;
        for (auto& i : org_edges_) {
            T t = f(i.second);
            if (check(t)) {
                es.push_back(t);
            }
        }
        return es;
    }

    
    SgNode* ReverseNode(SgNode* n) {
        return org_nodes_[SgNode::ID::Reverse(n->Id())];
    }
    const SgNode* ReverseNode(const SgNode* n) {
        return org_nodes_[SgNode::ID::Reverse(n->Id())];
    }


    SgEdge* ReverseEdge(SgEdge* e) {
        return org_edges_[SgEdge::ID::Reverse(e->Id())];
    }
    
    std::vector<SgEdge*> ReversePath(const std::vector<SgEdge*> &path) {
        std::vector<SgEdge*> vpath(path.size(), nullptr);
        std::transform(path.rbegin(), path.rend(), vpath.begin(), [&](SgEdge* e) { return ReverseEdge(e); });
        return vpath;
    }

    std::vector<SgEdge*> GetLinearPath(SgEdge* start, size_t max_length, size_t max_nodesize);
    size_t PathLength(const std::vector<SgEdge*> &path);



    SgNode* QueryNode(const SgNode::ID &id) {
        auto iter = org_nodes_.find(id);
        return iter != org_nodes_.end() ? iter->second : nullptr;
    }

    SgEdge* QueryEdge(const SgEdge::ID &id) {
        auto iter = org_edges_.find(id);
        return iter != org_edges_.end() ? iter->second : nullptr;
    }

    void InsertEdge(SgEdge* e) {
        assert(QueryEdge(e->Id()) == nullptr);

        org_edges_[e->Id()] = e;
        assert(e->InNode() != nullptr);
        assert(e->OutNode() != nullptr);
        e->InNode()->AddOutEdge(e);
        e->OutNode()->AddInEdge(e);

        e->ReduceSubEdges();
    }

    void InsertNode(SgNode* n) {
        org_nodes_[n->Id()] = n;
    }
    
    void Simplify(const std::string &strategy, const std::string &reducer="");

    std::list<SgNode*> GetEgoNodes(SgNode* n, int max_depth);
    std::list<SgNode*> GetEgoNodes(SgNode* n, int max_depth, int max_len) ;
    std::list<SgNode*> GetEgoNodes(SgNode* n, int max_depth, int max_len, int max_nodesize) ;

    std::list<SgNode*> GetNeighborNodes(SgNode* n, int max_depth);

    static std::vector<std::vector<SgEdge*>> MaximumFlow(const SgNode* src, const SgNode *dst, const std::unordered_set<const SgEdge*> &edges);

    template<typename C>
    static size_t PathNodeSize(const C& path)  {
        return std::accumulate(path.begin(), path.end(), 0, [](size_t a, const typename C::value_type& b) {
            return a + b->NodeSize();
        });   
    }
    AsmDataset &asmdata_;
    
    std::unordered_map<SgNode::ID, SgNode*, SgNode::ID::Hash> org_nodes_;
    std::unordered_map<SgEdge::ID, SgEdge*, SgEdge::ID::Hash> org_edges_;

    std::unordered_map<std::string,  std::shared_ptr<class Simplifier>> simplifiers_;

};

class StringGraph : public SgGraph {
public:
    StringGraph(class AsmDataset &asmdata);
    virtual ~StringGraph();
public:
    std::string IdToString(BaseNode::ID id) const;

    static BaseNode::ID ReverseNode(BaseNode::ID id) {
        return SgNodeID::Reverse(id);
    }
    static BaseEdge::ID ReverseEdge(BaseEdge::ID id) {
        return BaseEdge::ID::Reverse(id);
    }
    BaseNode* ReverseNode(BaseNode* n) {
        return nodes_[ReverseNode(n->Id())];
    }

    BaseEdge* ReverseEdge(BaseEdge* e) {
        return edges_[ReverseEdge(e->Id())];
    }
    
    std::list<BaseEdge*> Reverse(const std::list<BaseEdge*>& path);

    BaseNode* GetNode(BaseNode::ID id) {
        auto i = nodes_.find(id);
        return i != nodes_.end() ? i->second : nullptr;
    }

    BaseEdge* GetEdge(BaseNode::ID iid, BaseNode::ID oid) {
        auto i = edges_.find(BaseEdge::ID(0, iid, oid));
        return i != edges_.end() ? i->second : nullptr;
    }

    void AddOverlap(const Overlap* overlap);
    void AddOverlaps(const class OverlapStore &ol_store);
    void AddOverlaps(const class OverlapStore &ol_store, int min_length, int min_aligned_lenght, float min_identity);
    bool FilterOverlap(const Overlap &ovlp, const std::unordered_set<BaseNode::ID, BaseNode::ID::Hash> &contained, int min_length, int min_aligned_length, float min_identity);    
    void AddEdge(const Overlap* ol, int in_node, int out_node, int read);


    std::unordered_set<BaseNode*> BfsNodes(BaseNode* n, BaseNode* exclude=nullptr, int depth=5);
    
    template<typename C>
    void ReduceEdges(const C& c,  BaseEdge::ReduceType type) {
        for (auto& e : c) {
            if (!e->IsReduce()) {
                e->Reduce(type);
                ReverseEdge(e)->Reduce(type);
            }
        }
    }

    template<typename C>
    void ReactiveEdges(const C& c) {
        for (auto& e : c) {
            if (e->IsReduce()) {
                e->Reactivate();
                ReverseEdge(e)->Reactivate();
            }
        }
    }


    struct LinearPath {
        std::vector<BaseEdge*> path;
        int length;
    };

    template<typename F>
    std::vector<BaseNode*> CollectNodes(F f) {
        std::vector<BaseNode*> ns;
        for (auto& i : nodes_) {
            if (f(i.second)) ns.push_back(i.second);
        }
        return ns;
    }

    template<typename F> 
    std::vector<BaseEdge*> CollectEdges(F f) {
        std::vector<BaseEdge*> es;
        for (auto& i : edges_) {
            if (f(i.second)) es.push_back(i.second);
        }
        return es;
    }

    LinearPath FindBridgePath(BaseEdge* start, int max_length, int max_nodesize);

 //   void PhaseCross();
    bool IsDiffType(const BaseNode* n0, const BaseNode* n1);

    void IdentifySimplePaths();
    std::list<BaseEdge*> ExtendSimplePath(BaseEdge*n, std::unordered_set<BaseEdge*> &visited);

        
    std::list<std::list<BaseEdge*>>& GetPaths() {  return paths_; }

    std::vector<BaseEdge*> ShortestPath(const BaseNode* src, const BaseNode *dst, std::unordered_set<BaseEdge*> edges, int(*score)(BaseEdge*) = [](BaseEdge*) {return 1; });
    std::vector<BaseEdge*> GetPath(BaseNode* src, BaseNode *dst, size_t max_depth);

    // To check whether path corresponds to a reversed path.
    bool Assert_PathDual(const std::list<std::list<BaseEdge*>> paths);
    void SaveEdges(const std::string &fname);


    double GetOverlapQuality(const Overlap& ol);
    void ReduceOtherEdges(const std::unordered_set<BaseEdge*> reserved, BaseEdge::ReduceType type);
protected:
    std::list<std::list<BaseEdge*>> paths_;
    
    std::unordered_map<BaseNode::ID, BaseNode*, BaseNode::ID::Hash>& nodes_;
    std::unordered_map<BaseEdge::ID, BaseEdge*, BaseEdge::ID::Hash>& edges_;


friend class PhaseCrossSimplifier;

public:
    double desired_edge_quality { 0.0 };
    double actual_min_identity { 0.0 };
};



class PathGraph : public SgGraph {
public:
    PathGraph(class AsmDataset &asmdata);
    virtual ~PathGraph();

    void BuildFrom(StringGraph& sg);
    void AddEdge(std::list<BaseEdge*> &path);

    PathNode* ReverseNode(PathNode* n) {
        return nodes_[PathNode::ReverseId(n->Id())];
    }

    PathEdge* ReverseEdge(PathEdge* e) {
        return edges_[PathEdge::ReverseId(e->Id())];
    }

    std::vector<PathEdge*> ReversePath(const std::vector<PathEdge*> &path) {
        std::vector<PathEdge*> vpath(path.size(), nullptr);
        std::transform(path.rbegin(), path.rend(), vpath.begin(), [&](PathEdge* e) { return ReverseEdge(e); });
        return vpath;
    }

    //void IdentifyPathSpur();
    void IdentifyPathSpur2();
    std::vector<PathEdge*> ShortestPath(const PathNode* src, const PathNode *dst, 
        std::unordered_set<PathNode*> candnodes, int(*score)(PathEdge*) = [](PathEdge*) {return 1; });

    void FindBubbles(bool check);

    struct LinearPath {
        LinearPath() {}
        LinearPath(const std::vector<PathEdge*>& p, int l = 0, int n = 0, int s = 0)
            : path(p), length(l), nodesize(n), score(s) {}
        std::vector<PathEdge*> path;
        int length{0};
        int nodesize {0};
        int score {0};
    };

    LinearPath FindBridgePath(PathEdge* start, int max_length, int max_nodesize);
    LinearPath FindBridgePath1(PathEdge* start, int max_length, int max_nodesize);
    bool HasBridgeJunction(const LinearPath& path, int max_depth);

    void RemoveChimeric();
    void IdentifyPaths(const std::string &method="no");
    std::list<PathEdge*> ExtendPath(PathEdge* e, std::unordered_set<PathEdge*> &visited, const std::string &method);
    
    template<typename TI, typename TO>
    std::list<PathEdge*> ExtendPathWithMethod(PathEdge* e, std::unordered_set<PathEdge*> &visited, TI get_in_edge, TO get_out_edge);

    std::vector<std::list<PathEdge*>>& GetPaths() { return paths_; }
 
    void SaveEdges(const std::string &fname);
    void SortPaths();
    void SavePaths(const std::string &fname);

    size_t PathLength(const std::list<PathEdge*> &path);

    const class AsmOptions& Options() const { return asmdata_.opts_; }

    void TestGraph();
    void InsertLoopNode(LoopNode* n);
    void InsertCrossNode(CrossNode* n);
protected:

    std::vector<std::list<PathEdge*>> paths_;
    std::unordered_map<std::string,  void (PathGraph::*)()> simplifications_;

public:
    std::unordered_map<SgNode::ID, PathNode*, SgNode::ID::Hash> &nodes_;
    std::unordered_map<PathEdge::ID, PathEdge*, PathEdge::ID::Hash> &edges_;

    StringGraph * string_graph_ { nullptr };
    

    struct Cluster {
        std::unordered_set<PathEdge*> edges;
        std::unordered_set<PathEdge*> nontrivial;;
        
        using LPath = std::pair<size_t, std::vector<PathEdge*>>;
        size_t length {0};
        size_t Size() const {
            size_t sz = 0;
            for (auto e : edges) {
                sz += e->Length();
            }
            return sz;
        }
        size_t Length() const { return length; }
        void FindLongest();
        void FindLongest0(PathEdge* e, std::unordered_set<PathEdge*> &visited, std::unordered_map<PathEdge*, LPath>& longests);
        void FindLongest1(PathEdge* e, std::unordered_set<PathEdge*> &visited, std::unordered_map<PathEdge*, LPath>& longests);
        std::vector<PathEdge*> GetStarts() const;
        std::vector<PathEdge*> GetEnds() const;
    };

    std::vector<Cluster> clusters_;
    std::unordered_map<PathEdge*, size_t> edge2cluster_;
};

} //namespace fsa {

