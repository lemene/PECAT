#pragma once

#include <algorithm>
#include <set>
#include <vector>
#include <unordered_map>

#include "sequence.hpp"

namespace fsa {

template<typename E>
class TNode {
public:
    using ID = Seq::EndId;

public:
    TNode(ID id) : id_(id) {}
    ID Id() const { return id_; }
    size_t  InDegree() const { return  in_edges_.size(); }
    size_t OutDegree() const { return out_edges_.size(); }

    static ID ReverseId(ID id)  { return Seq::ReverseEndId(id); }
    void ReduceInEdge(const E *e) { ReduceEdge(in_edges_, reduced_in_edges_, e); }
    void ReduceOutEdge(const E *e) { ReduceEdge(out_edges_, reduced_out_edges_, e); }
    void ReduceEdge(std::vector<E*> &src, std::vector<E*> &dst, const E *e) {
        auto it = std::find(src.begin(), src.end(), e);    
        assert(it != src.end());
        dst.push_back(*it);
        src.erase(it);
    }

    ID id_;
    std::vector<E*> in_edges_;
    std::vector<E*> out_edges_;

    std::vector<E*> reduced_out_edges_;
    std::vector<E*> reduced_in_edges_;

    E* best_in_edge_{ nullptr };
    E* best_out_edge_{ nullptr };
    int mark_;

};

template <typename N>
class TEdge {

public:
    TEdge(N* in, N* out) : in_node_(in), out_node_(out) {}
    N* in_node_ { nullptr };
    N* out_node_ { nullptr };
};

template<typename N, typename E>
class TGraph {
public:
    std::unordered_map<int, N*> nodes_;
    std::unordered_map<typename E::ID, E*, typename E::Hash, typename E::Compare> edges_;
};


class MatrixGraph {
public:
    MatrixGraph(size_t sz) : size_(sz), edges_(sz*sz, 0) {}

    void AddEdge(int i, int j, int v) {
        edges_[i*size_+j] = v;
        edges_[j*size_+i]  = v;
    }
    size_t Size() const { return size_; }
    size_t Degree(int i) const;
    std::vector<int> Partition() const;
    std::vector<std::set<int>> Cluster() const;

    void Print() const;
protected:
    size_t size_;
    std::vector<int> edges_;
};

} // namespace fsa {
