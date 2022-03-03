#ifndef FSA_CONTIG_GRAPH_HPP
#define FSA_CONTIG_GRAPH_HPP

#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <deque>

#include "../graph.hpp"
#include "utility.hpp"

#include "contig_link_store.hpp"
#include "sequence.hpp"

namespace fsa {

class ContigEdge;
class ContigNode : public TNode<ContigEdge>{
public:
    ContigNode(ID id) : TNode(id) {}
};


class ContigEdge : public TEdge<ContigNode> {
public:
    using Hash = ArrayHash<int,2>;
    using Compare = ArrayEqual<int,2>;
    using ID = std::array<ContigNode::ID, 2>;   

public:
    ContigEdge(ContigNode* in, ContigNode* out) : TEdge(in, out) {}

    ID Id() const { return std::array<int, 2>{in_node_->id_, out_node_->id_}; }
    static ID ReverseId(ID id)  {
        return std::array<int, 2>{ContigNode::ReverseId(id[1]), ContigNode::ReverseId(id[0]) };
    }

    std::vector<Seq::Tile> GetSeqTiles() const;
    size_t LinkLength() const { return link_->TotalLength(); }

    size_t Length() const { 
        //printf("length %d %d\n",link_->Target() , Seq::EndIdToId(in_node_->id_));
        if (link_->Target() == Seq::EndIdToId(in_node_->id_))  {
            return link_->TotalLength() - (int)link_->TargetLength();
        } else {
            return link_->TotalLength() - (int)link_->SourceLength();
        }
    }
    int GapLength() const {
        return (int)link_->TotalLength() - (int)link_->TargetLength() - (int)link_->SourceLength();
    }

    double Support() const { return link_->Score(); }

    void SetCovered(const std::array<const ContigEdge*,2> &c) {covered_ = c;}
    void Remove() {
        removed  = true; 
        auto it = std::find(in_node_->out_edges_.begin(), in_node_->out_edges_.end(), this);
    
        assert(it != in_node_->out_edges_.end());
        //dst.push_back(*it);
        in_node_->out_edges_.erase(it);

        
        it = std::find(out_node_->in_edges_.begin(), out_node_->in_edges_.end(), this);
    
        assert(it != out_node_->in_edges_.end());
        //dst.push_back(*it);
        out_node_->in_edges_.erase(it);
    }

    void Reduce() {
        out_node_->ReduceInEdge(this);
        in_node_->ReduceOutEdge(this);
        //reduce_type_ = t;
        removed = true;
    }

    ContigLink* link_;
    std::array<const ContigEdge*,2> covered_ { {nullptr, nullptr} };
    bool removed = {false};
};

class ContigGraph : public TGraph<ContigNode, ContigEdge>  {
public:
    ContigGraph(ContigLinkStore &links);
    ~ContigGraph();

    void Create();

    void AddLink(ContigLink& link);
    void AddEdge(Seq::EndId s, Seq::EndId t, ContigLink& link);

    void Transitive();
    void CheckRepeat();
    void RemoveCoveredEdges();
    void SolveRepeat();

    void IdentifyPaths(const std::string &method);
    std::unordered_set<Seq::Id> CollectSeqIdsInPaths() const;
    const std::unordered_set<Seq::Id>& GetContained() const { return contained_; }

    std::deque<ContigNode*> ExtendPath(ContigNode* e, std::unordered_set<ContigNode*> &visitied, const std::string &method);

    template<typename I, typename O>
    std::deque<ContigNode*> ExtendPathWithMethod(ContigNode* e, std::unordered_set<ContigNode*> &visited, I get_in_node, O get_out_node);

    ContigEdge* GetEdge(ContigNode* in, ContigNode* out) const {
        for (ContigEdge* e : in->out_edges_) {
            if (e->out_node_ == out) return e;
        }
        return nullptr;
    }

    ContigEdge* ReverseEdge(ContigEdge* e) { return edges_[ContigEdge::ReverseId(e->Id())]; }
    ContigNode* ReverseNode(ContigNode* n) { return nodes_[ContigNode::ReverseId(n->Id())]; }

    void CalucateBest(const std::string &method);
    void Dump(const std::string &fname, const ReadStore& rd_store);

    const std::list<std::deque<ContigNode*>>& GetPaths() { return paths_; }
    std::string ConstructContig(const std::list<ContigEdge*> &path);

protected:
    ContigLinkStore &links_;
  
    std::list<std::deque<ContigNode*>> paths_;
    std::unordered_set<Seq::Id> contained_;
};

#endif // FSA_CONTIG_GRAPH_HPP  


} // namespace fsa {