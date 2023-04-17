#pragma once

#include <list>

#include "../sequence.hpp"
#include "utils/logger.hpp"

namespace fsa {

static const int ID_NODE_TYPE_STRING = 0;
static const int ID_NODE_TYPE_BASE = 1;
static const int ID_NODE_TYPE_MULT = 2;
static const int ID_NODE_TYPE_BUBBLE = 3;
static const int ID_NODE_TYPE_SEMI_BUBBLE = 4;
static const int ID_NODE_TYPE_LOOP = 5;
static const int ID_NODE_TYPE_CROSS = 6;


static const int ID_EDGE_TYPE_STRING = 0;
static const int ID_EDGE_TYPE_BASE = 1;
static const int ID_EDGE_TYPE_MULT = 2;
static const int ID_EDGE_TYPE_BUBBLE = 3;
static const int ID_EDGE_TYPE_SEMI_BUBBLE = 4;
static const int ID_EDGE_TYPE_LOOP = 5;
static const int ID_EDGE_TYPE_CROSS = 6;

static const int ID_EDGE_TYPE_ALT = 7;

struct SgNodeID {
    SgNodeID() { t = -1; }
    SgNodeID(int id, int type=0) {
        v.push_back(id);
        v.shrink_to_fit();
        t = type;
    }

    SgNodeID(int type, int v0, int v1) {
        t = type;
        v.push_back(v0);
        v.push_back(v1);
        v.shrink_to_fit();
    }

    SgNodeID(int type, const std::vector<int> &values) {
        t = type;
        v = values;
        v.shrink_to_fit();
    }
    
    SgNodeID(int type, const SgNodeID& id) {
        v = id.v;
        v.shrink_to_fit();
        t = type;
    }

    bool operator != (const SgNodeID& a) const {
        return ! operator==(a);
    }

    bool operator == (const SgNodeID& a) const {
        if (t != a.t) return false;
        if (v.size() != a.v.size()) return false;

        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i] != a.v[i]) return false;
        } 
        return true;
    }

    bool operator < (const SgNodeID &a) const {
        if (t < a.t) return true;
        for (size_t i = 0; i < std::min(v.size(), a.v.size()); ++i) {
            if (v[i] < a.v[i]) return false;
        } 
        return v.size() < a.v.size();
    }

    struct Hash {
        size_t operator() (const SgNodeID &a) const {
            return a.hash();
        }
    };

    size_t hash() const {
        size_t h = 17;
        h = h * 31 + std::hash<int>()(t);
        
        for (size_t i = 0; i < v.size(); ++i) {
            h = h * 31 + std::hash<int>()(v[i]);
        }
        return h;
    }

    int End() const { return v[0] > 0 ? 0 : 1; }
    int MainNode() const { return v[0]; }
    int Type() const { return t; }
    const std::vector<int>& Values() const { return v; }
    int Value(size_t i) const { return v[i]; }
    size_t ValueSize() const { return v.size(); }
    static SgNodeID Reverse(SgNodeID id) {
        SgNodeID r = id;
        r.Reverse();
        return r;
    }

    void Reverse() {
        std::reverse(v.begin(), v.end());
        std::for_each(v.begin(), v.end(), [](int& c){ c = Seq::ReverseEndId(c); });
    }

    std::string ToString(const class StringPool &sp) const;
    std::string ToString() const;

    //int v;
protected:
    int t;
    std::vector<int> v;

};

class SgEdge;

class SgNode {
public:
    virtual ~SgNode() {}
    using ID = SgNodeID;
    
    const ID& Id() const { return id_; }
    ID RId() const { return ID::Reverse(id_); }
    size_t  InDegree() const { return  org_in_edges_.size(); }
    size_t OutDegree() const { return org_out_edges_.size(); }

    void AddInEdge(SgEdge* e) { org_in_edges_.push_back(e); }
    void AddOutEdge(SgEdge* e) { org_out_edges_.push_back(e); }

    SgEdge* InEdge(size_t i) { return org_in_edges_[i]; }
    const SgEdge* InEdge(size_t i) const { return org_in_edges_[i]; }
    SgEdge* OutEdge(size_t i) { return org_out_edges_[i]; }
    const SgEdge* OutEdge(size_t i) const { return org_out_edges_[i]; }
    template<typename T>
    T* InEdge(size_t i) { return static_cast<T*>(org_in_edges_[i]);}
    template<typename T>
    T* OutEdge(size_t i) { return static_cast<T*>(org_out_edges_[i]);}

    void ReduceInEdge(const SgEdge *e) { MoveEdge(org_in_edges_, org_reduced_in_edges_, e); }
    void ReduceOutEdge(const SgEdge *e) { MoveEdge(org_out_edges_, org_reduced_out_edges_, e); }
    void ReactivateInEdge(const SgEdge *e) { MoveEdge(org_reduced_in_edges_, org_in_edges_, e); }
    void ReactivateOutEdge(const SgEdge *e) { MoveEdge(org_reduced_out_edges_, org_out_edges_, e); }

    SgNode* InNode(size_t i);
    const SgNode* InNode(size_t i) const;
    SgNode* OutNode(size_t i);
    const SgNode* OutNode(size_t i) const ;
    template<typename T>
    T* InNode(size_t i) { return static_cast<T*>(InNode(i));}
    template<typename T>
    T* OutNode(size_t i) { return static_cast<T*>(OutNode(i));}

    virtual bool IsConsistOf(const SgNode* n) const { return false;}


    void AttachType(const std::string& t) { type_ += ":"; type_ += t; }
    const std::string& Type() const { return  type_; }
    bool IsType(const std::string &t) const { return type_.find(t) != std::string::npos; }

    bool operator < (const SgNode& n) const { return id_ < n.id_; }
protected: 
    void MoveEdge(std::vector<SgEdge*> &src, std::vector<SgEdge*> &dst, const SgEdge *e) {
        auto it = std::find(src.begin(), src.end(), e);    
        assert(it != src.end());
        dst.push_back(*it);
        src.erase(it);
    }
protected:
    ID id_;
    std::string type_ { "sg" };

    std::vector<SgEdge*> org_in_edges_;
    std::vector<SgEdge*> org_out_edges_;

    std::vector<SgEdge*> org_reduced_in_edges_;
    std::vector<SgEdge*> org_reduced_out_edges_;
};

class BaseEdge;


class BaseNode : public SgNode {
public:
    BaseNode(ID id) 
        : in_edges_(reinterpret_cast<std::vector<BaseEdge*>&>(org_in_edges_))
        ,  out_edges_(reinterpret_cast<std::vector<BaseEdge*>&>(org_out_edges_))
        ,  reduced_in_edges_(reinterpret_cast<std::vector<BaseEdge*>&>(org_reduced_in_edges_))
        ,  reduced_out_edges_(reinterpret_cast<std::vector<BaseEdge*>&>(org_reduced_out_edges_))
    { id_ = id; AttachType("base");}

    std::vector<BaseNode*> GetAllOutNodes() const; 
    std::vector<BaseNode*> GetAllInNodes() const;
    const std::vector<BaseEdge*>& GetInEdges() const { return in_edges_; }
    std::vector<BaseEdge*>& GetInEdges() { return in_edges_; }
    const std::vector<BaseEdge*>& GetOutEdges() const { return out_edges_; }
    std::vector<BaseEdge*>& GetOutEdges() { return out_edges_; }

    const BaseEdge* GetInEdge(size_t i) const { return in_edges_[i]; }
    BaseEdge* GetInEdge(size_t i) { return in_edges_[i]; }
    const BaseEdge* GetOutEdge(size_t i) const { return out_edges_[i]; }
    BaseEdge* GetOutEdge(size_t i) { return out_edges_[i]; }

    static ID ReverseId(ID id)  { return SgNodeID::Reverse(id); }

    int ReadId() const { return Seq::EndIdToId(id_.MainNode());  }
    BaseEdge* GetBestInEdge() const ;
    BaseEdge* GetBestOutEdge() const ;
    const std::vector<BaseEdge*> & GetReducedOutEdge() const {
        return reduced_out_edges_;
    }

    static const std::string TypeName() { return "base"; };
    int mark_;
protected:

protected:

    std::vector<BaseEdge*>& in_edges_;
    std::vector<BaseEdge*>& out_edges_;

    std::vector<BaseEdge*>& reduced_in_edges_;
    std::vector<BaseEdge*>& reduced_out_edges_;
};


class PathEdge;
class PathGraph;
class PathNode : public SgNode {
    friend class PathEdge;
    friend class PathGraph;
public:
    PathNode(BaseNode *n, int type=0) : SgNode(), string_node_(n) 
        ,  in_edges_(reinterpret_cast<std::vector<PathEdge*>&>(org_in_edges_))
        ,  out_edges_(reinterpret_cast<std::vector<PathEdge*>&>(org_out_edges_))
        ,  reduced_in_edges_(reinterpret_cast<std::vector<PathEdge*>&>(org_reduced_in_edges_))
        ,  reduced_out_edges_(reinterpret_cast<std::vector<PathEdge*>&>(org_reduced_out_edges_))
    { 
        if (n != nullptr) id_ = CreateId(n->Id(), type); 
        AttachType(TypeName()); 
    }
    virtual ~PathNode() {}

    static std::string TypeName() { return "path"; }

    static ID CreateId(BaseNode::ID id, int type = 0) {
        return BaseNode::ID(type, id);
    }

    static int GetIdValue(ID id) {
        return id.MainNode();
    }

    static ID  ReverseId(ID id) {
        return ID::Reverse(id);
    }

    const std::vector<PathEdge*>& GetInEdges() const { return in_edges_; }
    std::vector<PathEdge*>& GetInEdges() { return in_edges_; }
    const std::vector<PathEdge*>& GetOutEdges() const { return out_edges_; }
    std::vector<PathEdge*>& GetOutEdges() { return out_edges_; }

    const PathEdge* GetInEdge(size_t i) const { return in_edges_[i]; }
    PathEdge* GetInEdge(size_t i) { return in_edges_[i]; }
    const PathEdge* GetOutEdge(size_t i) const { return out_edges_[i]; }
    PathEdge* GetOutEdge(size_t i) { return out_edges_[i]; }

    PathEdge* GetBestInEdge() const;
    PathEdge* GetBestOutEdge() const;

    const BaseNode* GetBaseNode() const { return string_node_; }

protected:
    BaseNode* string_node_ { nullptr };
    
    std::vector<PathEdge*>& in_edges_;
    std::vector<PathEdge*>& out_edges_;

    std::vector<PathEdge*>& reduced_in_edges_;
    std::vector<PathEdge*>& reduced_out_edges_;

    PathEdge* best_in_edge_{ nullptr };
    PathEdge* best_out_edge_{ nullptr };

};

class LoopNode : public PathNode {
public:
    LoopNode(const std::vector<SgEdge*> &backward);
    LoopNode(const std::vector<SgEdge*> &backward, const std::vector<SgEdge*> &forward);

    LoopNode* Reverse(PathGraph& sg) const;

    size_t OriginOutDegree() const;
    size_t OriginInDegree() const;
    const SgEdge* OriginOutEdge(size_t i) const;
    const SgEdge* OriginInEdge(size_t i) const;
    virtual bool IsConsistOf(const SgNode* n) const ;
    void ReduceSub();
    std::vector<SgEdge*>& GetBackward() { return backward_; }
protected:
    SgNode::ID CreateID() const;
    std::vector<SgEdge*> backward_;
    std::vector<SgEdge*> forward_;
};

}

