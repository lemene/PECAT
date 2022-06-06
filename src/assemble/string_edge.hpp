#pragma once

#include <array>
#include <list>

#include "../overlap.hpp"
#include "string_node.hpp"
#include "../utility.hpp"
#include "../utils/logger.hpp"

namespace fsa {

class ReadStore;
class PhaseInfoFile;
class StringGraph;
class StringPool;

struct SgEdgeState {

};

struct SgEdgeID {
    SgEdgeID() {}
    SgEdgeID(int t, const SgNode::ID &a, const SgNode::ID& b) : type_(t), values_({a, b}){
        values_.shrink_to_fit();
    }
    
    SgEdgeID(int t, const SgNode::ID &a0, const SgNode::ID &a1, const SgNode::ID &b1, const SgNode::ID& b0) 
     : type_(t), values_({a0, a1, b1, b0}) {

        values_.shrink_to_fit();
    }

    SgEdgeID(int t, const std::vector<BaseNode::ID> &vs) : type_(t), values_(vs) {
        assert(vs.size() == 2 || vs.size() == 4);
        values_.shrink_to_fit();
    }

    bool IsMultiple() const { assert(values_.size() == 2 || values_.size() == 4); return values_.size() == 4;}
    int Type() const { return type_; }
    const SgNodeID& Value(size_t i) const { return values_[i]; }
    size_t ValueSize() const { return values_.size(); }
    struct Hash {
        size_t operator()(const SgEdgeID& r) const { 
            return r.hash();
        }
    };

    static SgEdgeID Reverse(const SgEdgeID& a) {
        SgEdgeID b = a;
        b.Reverse();
        return b;
    }

    void Reverse() {
        std::for_each(values_.begin(), values_.end(), [](BaseNode::ID& i) {
            i.Reverse();
        });
        std::reverse(values_.begin(), values_.end());
    }

    bool operator == (const SgEdgeID &a) const {
        if (Type() != a.Type()) return false;
        if (ValueSize() != a.ValueSize()) return false;
        for (size_t i = 0; i < ValueSize(); ++i) {
            if (Value(i) != a.Value(i)) return false;
        }

        return true;
    }
    size_t hash() const {
        size_t h = 17;
        h = h * 31 + std::hash<int>()(type_);
        
        for (size_t i = 0; i < values_.size(); ++i) {
            h = h * 31 + values_[i].hash();
        }
        return h;
    }

    std::string ToString(const class StringPool &sp) const;
    std::string ToString() const;

protected:
    int type_ { 0 };
    std::vector<SgNodeID> values_; // 
};


class SgEdge {
public:
    using ID = SgEdgeID;
public:
    SgEdge()  { }
    virtual ~SgEdge() {}
    SgEdge(SgNode* in_node, SgNode *out_node) : in_node_(in_node), out_node_(out_node) { }

    const ID& Id() const { return id_; }
    void SetInNode(SgNode* node) { in_node_ = node;}
    void SetOutNode(SgNode* node) { out_node_ = node;}

    SgNode* InNode() { return in_node_; }
    const SgNode* InNode() const { return in_node_; }
    SgNode* OutNode() { return out_node_; }
    const SgNode* OutNode() const { return out_node_; }

    template<typename T>
    T* OutNode() { return static_cast<T*>(out_node_); }
    template<typename T>
    const T* OutNode() const { return static_cast<T*>(out_node_); }
    template<typename T>
    T* InNode() { return static_cast<T*>(in_node_); }
    template<typename T>
    const T* InNode() const { return static_cast<T*>(in_node_); }

    virtual size_t Score() const { return 0; }
    virtual size_t Length() const { return 0; }
    virtual size_t NodeSize() const { return 0; } 

    virtual void ReduceSubEdges() {}

    virtual std::string ToString(const StringPool &sp) const { return ""; }

    void AttachType(const std::string& t) { type += ":"; type += t; }
    const std::string& Type() const { return  type; }
    bool IsType(const std::string &t) { return type.find(t) != std::string::npos; }

protected:
    std::string type {"sg"};
    ID id_;
    SgNode* in_node_ { nullptr };
    SgNode* out_node_ { nullptr };
};

class BaseEdge : public SgEdge {
public:
    using ID = SgEdgeID;
public:
    enum  ReduceType {
        RT_ACTIVE,
        RT_SPUR,
        RT_CHIMER,
        RT_REMOVED,
        RT_REPEAT,
        RT_TRANSITIVE,
        RT_NO_BEST,
        RT_UNRELIABLE,
        RT_BRIDGED,
        RT_TYPE1,
        RT_PHASED,
        RT_LOW_QUALITY,
        RT_UNKNOWN,
    };

    static const char* ReduceTypeName(ReduceType t) {
        switch (t) {
        case RT_ACTIVE:     return "active";
        case RT_SPUR:       return "spur";
        case RT_CHIMER:     return "chimer";
        case RT_REMOVED:    return "removed";
        case RT_TRANSITIVE: return "transitive";
        case RT_NO_BEST:    return "no_best";
        case RT_UNRELIABLE: return "unreliable";
        case RT_BRIDGED:    return "bridged";
        case RT_TYPE1:      return "type1";
        case RT_PHASED:      return "phased";
        case RT_LOW_QUALITY:      return "low_quality";
        case RT_REPEAT:      return "repeat";
        case RT_UNKNOWN: default:   return "unknown";
        }
    }

public:
    BaseEdge(BaseNode* in_node, BaseNode *out_node) : SgEdge(in_node, out_node) { 
        id_ = ID(0, InNode()->Id(), OutNode()->Id());
        AttachType("base"); 
    }

    BaseNode* InNode() { return (BaseNode*)in_node_; }
    const BaseNode* InNode() const { return (BaseNode*)in_node_; }
    BaseNode* OutNode() { return (BaseNode*)out_node_; }
    const BaseNode* OutNode() const { return (BaseNode*)out_node_; }

    static bool Less(const BaseEdge* a, const BaseEdge *b) { return a->Score() < b->Score(); }
    
    static ID ReverseID(ID id) { return ID::Reverse(id); }

    static ID CreateID(BaseNode::ID id0, BaseNode::ID id1) {  return ID(0, id0, id1); }

    ReduceType GetReduceType() const { return reduce_type_; }
    const char* GetReduceTypeName() const { return ReduceTypeName(reduce_type_); }
    

    Seq::Tile GetTile() const;
    bool IsReduce() const { return reduce_type_ != RT_ACTIVE; }
    void Reduce(ReduceType t) {
        OutNode()->ReduceInEdge(this);
        InNode()->ReduceOutEdge(this);
        reduce_type_ = t;
    }

    void Reactivate();

    double Identity() const;
    // edge score, the larger the score, the better the link 
    virtual size_t Score() const { return ol_->AlignedLength(); }
    virtual size_t  Length() const;
    virtual size_t NodeSize() const { return 1; }

    bool TestOutExtend(int minlen, int minnode, bool both=true);
    bool TestInExtend(int minlen, int minnode, bool both=true);

    const Overlap::Read* InRead() const { return read_ == ol_->a_.id ? &(ol_->b_) : &(ol_->a_); }
    const Overlap::Read* OutRead() const { return read_ == ol_->a_.id ? &(ol_->a_) : &(ol_->b_); }

    int read_ { 0 };
    ReduceType reduce_type_{ RT_ACTIVE };
    const Overlap* ol_ { nullptr } ;
    std::array<int, 2> score_;
    bool subject_ { true };

protected:
};


class PathEdge : public SgEdge {
public:
    using ID = SgEdgeID;
    enum ReduceType {
        RT_OK,
        RT_CONTAINED,
        RT_UNKNOWN,
        REPEAT_BRIDGE,
        SPUR,
        SIMPLE_DUP,
        CONTAINED,
        CROSS    
    };
public:
    PathEdge(SgNode* in, SgNode* out) : SgEdge(in, out) {
        AttachType("path"); 
    }
    virtual ~PathEdge() {}

    PathNode* InNode() { return (PathNode*)in_node_; }
    const PathNode* InNode() const { return (PathNode*)in_node_; }
    PathNode* OutNode() { return (PathNode*)out_node_; }
    const PathNode* OutNode() const { return (PathNode*)out_node_; }

    static ID ReverseId(ID id) {
        return ID::Reverse(id);
    }

    ID ReverseId() const { return ReverseId(Id()); }
    virtual void ReduceStringEdge() = 0;

    virtual size_t SimplePathSize() const = 0;
    virtual const std::list<BaseEdge*>& GetSimplePath(size_t i) const = 0;
    virtual size_t NodeSize() const  = 0;
    virtual size_t MinPathSize() const = 0;
    virtual bool IsPathValidated() const { return validated_; }

    virtual size_t Score() const { assert(score_ >= 0); return score_; }
    virtual size_t Length() const { return length_; }
    void Reduce(const char* reduce_type, bool remove_string_edge=false) {
        if (!IsReduced()) { 
            reduce_type_ = RT_UNKNOWN;
            type_ = reduce_type;

            if (remove_string_edge) ReduceStringEdge();
            OutNode()->ReduceInEdge(this);
            InNode()->ReduceOutEdge(this);
        }
    }

    bool IsReduced() const { return reduce_type_ != RT_OK ;}

    // identify paths which can constract contigs
    virtual void IdentifySimplePaths(StringGraph& string_graph) {};
    virtual bool Contain(const BaseEdge* e) const = 0;
    virtual std::unordered_set<Seq::Id> GetReads() = 0;


    int length_{ 0 };
    double width_{ 0.0 };
    int score_{ 0 };
    std::string type_{ "" };
    ReduceType reduce_type_ { RT_OK };
    bool validated_ { false };

protected:
};

class SimplePathEdge : public PathEdge {
public:
    SimplePathEdge(SgNode* in, SgNode* out, const std::list<BaseEdge*>& path)
        : PathEdge(in, out), path_(path) {
        score_ = std::accumulate(path.begin(), path.end(), 0, [](int a, BaseEdge* b) { return a + b->Score(); });
        length_ = std::accumulate(path.begin(), path.end(), 0, [](int a, BaseEdge* b) { return a + b->Length(); });
        width_ = 1.0;
        AttachType("simple");
        id_ = ID(ID_EDGE_TYPE_MULT, { InNode()->Id(), PathNode::CreateId(path_.front()->OutNode()->Id()), 
            PathNode::CreateId(path_.back()->InNode()->Id()), OutNode()->Id()});
    }
    virtual ~SimplePathEdge()  {}

    virtual bool Contain(const BaseEdge* e) const  {
        return std::find(path_.begin(), path_.end(), e) != path_.end();
    }

    virtual size_t SimplePathSize() const  { return 1; }
    virtual const std::list<BaseEdge*>& GetSimplePath(size_t i) const { return path_; }
    virtual size_t NodeSize() const { return path_.size(); }
    virtual size_t MinPathSize() const { return path_.size(); }
    virtual std::unordered_set<Seq::Id> GetReads() {
        std::unordered_set<Seq::Id> reads;
        for (const auto p : path_) {
            reads.insert(p->InNode()->ReadId());
            reads.insert(p->OutNode()->ReadId());
        }
        return reads;
    }

    virtual void ReduceStringEdge() {
        for (auto e : path_) {
            e->Reduce(BaseEdge::RT_TYPE1);
        }
    }
    virtual std::string ToString(const StringPool &sp) const;
    std::list<BaseEdge*> path_;
};

class BubbleEdge : public PathEdge {
public:
    BubbleEdge(SgNode *in, SgNode *out, std::list<PathEdge*> &path, int length, double width, int score)
        : PathEdge(in, out), simple_paths_(path)
    {
        AttachType("bubble");
        length_ = length;
        score_ = score;
        width_ = width;
        id_ = ID(ID_EDGE_TYPE_BUBBLE, InNode()->Id(), OutNode()->Id());
    }
    virtual ~BubbleEdge()  {}
    virtual std::string ToString(const StringPool &sp) const;

    virtual bool Contain(const BaseEdge* e) const {
        for (auto p : simple_paths_) {
            if (p->Contain(e)) return true;
        }
        return false;
    }

    virtual size_t SimplePathSize() const  { return string_edges_.size(); }
    virtual const std::list<BaseEdge*>& GetSimplePath(size_t i) const { return string_edges_[i]; }
    virtual void IdentifySimplePaths(StringGraph& string_graph); 
    virtual size_t NodeSize() const { 
        return std::accumulate(string_edges_.begin(), string_edges_.end(), 0, [](size_t a, const std::list<BaseEdge*>& b) {
            return a + b.size();
        }); 
    }
    virtual size_t MinPathSize() const {
        std::vector<size_t> sz(string_edges_.size());
        std::transform(string_edges_.begin(), string_edges_.end(), sz.begin(), [](const std::list<BaseEdge*>& path) {
            return path.size();
        });
        auto m = std::min_element(sz.begin(), sz.end());
        return m == sz.end() ? 0 : *m;
    }

    bool Validate(PhaseInfoFile *ignored, const ReadStore &rs, StringGraph *sg, int max_length);


    virtual std::unordered_set<Seq::Id> GetReads() {
        std::unordered_set<Seq::Id> reads;
        for (const auto p : simple_paths_) {
            auto r = p->GetReads();
            reads.insert(r.begin(), r.end());
        }
        return reads;
    }

    virtual void ReduceStringEdge() {
        for (auto p : simple_paths_) {
            p->ReduceStringEdge();
        }
    }

    virtual void ReduceSubEdges() {
        for (auto p : simple_paths_) {
            p->Reduce("Contained");
        }
    }

    std::list<PathEdge*> simple_paths_;
    std::vector<std::list<BaseEdge*>> string_edges_;

};

class SemiBubbleEdge : public PathEdge {
public:
    SemiBubbleEdge(const std::vector<std::vector<PathEdge*>> &paths)
     : PathEdge(paths[0].front()->InNode(), paths[0].back()->OutNode()), paths_(paths) {
        AttachType("semi");
        id_ = ID(ID_EDGE_TYPE_SEMI_BUBBLE, InNode()->Id(), OutNode()->Id());

        for (auto &p : paths) {
            auto len = std::accumulate(p.begin(), p.end(), 0, [](int a, const PathEdge* b) { return a + b->Length(); });
            length_ = std::max(len, length_);
        }
    }

    virtual ~SemiBubbleEdge()  {}    
    virtual std::string ToString(const StringPool &sp) const;

    virtual void ReduceStringEdge();
    virtual size_t SimplePathSize() const { return string_edges_.size(); }
    virtual const std::list<BaseEdge*>& GetSimplePath(size_t i) const { return string_edges_[i]; };
    virtual size_t NodeSize() const { 
        return std::accumulate(string_edges_.begin(), string_edges_.end(), 0, [](size_t a, const std::list<BaseEdge*>& b) {
            return a + b.size();
        }); 
    }

    virtual size_t MinPathSize() const {
        std::vector<size_t> sz(string_edges_.size());
        std::transform(string_edges_.begin(), string_edges_.end(), sz.begin(), [](const std::list<BaseEdge*>& path) {
            return path.size();
        });
        auto m = std::min_element(sz.begin(), sz.end());
        return m == sz.end() ? 0 : *m;
    }

    virtual bool Contain(const BaseEdge* e) const;
    virtual std::unordered_set<Seq::Id> GetReads();

    bool Validate(PhaseInfoFile *ignored, const ReadStore &rs, StringGraph *sg, int max_length);

    virtual void ReduceSubEdges() {
        for (const auto &path : paths_) {
            for (const auto &p : path) {
                p->Reduce("Contained");
            }
        }
    }    
    virtual void IdentifySimplePaths(StringGraph& string_graph); 

    uint8_t Branchs() const { // 1 head 2 tail 3 head+tail
        assert(paths_.size() >= 2);
        uint8_t r = 0;
        for (size_t i = 1; i < paths_.size(); ++i) {
            if (paths_[i].front()->InNode() == InNode()) {
                r |= 1;
            }

            if (paths_[i].back()->OutNode() == OutNode()) {
                r |= 2;
            }
        }
        return r;
    }

    // create new edges which is reversed
    SemiBubbleEdge* Reverse(PathGraph &graph);

    void Merge(SemiBubbleEdge& e) {
        auto b = e.Branchs();
        assert(b == 1 || b == 2);
        if ((b & Branchs()) == 0)  {
            paths_.insert(paths_.end(), e.paths_.begin()+1, e.paths_.end());
            for (const auto &path : paths_) {
                for (const auto &p : path) {
                    p->Reduce("Contained");
                }
            }
        } else {
            // Merged do nothing
        }

        length_ = std::max(length_, e.length_);

    }


    std::vector<std::vector<PathEdge*>> paths_;
    
    std::vector<std::list<BaseEdge*>> string_edges_;
};

class LoopEdge : public PathEdge {

public:
    LoopEdge(SgNode* in_node, SgNode* out_node, const std::vector<SgEdge*> &forward, const std::vector<SgEdge*>& backward)
     : PathEdge(in_node, out_node), forward_(forward), backward_(backward) {
         // TODO remove duplicate
        AttachType("loop");

        id_ = ID(ID_EDGE_TYPE_LOOP, InNode()->Id(),  OutNode()->Id());
        length_ += std::accumulate(forward.begin(), forward.end(), 0, [](int a, const SgEdge* b) { return a + b->Length(); });
        length_ += std::accumulate(backward.begin(), backward.end(), 0, [](int a, const SgEdge* b) { return a + b->Length(); });
    }

    virtual ~LoopEdge()  {}
    
    virtual std::string ToString(const StringPool &sp) const;

    virtual void ReduceStringEdge();
    virtual size_t SimplePathSize() const { return string_edges_.size(); }
    virtual const std::list<BaseEdge*>& GetSimplePath(size_t i) const { return string_edges_[i]; };
    virtual size_t NodeSize() const { 
        return std::accumulate(string_edges_.begin(), string_edges_.end(), 0, [](size_t a, const std::list<BaseEdge*>& b) {
            return a + b.size();
        }); 
    }
    
    virtual size_t MinPathSize() const {
        std::vector<size_t> sz(string_edges_.size());
        std::transform(string_edges_.begin(), string_edges_.end(), sz.begin(), [](const std::list<BaseEdge*>& path) {
            return path.size();
        });
        auto m = std::min_element(sz.begin(), sz.end());
        return m == sz.end() ? 0 : *m;
    }

    virtual bool Contain(const BaseEdge* e) const;
    virtual std::unordered_set<Seq::Id> GetReads();

    virtual void ReduceSubEdges() {
        for (const auto &p : forward_) {
            static_cast<PathEdge*>(p)->Reduce("Contained");
        }
        for (const auto &p : backward_) {
            static_cast<PathEdge*>(p)->Reduce("Contained");
        }
    }    
     
    virtual void IdentifySimplePaths(StringGraph& string_graph); 

    // create new edges which is reversed
    LoopEdge* Reverse(PathGraph &graph);
    void AddToGraph(PathGraph &path);
    
    std::vector<SgEdge*> forward_;
    std::vector<SgEdge*> backward_;
    
    std::vector<std::list<BaseEdge*>> string_edges_;
};

class AAAEdge : public PathEdge {
public:
    AAAEdge(PathEdge* e, SgNode* innode, SgNode* outnode)
        : PathEdge(innode, outnode), origin(e) {
        
        AttachType("alt");
        assert(origin != nullptr);

        if (out_node_ == nullptr) {
            out_node_ = origin->OutNode();
        }
        if (in_node_ == nullptr) {
            in_node_ = origin->InNode();
        }
        if (e->Id().IsMultiple()) {
            id_ = ID(ID_EDGE_TYPE_ALT, in_node_->Id(), e->Id().Value(1), e->Id().Value(2), out_node_->Id());
        } else {
            id_ = ID(ID_EDGE_TYPE_ALT, in_node_->Id(), out_node_->Id());

        }
    }
    virtual ~AAAEdge() {}

    virtual void ReduceSubEdges() {
        origin->Reduce("Contained");
    }   
    virtual void ReduceStringEdge() { origin->ReduceStringEdge(); };    
    virtual std::string ToString(const StringPool &sp) const { return "alt"; }

    virtual void IdentifySimplePaths(StringGraph& sg);


    virtual size_t SimplePathSize() const { return string_edges_.size(); };
    virtual const std::list<BaseEdge*>& GetSimplePath(size_t i) const { return string_edges_[0]; };
    virtual size_t NodeSize() const  { return origin->NodeSize(); };
    virtual size_t Length() const { return origin->Length(); }
    virtual size_t MinPathSize() const { return origin->MinPathSize(); }
    virtual bool Contain(const BaseEdge* e) const  { return origin->Contain(e); };
    virtual std::unordered_set<Seq::Id> GetReads() { return origin->GetReads(); }

    PathEdge* origin;
    std::vector<std::list<BaseEdge*>> string_edges_;
};
}