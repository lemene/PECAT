#include "minimizer_graph.hpp"
#include <inttypes.h>
#include <fstream>

#include "mkseq_store.hpp"

namespace fsa {

void MinimizerGraph::Build(const MkseqStore& mkseqs) {
    mkseqs_ = &mkseqs;
    for (size_t i = 0; i < mkseqs.Size(); ++i) {
        auto mks = mkseqs.Get(i);
        //assert(mks.Size() > 0);
        if (mks.Size() >= 2) {
            for (size_t j = 0; j + 1< mks.Size(); ++j) {
                redge_list_.push_back(new RawEdge(i, j, 0));
                redge_list_.push_back(new RawEdge(i, j, 1));
            }

        }
    }
    LOG(INFO)("RawEdgeSize=%zd", redge_list_.size());
    std::sort(redge_list_.begin(), redge_list_.end(), [this](const RawEdge* a, const RawEdge* b) {
        const auto a0 = a->InNode(*mkseqs_);
        const auto a1 = a->OutNode(*mkseqs_);
        const auto b0 = b->InNode(*mkseqs_);
        const auto b1 = b->OutNode(*mkseqs_); 

        return   a0 < b0 || (a0 == b0 && a1 < b1) || (a0 == b0 && a1 == b1 && a->seq_id_ < b->seq_id_);
    });

    size_t start = 0;
    while (start < redge_list_.size()) {
        size_t end = start + 1;
        while (end < redge_list_.size() && redge_list_[start]->Equal(*mkseqs_, *redge_list_[end])) {
            end ++;
        }
        edge_list_.push_back(new Edge(start, end));
        start = end;
    }
    LOG(INFO)("EdgeSize=%zd", edge_list_.size());

    // build index node_to_edge
    for (size_t i = 0; i < edge_list_.size(); ++i) {
        Node n = edge_list_[i]->InNode(*this);
        if (node_to_edge.find(n) == node_to_edge.end()) {
            node_to_edge[n] = i;
        }
    }
}

void MinimizerGraph::Save(const std::string& fname) const {
    std::ofstream ofile(fname);
    size_t k = mkseqs_->K();

    ofile << "Source,Target,Weight\n";

    for (auto &e : edge_list_) {
        const auto& n0 = redge_list_[e->start_]->GetInNode(*mkseqs_);
        const auto& n1 = redge_list_[e->start_]->GetOutNode(*mkseqs_);
        ofile << KmerId2String(n0.kmer, k) << "_" << (n0.dir^redge_list_[e->start_]->dir_) << "," << 
                 KmerId2String(n1.kmer, k) << "_" << (n1.dir^redge_list_[e->start_]->dir_) << "," << e->Count() << "\n";
    }
}

void MinimizerGraph::Neighbor(uint64_t h, uint8_t d, uint8_t deep, std::unordered_set<uint64_t>& hashs) {
    std::unordered_set<uint64_t> hash_values;
    //LOG(INFO)("Neighbor: %lld, %d, %d", h, d, deep);
    hashs.insert(h);
    if (deep > 0) {
        
        Node n(h, d, 0);
        auto es = OutEdges(n);
        LOG(INFO)("Neighbor(out_edge): %zd", es.size());
        for (size_t i = 0; i < es.size(); ++i) {
            if (es[i]->Count() > 3) {   // TODO PARAM
                auto on = es[i]->OutNode(*this);
                if (hashs.find(on.Hash()) == hashs.end() ) {
                    Neighbor(on.Hash(), on.Direction(), deep-1, hashs);
                    Neighbor(on.Hash(), !on.Direction(), deep-1, hashs);
                }
            }
        }
    } 

}

} // namespace fsa