#include "contig_graph.hpp"

#include <unordered_set>

#include "../align/match_info.hpp"
#include "contig_analyzer.hpp"

namespace fsa {

ContigGraph::~ContigGraph() {
    for (auto it : nodes_) {
        delete it.second;
    }
    nodes_.clear();

    for (auto it : edges_) {
        delete it.second;
    }
    edges_.clear();
}

void ContigGraph::Build() {
    // 创建图的节点
    for (size_t i = 0; i < fragments_.size(); ++i) {
        auto n0_id = MakeNodeId(i, 0);
        auto n1_id = MakeNodeId(i, 1);
        nodes_[n0_id] = new Node(n0_id, &fragments_[i], 0);
        nodes_[n1_id] = new Node(n1_id, &fragments_[i], 1);
        LOG(INFO)("Fragments %zd %zd-%zd", i, fragments_[i].Start(), fragments_[i].End());
    }

    // 创建默认边
    for (size_t i = 1; i < fragments_.size(); ++i) {
        auto& frg0 = fragments_[i-1];
        auto& frg1 = fragments_[i];

        if (frg0.ContigId() == frg1.ContigId()) {
            auto n0_id = MakeNodeId(i-1, 1);
            auto n1_id = MakeNodeId(i, 0);
            auto n0 = nodes_[n0_id];
            auto n1 = nodes_[n1_id];

            Edge* e = new Edge(n0, n1);
            edges_[MakeEdgeId(n0_id, n1_id)] = e;

            n0->AddEdge(e);
            n1->AddEdge(e);
        }
    }

    // 检查节点是否通过读数直接连接
    std::vector<std::vector<const MatchInfo*>> end_covs(fragments_.size()*2);
    for (size_t i = 0; i < fragments_.size()*2; ++i) {
        auto n_id = MakeNodeId(i/2, i%2);
        assert(nodes_.find(n_id) != nodes_.end());

        const int FLAKING = 3000;
        auto& n = nodes_[n_id];
        if (i % 2 == 0) {
            //end_covs[i] = n->Fragment()->Analyzer()->GetCoverage(n->Fragment()->Start(), FLAKING);
        } else {
            //end_covs[i] = n->Fragment()->Analyzer()->GetCoverage(n->Fragment()->End(), -FLAKING);
        }
        LOG(INFO)("End-Cov %zd %08X %zd-%zd %zd", i, n_id, n->Fragment()->Start(), n->Fragment()->End(), end_covs[i].size());
    }

    for (size_t i0 = 0; i0 < fragments_.size()*2; ++i0) {
        auto n0_id = MakeNodeId(i0/2, i0%2);
        auto& n0 = nodes_[n0_id];
        for (size_t i1 = i0+1; i1 < fragments_.size()*2; ++i1) {
            if (i0 / 2 == i1 / 2) continue; 
            auto n1_id = MakeNodeId(i1/2, i1%2);
            auto& n1 = nodes_[n1_id];
            assert(n0_id < n1_id);

            std::vector<const MatchInfo*> shared_left;
            std::vector<const MatchInfo*> shared_right;
            for (auto& m0 : end_covs[i0]) {
                for (auto& m1 : end_covs[i1]) {
                    if (m0->GetOverlap()->a_.id == m1->GetOverlap()->a_.id && 
                       (i0 % 2 != i1 % 2 && m0->GetOverlap()->SameDirect() == m1->GetOverlap()->SameDirect() ||
                       (i0 % 2 == i1 % 2 && m0->GetOverlap()->SameDirect() != m1->GetOverlap()->SameDirect()))) {
                        shared_left.push_back(m0);
                        shared_right.push_back(m1);
                        break;
                    }
                }
            }
            assert(shared_left.size() == shared_right.size());
            if (shared_left.size() >= 3) {
                Edge* e = n0->FindEdge(n1);
                if (e != nullptr) {
                    e->AddMatch(shared_left, shared_right );
                } else {
                    auto e_id = MakeEdgeId(n0_id, n1_id);
                    assert(edges_.find(e_id) == edges_.end());

                    Edge* e = new Edge(n0, n1);
                    edges_[e_id] = e;
                    n0->AddEdge(e);
                    n1->AddEdge(e);
                    e->AddMatch(shared_left, shared_right);
                }
            }
        }
    }
}

void ContigGraph::Simplify() {
    // 删除没有边的节点
    // for (const auto& it : nodes_) {
    //     const Node* n = it.second;
    //     if (n->Edges().size() > 1) {

    //         to_delete.insert(n->Id());
    //     }
    // }

}

auto ContigGraph::GetChains() -> std::vector<Chain> {
    std::vector<Chain> chains;

    std::vector<const Node*> path;
    std::unordered_set<const Node*> visited;

    auto extend = [&visited, this](const Node* n) {
        assert(visited.find(n) != visited.end());
        std::vector<const Node*> stack;

        const Node* curr = n;
        while (curr != nullptr && curr->Edges().size() == 1 && visited.find(curr->Edges()[0]->OtherNode(curr)) == visited.end()) {
            
            stack.push_back(curr->Edges()[0]->OtherNode(curr));
            visited.insert(stack.back());

            stack.push_back(GetPairNode(stack.back()));
            visited.insert(stack.back());
            curr = stack.back();
            LOG(INFO)("Extend: %08X %08X", curr->Id(), GetPairNode(curr)->Id());
        }
        return stack;

    };


    // 遍历所有节点，构建链
    for (const auto& it : nodes_) {
        const Node* node = it.second;
        if (visited.find(node) != visited.end()) continue;
        visited.insert(node);
        auto path0 = extend(node);
        path.insert(path.end(), path0.rbegin(), path0.rend());
        path.push_back(node);
        path.push_back(GetPairNode(node));
        visited.insert(path.back());

        auto path1 = extend(path.back());
        path.insert(path.end(), path1.begin(), path1.end());
        // 填充 chain 的数据
        chains.push_back(Chain(std::move(path)));
    }

    return chains;  
}

std::string ContigGraph::Chain::Polish() {
    std::string seq;

    return seq;
}

void ContigGraph::Save(const std::string& fname) const {
    LOG(INFO)("Save graph to %s", fname.c_str());
    std::ofstream of(fname);
    if (of.is_open()) {
        of << "Source,Target,type\n";

        for (const auto& node : nodes_) {
            const Node* n = node.second;
            if (n->Id() % 2 == 0) {
                of << std::uppercase << std::hex << "0x" << n->Id() << ",0x" << GetPairNode(n)->Id() << ",Undirected\n";
            } 
            
            for (const auto& e : n->Edges()) {
                auto o_n = e->OtherNode(n);
                if (o_n->Id() > n->Id()) {
                    of << std::uppercase << std::hex << "0x" << n->Id() << ",0x" << o_n->Id() << ",Undirected\n";
                }
            }
        }
    }

    of.close();
}
}