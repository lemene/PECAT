#include "contig_graph.hpp"

#include <unordered_set>

namespace fsa {

ContigGraph::~ContigGraph() {
    for (auto it : nodes_) {
        delete it.second;
    }
    nodes_.clear();
}

void ContigGraph::Build() {
    // 创建图的节点
    for (size_t i = 0; i < fragments_.size(); ++i) {
        auto n0_id = MakeNodeId(i, 0);
        auto n1_id = MakeNodeId(i, 1);
        nodes_[n0_id] = new Node(n0_id, &fragments_[i], 0);
        nodes_[n1_id] = new Node(n1_id, &fragments_[i], 1);
    }

    // 创建默认边
    for (size_t i = 1; i < fragments_.size(); ++i) {
        auto& frg0 = fragments_[i-1];
        auto& frg1 = fragments_[i];

        if (frg0.ContigId() == frg1.ContigId()) {
            auto n0_id = MakeNodeId(i, 1);
            auto n1_id = MakeNodeId(i+1, 0);
            nodes_[n0_id]->AddLink(nodes_[n1_id]);
            nodes_[n1_id]->AddLink(nodes_[n0_id]);
        }
    }

    // 检查节点是否通过读数直接连接


}

auto ContigGraph::GetChains() -> std::vector<Chain> {
    std::vector<Chain> chains;

    std::vector<const Node*> path;
    std::unordered_set<const Node*> visited;

    auto extend = [&visited, this](const Node* n) {
        assert(visited.find(n) != visited.end());
        std::vector<const Node*> stack;

        const Node* curr = n;
        while (curr != nullptr && visited.find(curr) == visited.end() && !curr->Links().size() == 1) {
            stack.push_back(curr->Links()[0].To());
            visited.insert(stack.back());

            stack.push_back(GetPairNode(stack.back()));
            visited.insert(stack.back());
            curr = stack.back();
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

        auto path1 = extend(GetPairNode(path.back()));
        path.insert(path.end(), path1.begin(), path1.end());

        // 填充 chain 的数据
        chains.push_back(Chain(std::move(path)));
    }

    return chains;  
}

std::string ContigGraph::Chain::Polish() {
    for (auto& n : paths_) {
        printf("Node ID: %u, Fragment: %u\n", n->Id(), n->Fragment()->ContigId());
    }
    return "";
}

}