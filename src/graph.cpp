#include "graph.hpp"

#include <map>
#include <numeric>

#include "./utils/logger.hpp"
namespace fsa {


size_t MatrixGraph::Degree(int i) const {
    return std::accumulate(edges_.begin()+i*size_, edges_.begin()+(i+1)*size_, 0);
}

std::vector<int> MatrixGraph::Partition() const {
    std::vector<int> clusters (size_);
    for (size_t i=0; i<clusters.size(); ++i) {
        clusters[i] = i;
    }

    for (size_t i = 0; i<clusters.size(); ++i) {
        for (size_t j=i+1; j<clusters.size(); ++j) {
            if (edges_[i*size_ + j] == 1) {
                if (clusters[j] > clusters[i])
                    clusters[j] = clusters[i];
                else 
                    clusters[i] = clusters[j];
            }
        }
    }

    return clusters;
}

std::vector<std::set<int>> MatrixGraph::Cluster() const {
    
    std::vector<int> types = Partition();

    std::map<int, std::set<int>> cluster;
    for (int i=0; i<(int)types.size(); ++i) {
        cluster[types[i]].insert(i);
    }

    std::vector<std::set<int>> r(cluster.size());
    std::transform(cluster.begin(), cluster.end(), r.begin(), [](const std::pair<int, std::set<int>>& a){ return a.second; });
    return r;
}

void MatrixGraph::Print() const {
    for (size_t i=0;i<size_; ++i) {
        printf("Graph");
        for (size_t j=0; j<size_; ++j) {
            printf(" %d", edges_[i*size_ +j]);
        }
        printf("\n");
    }
    printf("\n");
}
/*

std::vector<int> OverlapFilter::PartitionGraph(size_t sz, const std::vector<uint8_t> & graph) {
    std::vector<int> clusters (sz, 0);
    for (size_t i=0; i<clusters.size(); ++i) {
        clusters[i] = i;
    }

    std::vector<int> k(sz, 0);
    int m = 0;
    for (size_t i=0; i<sz; ++i) {
        for (size_t j=0; j<sz; ++j) {
            if (graph[i*sz+j] > 0) {
                k[i] += 1;
            }
        }
        m += k[i];
    }
    m /= 2;


    auto modularity = [&k, m](size_t sz, const std::vector<uint8_t> & graph, std::vector<int>& clusters) {
        double Q = 0;
        for (size_t i=0; i<sz; ++i) {
            for (size_t j=i+1; j<sz; ++j) {
                if (clusters[i] == clusters[j]) {
                    Q += graph[i*sz+j] - (k[i]*k[j]*1.0)/ (2*m);
                }
            }
        }
        return Q/(2*m);
    };

    auto calc_delta = [&k, m](size_t sz, const std::vector<uint8_t> & graph, std::vector<int>& clusters, int c0, int c1) {
        double d = 0;
        for (size_t i=0; i<sz; ++i) {
            if (clusters[i] == c0 || clusters[i] == c1) {
                for (size_t j=i+1; j<sz; ++j) {
                    if (clusters[j] != clusters[i] &&(clusters[j] == c0 || clusters[j] == c1)) {
                        d += graph[i*sz+j] - (k[i]*k[j]*1.0)/ (2*m);
                    }
                }

            }
        }
        return d/(2*m);
    };
        auto Q = modularity(sz, graph, clusters);
        printf("Q %f\n", Q);

    for (size_t _ = 0; _<clusters.size(); ++_) {
        std::set<int> clu_set(clusters.begin(), clusters.end());
        std::vector<int> clu_vec(clu_set.begin(), clu_set.end());
        
        printf("-----");
        for (auto s : clu_set) printf("%d, ", s);
        printf("\n");
        double max_delta = 0;
        std::array<int,2> max_pair {-1, -1};
        for (size_t i=0; i<clu_vec.size(); ++i) {
            for (size_t j=i+1; j<clu_vec.size(); ++j) {
                double d = calc_delta(sz, graph, clusters, clu_vec[i], clu_vec[j]);
                printf("delta %d %d %f\n", clu_vec[i], clu_vec[j], d);
                if (d > max_delta) {
                    max_delta = d;
                    max_pair = {clu_vec[i], clu_vec[j]};
                }
            }
        }
        if (max_delta > 0) {
            for (size_t i=0; i<clusters.size(); ++i) {
                if (clusters[i] == max_pair[1]) {
                    clusters[i] = max_pair[0];
                }
            }
        }
        auto Q = modularity(sz, graph, clusters);
        printf("Q %f\n", Q);
    }
    
    return clusters;
}
*/
} // namespace fsa {