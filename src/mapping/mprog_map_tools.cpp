#include "mprog_map_tools.hpp"

#include "overlap_store.hpp"
#include "read_store.hpp"

#include "overlap/mapping.hpp"
namespace fsa {

void Program_NextBaseAccuracy::Running() {
    std::mutex mutex;
    
    std::array<size_t, 32> total = {0};
    std::array<size_t, 32> corr = {0};

    TraversePafFile(ifname_, [&mutex, &total, &corr](const Overlap& ol) {
        std::array<size_t, 32> ttt = {0};
        std::array<size_t, 32> ccc = {0};
        for (auto& d : ol.detail_) {
            if (d.type == '=') {
                for (int i = 0; i < ttt.size(); ++i) {
                    if (i + 1 <= d.len) {
                        ttt[i] += d.len - (i+1) + 1;
                        ccc[i] += d.len - (i+1);
                    } else {
                        break;
                    }
                }
                
            }
        }
        {
            std::lock_guard<std::mutex> lock(mutex);
            for (int j = 0; j < ttt.size(); ++j) {
                total[j] += ttt[j];
                corr[j] += ccc[j];
            }
        }
    }, (size_t)thread_size_);
    for (size_t i = 0; i < total.size(); ++i) {
        printf("first_error %zd %zd %zd %.06f\n", i+1, corr[i], total[i], total[i] > 0 ? corr[i]*1.0/total[i] : 0.0);
    }
}

} // namespace fsa
