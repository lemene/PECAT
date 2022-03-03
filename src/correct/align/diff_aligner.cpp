#include "diff_aligner.hpp"

#include "../../utility.hpp"

namespace fsa {

DiffAligner::DiffAligner(const std::vector<std::string>& paramters) {

    for (size_t i=1; i<paramters.size(); ++i) {
        auto kv = SplitStringByChar(paramters[i], '=');
        if (kv[0] == "s") {
            segment_size_ = std::stoi(kv[1]);
        } else if (kv[0] == "e") {
            error_rate_ = std::stod(kv[1]);
        }
    }
    drd = new test::DiffRunningData(segment_size_);
}

DiffAligner::~DiffAligner() {
    delete drd;
}

bool DiffAligner::Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t, 2> qrange, std::array<size_t, 2> trange, Alignment& al) {
    auto r = test::dw(qseq, qsize, qrange[0], tseq, tsize, trange[0], *drd, error_rate_);
    if (r) {
        al.target_start = drd->result.target_start;
        al.target_end = drd->result.target_end;
        al.query_start = drd->result.query_start;
        al.query_end = drd->result.query_end;
        al.aligned_query = &drd->result.out_store1[0];
        al.aligned_target = &drd->result.out_store2[0];
        al.distance = drd->result.mis*2 + drd->result.ins + drd->result.del;
        return true;
    } else {
        return false;
    }

}


// bool DiffAligner::Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t, 2> qrange, std::array<size_t, 2> trange, Alignment& al) {
//     auto r = test::dw(qseq, qsize, qrange[0], tseq, tsize, trange[0], *drd, 0.15, 1000);
//     if (r) {
        
//         al.target_start = drd->result.target_start;
//         al.target_end = drd->result.target_end;
//         al.query_start = drd->result.query_start;
//         al.query_end = drd->result.query_end;
//         al.aligned_query = drd->result.out_store1;
//         al.aligned_target = drd->result.out_store2;
//         al.distance = drd->result.mis*2 + drd->result.ins + drd->result.del;
//         return true;
//     } else {
//         return false;
//     }

// }


}