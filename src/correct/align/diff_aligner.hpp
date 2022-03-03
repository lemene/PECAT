#ifndef FSA_CORRECT_DIFF_ALIGNER_HPP
#define FSA_CORRECT_DIFF_ALIGNER_HPP

#include "tool_aligner.hpp"
#include "../diff.hpp"

namespace fsa {
class DiffAligner : public ToolAligner {
public:
    DiffAligner(const std::vector<std::string> &paramers);
    virtual ~DiffAligner();
    
    virtual bool Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t,2> qrange, std::array<size_t, 2> trange, Alignment& al);
protected:
    double error_rate_ { 0.15 };
    int segment_size_ { 500 };
    test::DiffRunningData* drd;
};

} // namespace fsa

#endif  // FSA_CORRECT_DIFF_ALIGNER_HPP