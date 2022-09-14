#ifndef FSA_CORRECT_EDLIB_ALIGNER_HPP
#define FSA_CORRECT_EDLIB_ALIGNER_HPP

#include "tool_aligner.hpp"
namespace fsa {

class EdlibAligner : public ToolAligner {
public:
    EdlibAligner(const std::vector<std::string> &parameters);
    virtual ~EdlibAligner();
    
    virtual bool Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t,2> qrange, std::array<size_t, 2> trange, Alignment& al);
    virtual bool Align1(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t,2> qrange, std::array<size_t, 2> trange, Alignment& al);

protected:
    struct AlignResult {
        int distance {0};
        size_t qstart {0};
        size_t qend {0};
        size_t tstart {0};
        size_t tend {0};
        std::vector<unsigned char> alignment;
    };
    
    void AlignedString(const unsigned char* alignment, int alignmentLenght, const char* query, const char* target, std::string& aligned_query, std::string& aligned_target);
    bool ExtendRight(const char* qseq, size_t qsize, const char* tseq, size_t tsize, size_t qstart, size_t tstart, AlignResult &result);
protected:
    size_t block_size_ { 1000 };
    size_t match_count_ { 6 };
    double block_error_ { 0.50 };
};
















} // namespace fsa {

#endif // FSA_CORRECT_EDLIB_ALIGNER_HPP