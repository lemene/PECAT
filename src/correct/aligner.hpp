#ifndef FSA_ALIGNER_HPP
#define FSA_ALIGNER_HPP


#include <array>
#include <string>
#include <atomic>
#include <iostream>
#include <memory>

#include "sequence.hpp"
#include "simple_align.hpp"
#include "utility.hpp"
#include "./utils/logger.hpp"

#include "align/tool_aligner.hpp"
#include "align/alignment.hpp"

namespace fsa {


class Aligner {
public:
    Aligner();
    ~Aligner();

    void SetParameter(const std::string& name, double v);
    void SetParameter(const std::string& name, const std::string &v);

    void SetTarget(const DnaSeq& tseq);
    void SetTarget(const DnaSeq& tseq, const std::array<size_t, 2> trange);

    bool Align(const DnaSeq &qseq, bool rc /* reverse complement */, const std::array<int, 4> &range, Alignment &al);
    std::array<int, 2> FindExactMatch(const std::vector<uint8_t>& tseq, const std::vector<uint8_t>& qseq, const std::array<int , 2> &s);

    static void AppendAlignedString(const uint32_t * cigar, size_t cigarLen, const char* query, const char* target, std::string& aligned_query, std::string& aligned_target);

    static std::array<double,2> ComputeIdentity(const std::string& alq, const std::string& alt, size_t window_size);
    static bool CheckAlignedString(const std::string &q, const std::string &t);

protected:
    void SetAligner(const std::string &opts);
protected:
    const DnaSeq* target_;
    std::vector<uint8_t> target_0123_;
    std::array<size_t, 2> trange_ {{0, 0}};
    double min_identity_ { 70 };
    double min_local_identity_ { 50 };
    double local_window_size_ { 5000 };

    std::shared_ptr<ToolAligner> worker;
};


} // namespace fsa {

#endif // FSA_ALIGNER_HPP
