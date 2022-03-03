#pragma once

#include <array>
#include <string>
#include <memory>

#include "alignment.hpp"

namespace fsa {

class ToolAligner {
public:
    virtual ~ToolAligner() {}
    virtual void SetParameter(const std::string& name, const std::string& value="") {}
    virtual bool Align(const char* qseq, size_t qsize, const char* tseq, size_t tsize, std::array<size_t,2> qrange, std::array<size_t, 2> trange, Alignment& al) = 0;
    
    static std::shared_ptr<ToolAligner> Create(const std::string &opts);
};

} // namespace fsa
