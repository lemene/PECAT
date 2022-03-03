#include "tool_aligner.hpp"
#include "../../utility.hpp"
#include "../..//utils/logger.hpp"


#include "./diff_aligner.hpp"
#include "./edlib_aligner.hpp"
#include "./ksw2_aligner.hpp"

namespace fsa {

std::shared_ptr<ToolAligner> ToolAligner::Create(const std::string &opts) {
    std::shared_ptr<ToolAligner> aligner;
    auto ss = SplitStringByChar(opts, ':');
    if (ss.size() >= 1) {
        if (ss[0] == "diff") {
            aligner.reset(new DiffAligner(ss));
        } else if (ss[0] == "edlib") {
            aligner.reset(new EdlibAligner(ss));
        } else if (ss[0] == "ksw2") {
            aligner.reset(new Ksw2Aligner(ss));
        } else {
            LOG(ERROR)("Not support parameter: aligner=%s", opts.c_str());
        }
    } else {
        LOG(ERROR)("Not support parameter: aligner=%s", opts.c_str());

    }

    return aligner;
}

} // namespace fsa {
