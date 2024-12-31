#include "sam_tools.hpp"

#include "overlap_store.hpp"
#include "read_store.hpp"
#include "utils/project_file.hpp"
#include "assemble/read_variants.hpp"

#include "overlap/mapping.hpp"

namespace fsa {

void Program_Test::Running() {
    OverlapStore ol_store;
    ol_store.Load(ifname_);
    
}

} // namespace fsa
