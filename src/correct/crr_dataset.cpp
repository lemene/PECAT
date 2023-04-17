#include "crr_dataset.hpp"

#include "crr_options.hpp"


namespace fsa {

void CrrDataset::Load(const StringPool &sp) {
    
    if (!opts_.variants.empty()) {

        variants.reset(new SnpFile(sp));
        variants->Load(opts_.variants);
    }
}


}   // namespace fsa