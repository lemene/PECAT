#pragma

#include <vector>

namespace fsa {

class BaseStore {
public:
    BaseStore(size_t bs= 100*1024*1024) : bsize_(bs) {
        blocks_.push_back(new uint8_t[bs]);
        bindex_ = 0;
    }
    ~BaseStore() {
        for (auto b : blocks_) {
            delete[] b;
        }
        blocks_.clear();
    }
protected:
    std::vector<uint8_t*> blocks_;
    size_t bsize_;
    size_t bindex_;
};

}