#pragma once

#include <vector>

namespace fsa {

class ByteBuffer {
public:
    class View {
    public:
        
        const ByteBuffer* buf { nullptr };
    };
public:
    ByteBuffer() {}
    ~ByteBuffer() {
        
    }

protected:

    std::vector<uint8_t*> data;

};
}