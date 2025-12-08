#include "utility.hpp"

#include <algorithm>
#include <unistd.h>
#include <string.h>

namespace fsa {

std::vector<std::string> SplitStringBySpace(const std::string &str) {
    std::vector<std::string> substrs;

    auto IsSpace = [](char c) { return ::isspace(c); };
    auto IsNotSpace = [](char c) { return !::isspace(c); };

    auto begin = std::find_if(str.begin(), str.end(), IsNotSpace);

    while (begin != str.end()) {
        auto end = std::find_if(begin, str.end(), IsSpace);
        substrs.push_back(std::string(begin, end));
        begin = std::find_if(end, str.end(), IsNotSpace);
    }

    return substrs;
}



size_t FindLongestXHeap(std::vector<int> &lengths, long long base_size) {

    assert(base_size > 0);

    long long accu = 0;
    size_t index = 0;
    auto cmp = [](size_t a, size_t b) { return a > b; };
    for (size_t i = 0; i< lengths.size(); ++i) {
        if (accu < base_size || lengths[0] < lengths[i]) {
            accu += lengths[i];
            std::swap(lengths[index], lengths[i]);
            index++;
            std::push_heap(lengths.begin(), lengths.begin()+index, cmp);
        }

        while (accu - lengths[0] >= base_size) {
            accu -= lengths[0];
            std::pop_heap(lengths.begin(), lengths.begin()+index, cmp);
            index--;
        }
    }
    assert(index > 0);
    
    return index;
}

size_t FindLongestXSort(std::vector<int> &lengths, long long base_size) {

    long long accu = 0;

    std::sort(lengths.begin(), lengths.end(), [](int a, int b) { return a > b; });

    for (size_t i = 0; i< lengths.size(); ++i) {
        if (accu < base_size) {
            accu += lengths[i];
        } else {
            return i;
        }
    }

    return lengths.size();
}

size_t GetMemoryUsage() {
    
    int vmrss_num = 0;
    int pid = getpid();
    char fname[48];
    sprintf(fname, "/proc/%d/status", pid);
    FILE *file = fopen(fname, "r");
    if (file != NULL) {
        const int BUF_SIZE = 1024;
        char buf[BUF_SIZE];
        char* line = fgets(buf, BUF_SIZE, file);
        while (line != NULL) {
            if(strstr(line, "VmRSS:") != NULL) { // KB
                char vmrss_name[48];
                sscanf(line, "%s %d", vmrss_name, &vmrss_num);
                break;
            }
            line = fgets(buf, BUF_SIZE, file);
        }
        fclose(file);
    }
    return vmrss_num;
}


size_t CountLinesInFile(std::ifstream& ifs, size_t threads) {
    ifs.seekg(0);
    std::string line;
    std::mutex mutex;
    std::atomic<size_t> count { 0 };

    auto generate_func = [&ifs, &mutex](char* buf, size_t bufsize) -> size_t {
        std::lock_guard<std::mutex> lock(mutex);

        ifs.read(buf, bufsize);
        return ifs.gcount();
    };
    auto worker_func = [&count, &generate_func](size_t id) {
        char buf[1024*1000];
        size_t bsize = generate_func(buf, sizeof(buf));
        size_t cnt = 0;
        while (bsize > 0) {
            for (char* p = buf; p < buf+bsize; ++p) {
                if (*p == '\n') cnt ++;
            }
            bsize = generate_func(buf, sizeof(buf));
        }
        count.fetch_add(cnt);
    };

    MultiThreadRun(threads, worker_func);
    return count.load();
}

size_t CountLinesInFile(const std::string& fnane, size_t threads) {
    std::ifstream ifs(fnane);
    return CountLinesInFile(ifs, threads);
}

} // namespace fsa {
    