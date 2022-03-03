#ifndef FSA_UTILITY_HPP
#define FSA_UTILITY_HPP

#include <array>
#include <vector>
#include <thread>
#include <future>
#include <atomic>
#include <unordered_set>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <cmath>

namespace fsa {

template<typename T>
auto SplitConstIterater(size_t sz, const T& container) -> std::vector<std::array<typename T::const_iterator, 2>> {
    assert(sz >= 1);

    std::vector<std::array<typename T::const_iterator, 2>> result(sz);

    size_t sub_size = (container.size() + sz - 1) / sz;
    size_t index = 0;
    for (typename T::const_iterator it = container.begin(); it != container.end(); ++it, ++index) {
        if (index % sub_size == 0) {
            result[index/sub_size][0] = it;
        }
    }

    for (size_t i=0; i<result.size()-1; ++i) {
        result[i][1] = result[i+1][0];
    }
    result.back()[1] = container.end();

    return result;
}

template<typename M>
auto SplitMapKeys(size_t sz, const M& map) -> std::vector<std::vector<typename M::key_type>> {
    std::vector<std::vector<typename M::key_type>> results(sz);
    size_t sub_size = (map.size() + sz - 1) / sz;
    typename M::size_type index = 0;
    for (const auto &i : map) {
        results[index++ / sub_size].push_back(i.first);
    }
    return results;
}

template<typename S>
auto SplitSet(size_t sz, const S& set) -> std::vector<S> {
    std::vector<S> results(sz);

    size_t sub_size = (set.size() + sz - 1) / sz;
    typename S::size_type index = 0;
    for (const auto &i : set) {
        results[index++ / sub_size].insert(i);
    }
    return results;
}

template<typename V>
auto SplitVectorKeys(size_t sz, const V& container) ->std::vector<std::array<size_t, 2>> {
    std::vector<std::array<size_t, 2>> result;
    
    size_t bsize = container.size() / sz;
    for (size_t i = 1; i < sz; ++i) {
        result.push_back(std::array<size_t, 2>({ (i - 1)*bsize, i*bsize }));
    }
    result.push_back(std::array<size_t, 2>({ (sz-1)*bsize, container.size() }));

    return result;
}


template<typename T=int>
auto SplitRange(size_t sz, T low, T high) -> std::vector<std::array<T,2>>{
    assert(sz >= 1 && high >= low);
    std::vector<std::array<T,2>> result(sz);
    T inv = (high - low) / sz;

    for (size_t i=0; i<sz-1; ++i) {
        result[i][0] = low + inv*i;
        result[i][1] = low + inv*(i+1);
    }
    result[sz-1][0] = low + inv*(sz-1);
    result[sz-1][1] = high;
    return result;
}

template<typename T>
auto MoveCombineVector(const std::vector<T> &sub_output) -> T {
    T output;
    for (auto &so : sub_output) {
        // TODO output.insert(output.end(), std::make_move_iterator(so.begin()), std::make_move_iterator(so.end()));
        output.insert(output.end(), so.begin(), so.end());
    }
    return output;
}
template<typename T>
auto MoveCombineMapOrSet(const std::vector<T> &sub_output) -> T {
    T output;
    for (auto &so : sub_output) {
        // TODO output.insert(output.end(), std::make_move_iterator(so.begin()), std::make_move_iterator(so.end()));
        output.insert(so.begin(), so.end());
    }
    return output;
}

template<typename W>
void MultiThreadRun(size_t thread_size, W work_func) {
    std::vector<std::thread> workers;
    for (size_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(work_func, i));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
}


template<typename W, typename C>
void AutoThreadRun(W work_func, C check_new) {
    std::vector<std::thread> workers;

    std::future<void> main = std::async(std::launch::async, work_func, 0);
    while (main.wait_for(std::chrono::milliseconds(1)) == std::future_status::timeout) {
        if (check_new()) {
            workers.push_back(std::thread(work_func, workers.size()+1));
        }
    }
    

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
}

template<typename S, typename W>
void MultiThreadRun(size_t thread_size, S split_func, W work_func) {
    auto sub_inputs = split_func();
    assert(sub_inputs.size() == thread_size);

    std::vector<std::thread> workers;
    for (size_t i = 0; i < sub_inputs.size(); ++i) {
        workers.push_back(std::thread(work_func, std::ref(sub_inputs[i])));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }

}


// TODO simplify return type
// 
template<typename I, typename S, typename W, typename C>
auto MultiThreadRun(size_t thread_size, const I &inputs, S split_func, W work_func, C combine_func)
  -> decltype(combine_func(std::vector<decltype(work_func(split_func(thread_size, inputs)[0]))>(1, work_func(split_func(thread_size, inputs)[0])))) {
    auto sub_inputs = split_func(thread_size, inputs);
    assert(sub_inputs.size() == thread_size);
    std::vector<decltype(work_func(sub_inputs[0]))> sub_outputs(thread_size);

    // TODO "const decltype(sub_inputs[0])& in" get compiler error
    auto thread_work = [&work_func](decltype(sub_inputs[0])& in, decltype(sub_outputs[0])& out) {
        out = work_func(in);
    };

    std::vector<std::thread> workers;
    for (size_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(thread_work, std::ref(sub_inputs[i]), std::ref(sub_outputs[i])));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }

    return combine_func(sub_outputs);
}

template<typename I, typename S, typename W>
void MultiThreadRun(size_t thread_size, I &inputs, S split_func, W work_func) {
    auto sub_inputs = split_func(thread_size, inputs);
    assert(sub_inputs.size() == thread_size);

    std::vector<std::thread> workers;
    for (size_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(work_func, std::ref(sub_inputs[i])));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
}

template<typename I, typename O, typename W>
void MultiThreadMap(size_t thread_size, const std::vector<I>& input, std::vector<O>& output, W map_func) {
    assert(input.size() == output.size());

    std::atomic<size_t> index { 0 };
    auto work_func = [&index, &input, &output, map_func](int tid) {
        for (size_t curr = index.fetch_add(1); curr < input.size(); curr = index.fetch_add(1)) {
            output[curr] = map_func(input[curr]);
        }
    };
    MultiThreadRun(thread_size, work_func);
}


template<typename T, size_t N>
struct ArrayHash {
    size_t operator()(const std::array<T, N> & r) const {
        size_t h = 17;
        for (size_t i = 0; i < N; ++i) {
            // h ^= std::hash<T>()(r[i]*(T)(i+1)* (T)122777);   // TODO test the efficiency of number
            h = h * 31 + std::hash<T>()(r[i]);
        }
        return h;
    }
};


template<typename T, size_t N>
struct ArrayEqual {
    size_t operator()(const std::array<T, N> &a, const std::array<T, N> & b) const {
        for (size_t i = 0; i < N; ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }
};


template<int N, typename T>
int CompareArray(const std::array<T, N>& a, const std::array<T, N> &b) {
    for (size_t i = 0; i < N; ++i) {
        if (a[i] < b[i]) return -1;
        else if (a[i] > b[i]) return 1;
    }
    return 0;
}

template<size_t N>
struct ArrayHash<std::string, N> {
    size_t operator()(const std::array<std::string, N> & r) const {
        size_t h = 17;
        for (size_t i = 0; i < N; ++i) {
            //h ^= std::hash<T>()(r[i]*(T)(i+1)* (T)122777);   // TODO test the efficiency of number
            h = h * 31 + std::hash<std::string>()(r[i]);
        }
        return h;
    }
};



template<typename M>
auto FindIntersectKeys(const M &a, const M &b) -> std::unordered_set<typename M::key_type> {
    std::unordered_set<typename M::key_type> result;

    for (const auto &i : a) {
        if (b.find(i.first) != b.end()) result.insert(i.first);
    }

    return result;
}


std::vector<std::string> SplitStringBySpace(const std::string &str);

template<typename C>
std::vector<std::string> SplitString(const std::string &str, C check) {
    auto is_not_split_point = [check](char c) { return !check(c); };
    auto is_split_point = [check](char c) { return check(c); };

    std::vector<std::string> substrs;
    auto begin = std::find_if(str.begin(), str.end(), is_not_split_point);

    while (begin != str.end()) {
        auto end = std::find_if(begin, str.end(), is_split_point);
        substrs.push_back(std::string(begin, end));
        begin = std::find_if(end, str.end(), is_not_split_point);
    }

    return substrs;
}

inline std::vector<std::string> SplitStringByChar(const std::string &str, char c) {
    return SplitString(str, [c](char i) { return i == c; });
}

inline std::string JoinStrings(const std::vector<std::string> &strs, const std::string& sep) {
    std::string all = strs.size() > 0 ? strs[0] : "";
    for (size_t i = 1; i < strs.size(); ++i) {
        all += sep;
        all += strs[i];
    }
    return all;
}

template<typename T>
void DeletePtrContainer(T & c) {
    for (auto e : c) {
        delete e;
    }
    c.clear();
}

size_t FindLongestXHeap(std::vector<int> &lengths, long long base_size);
size_t FindLongestXSort(std::vector<int> &lengths, long long base_size);

template<typename T>
void ComputeMedianAbsoluteDeviation(const std::vector<T>& data_, T &median, T &mad) {
    std::vector<T> data = data_;
    assert(data.size() > 0);
    std::nth_element(data.begin(), data.begin() + data.size()/2, data.end());
    median = data[data.size()/2];

    for (size_t i=0; i<data.size(); ++i) {
        data[i] = std::abs(data[i]-median);
    }
    std::nth_element(data.begin(), data.begin() + data.size()/2, data.end());
    mad = data[data.size()/2];
}

template<typename T>
std::array<T, 2> ComputeMedianAbsoluteDeviation(const std::vector<T>& data_) {
    std::vector<T> data = data_;
    assert(data.size() > 0);
    std::nth_element(data.begin(), data.begin() + data.size()/2, data.end());
    T median = data[data.size()/2];

    for (size_t i=0; i<data.size(); ++i) {
        data[i] = std::abs(data[i]-median);
    }
    std::nth_element(data.begin(), data.begin() + data.size()/2, data.end());
    T mad = data[data.size()/2];
    return {median, mad};
}

template<typename T>
void ComputeMedianAbsoluteDeviation(const std::vector<std::array<T,2>>& data_, T &median, T &mad) {
    std::vector<std::array<T,2>> data = data_;
    assert(data.size() > 0);

    auto find_median = [](std::vector<std::array<T,2>> &data) {
        std::sort(data.begin(), data.end(), [](const std::array<T,2> &a, const std::array<T,2> &b) {
            return a[0] < b[0];
        });

        T total = std::accumulate(data.begin(), data.end(), 0, [](T a, const std::array<T,2> &b) {
            return a + b[1];
        });
        
        T accu = 0;
        size_t index = 0;
        for (; index < data.size(); ++index) {
             accu += data[index][1];
             if (accu >= total / 2) break;
        }

        return data[index][0];
    };

    median = find_median(data);
    for (size_t i=0; i<data.size(); ++i) {
        data[i][0] = std::abs(data[i][0]-median);
    }

    mad = find_median(data);
}



template<typename T>
void ComputeMeanAbsoluteDeviation(std::vector<T>& data, T &mean, T &mad) {
    assert(data.size() > 0);
    mean = std::accumulate(data.begin(), data.end(), (T)0) / data.size();

    for (size_t i=0; i<data.size(); ++i) {
        data[i] = std::abs(data[i]-mean);
    }
    
    mad = std::accumulate(data.begin(), data.end(), (T)0) / data.size();;
}


template<typename T>
void ComputeMeanAbsoluteDeviation(std::vector<std::array<T,2>>& data, T &mean, T &mad) {
    assert(data.size() > 0);

    auto find_mean = [](std::vector<std::array<T,2>> &data) {
        std::sort(data.begin(), data.end(), [](const std::array<T,2> &a, const std::array<T,2> &b) {
            return a[0] < b[0];
        });

        T total = std::accumulate(data.begin(), data.end(), (T)0, [](T a, const std::array<T,2> &b) {
            return a + b[0]*b[1];
        });
        
        T weight = std::accumulate(data.begin(), data.end(), (T)0, [](T a, const std::array<T,2> &b) {
            return a + b[1];
        }); 

        return total / weight;
    };

    mean = find_mean(data);
    for (size_t i=0; i<data.size(); ++i) {
        data[i][0] = std::abs(data[i][0]-mean);
    }

    mad = find_mean(data);
}

struct TimeCounter {
    TimeCounter(const std::string& n) : name(n) {}
    ~TimeCounter()  { 
        std::cerr << "Record Time(" << name << "):" << count << ", " << ticks / 10  << ", " << value << "\n";
    }

    void Inc(double s) {
        ticks.fetch_add(int(s*10000));
        count.fetch_add(1);
    }

    void Add(long long v) { value.fetch_add(v); }
    struct Mark {
        Mark(TimeCounter& rt) : recTime(rt) { 
            //start = clock();
            clock_gettime(CLOCK_MONOTONIC, &start);
        }
        ~Mark() { 
            
            clock_gettime(CLOCK_MONOTONIC, &finish);
            
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            recTime.Inc(elapsed); 
        }
        TimeCounter& recTime;
        
        //clock_t start;
        //struct timespec start;
        struct timespec start, finish;
        double elapsed;

    };

    std::atomic<long long> ticks {0};
    std::atomic<long long> count {0};
    std::atomic<long long> value {0};
    std::string name;
};


template <typename T, typename U>
bool HasCommon(const T &a, const U &b) 
{
    for (const auto &i : a) {
        if (std::find(b.begin(), b.end(), i) != b.end()) return true;
    }
    return false;
}

template<typename T>
void Intersection(const std::unordered_set<T> &a, std::unordered_set<T> &b, std::unordered_set<T> &c) {

    for (const auto &i : a) {
        if (b.find(i) != b.end())
            c.insert(i);
    }
}


} // namespace fsa {

#endif // FSA_UTILITY_HPP  
