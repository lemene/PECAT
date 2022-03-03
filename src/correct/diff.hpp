#ifndef FSA_DIFF_HPP
#define FSA_DIFF_HPP

#include <algorithm>
#include <stdint.h>
#include <string>
#include <vector>

typedef int64_t        idx_t;



namespace test {

#define MAX_SEQ_SIZE 5000000

struct SW_Parameters {

    idx_t segment_size { 500 }; // small = 500 large = 1000
    idx_t row_size { 4096 };
    idx_t column_size { 4096 };
    idx_t segment_aln_size { 4096 };
    idx_t max_seq_size { MAX_SEQ_SIZE };
    idx_t max_aln_size { MAX_SEQ_SIZE };
    idx_t d_path_size { 5000000 };
    idx_t aln_path_size { 5000000 };
};

struct Alignment
{
    int aln_str_size;
    int dist;
    int aln_q_s;
    int aln_q_e;
    int aln_t_s;
    int aln_t_e;
    std::vector<char> q_aln_str;
    std::vector<char> t_aln_str;
    

    void init()
    {
        aln_str_size = 0;
        aln_q_s = aln_q_e = 0;
        aln_t_s = aln_t_e = 0;
    }

    Alignment(const idx_t max_aln_size) : q_aln_str(max_aln_size), t_aln_str(max_aln_size) { }

};

struct OutputStore {

    std::vector<char> left_store1;
    std::vector<char> left_store2;
    std::vector<char> right_store1;
    std::vector<char> right_store2;
    std::vector<char> out_store1;
    std::vector<char> out_store2;

    int out_store_size;
    int query_start, query_end;
    int target_start, target_end;
    int mat, mis, ins, del;
    
    OutputStore(const idx_t max_aln_size) 
     : left_store1(max_aln_size)
     , left_store2(max_aln_size)
     , right_store1(max_aln_size)
     , right_store2(max_aln_size)
     , out_store1(max_aln_size)
     , out_store2(max_aln_size) {

    }

    void init()
    {
        right_store1.clear();
        right_store2.clear();
        left_store1.clear();
        left_store2.clear();
        out_store_size = 0;
        mat = 0; mis = 0; ins = 0; del =0;
    }
};


struct DiffRunningData {
    
    struct DPathData2 {
        int d, k, pre_k, x1, y1, x2, y2;
    };

    struct PathPoint {
        int x, y;
    };

	DiffRunningData(int seglen);

    int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        double error_rate);

    bool Align(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, 
		  const int right_extend, double error_rate);

          
    bool Align_old(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, 
		  const int right_extend, double error_rate);


    void dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size, const int right_extend, double error_rate);


    static DPathData2* GetDPathIdx(const int d, const int k, const unsigned int max_idx, DPathData2* base);
    void StatAlignment();


    idx_t segment_size_ { 500 }; // small = 500 large = 1000
    idx_t row_size { 4096 };
    idx_t column_size { 4096 };
    idx_t segment_aln_size { 4096 };
    idx_t max_seq_size { MAX_SEQ_SIZE };
    idx_t max_aln_size { MAX_SEQ_SIZE };
    idx_t d_path_size { 5000000 };
    idx_t aln_path_size { 5000000 };

    //SW_Parameters   swp;
    std::vector<int> U; //DynQ;
    std::vector<int> V; //DynT;
    Alignment      align;
    OutputStore    result;
    std::vector<DPathData2> d_path;
    std::vector<PathPoint> aln_path;
};


inline int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        DiffRunningData &drd, double error_rate) {

    return drd.dw(query, query_size, query_start, target, target_size, target_start, 
			error_rate);

}

} // 

#endif  // FSA_DIFF_HPP
