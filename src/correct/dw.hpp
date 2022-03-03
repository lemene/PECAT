#ifndef DW_H
#define DW_H

#include <algorithm>
#include <stdint.h>

typedef int8_t      i1_t;
typedef uint8_t     u1_t;
typedef int16_t     i2_t;
typedef uint8_t     u2_t;
typedef int32_t     i4_t;
typedef uint32_t    u4_t;
typedef int64_t     i8_t;
typedef uint64_t    u8_t;
typedef i8_t        idx_t;

typedef u1_t		uint1;
typedef idx_t		index_t;

#define safe_malloc(arr, type, count) \
do { \
    size_t __sm__sz__ = sizeof(type) * (count); \
    (arr) = (type *)malloc(__sm__sz__); \
} while(0)

#define safe_free(arr) free(arr)

//#include "../common/alignment.h"
//#include "../common/defs.h"
//#include "../common/packed_db.h"

namespace mecat {

struct SW_Parameters
{
    idx_t segment_size;
    idx_t row_size;
    idx_t column_size;
    idx_t segment_aln_size;
    idx_t max_seq_size;
    idx_t max_aln_size;
    idx_t d_path_size;
    idx_t aln_path_size;
};

SW_Parameters
get_sw_parameters_small();

SW_Parameters
get_sw_parameters_large();

struct Alignment
{
    int aln_str_size;
    int dist;
    int aln_q_s;
    int aln_q_e;
    int aln_t_s;
    int aln_t_e;
    char* q_aln_str;
    char* t_aln_str;

    void init()
    {
        aln_str_size = 0;
        aln_q_s = aln_q_e = 0;
        aln_t_s = aln_t_e = 0;
    }

    Alignment(const idx_t max_aln_size)
    {
        safe_malloc(q_aln_str, char, max_aln_size);
        safe_malloc(t_aln_str, char, max_aln_size);
    }
    ~Alignment()
    {
        safe_free(q_aln_str);
        safe_free(t_aln_str);
    }
};

struct OutputStore
{
    char* left_store1;
    char* left_store2;
    char* right_store1;
    char* right_store2;
    char* out_store1;
    char* out_store2;
    char* out_match_pattern;

    int left_store_size;
    int right_store_size;
    int out_store_size;
    int query_start, query_end;
    int target_start, target_end;
    int mat, mis, ins, del;
	double ident;
    
    OutputStore(const idx_t max_aln_size)
    {
        safe_malloc(left_store1, char, max_aln_size);
        safe_malloc(left_store2, char, max_aln_size);
        safe_malloc(right_store1, char, max_aln_size);
        safe_malloc(right_store2, char, max_aln_size);
        safe_malloc(out_store1, char, max_aln_size);
        safe_malloc(out_store2, char, max_aln_size);
        safe_malloc(out_match_pattern, char, max_aln_size);
    }

    ~OutputStore()
    {
        safe_free(left_store1);
        safe_free(left_store2);
        safe_free(right_store1);
        safe_free(right_store2);
        safe_free(out_store1);
        safe_free(out_store2);
        safe_free(out_match_pattern);
    }

    void init()
    {
        left_store_size = right_store_size = out_store_size = 0;
    }
};

struct DPathData
{
    int pre_k, x1, y1, x2, y2;
};

struct DPathData2
{
    int d, k, pre_k, x1, y1, x2, y2;
};

struct PathPoint
{
    int x, y;
};

struct DiffRunningData
{
    SW_Parameters   swp;
    char*           query;
    char*           target;
    int*            DynQ;
    int*            DynT;
    Alignment*      align;
    OutputStore*    result;
    DPathData2*     d_path;
    PathPoint*      aln_path;
	
	DiffRunningData(const SW_Parameters& swp_in);
	~DiffRunningData();
};



struct CandidateStartPosition
{
	idx_t qoff;
	idx_t toff;
	idx_t tstart;
	idx_t tsize;
	idx_t tid;
	int left_q, left_t;
	int right_q, right_t;
	int num1, num2;
	int score;
	idx_t toff_in_aln;
	char chain;
};


int Align(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, Alignment* align, 
		  int* V, int* U, DPathData2* d_path, PathPoint* aln_path, const int right_extend);

void dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size,
						 int* U, int* V, Alignment* align, DPathData2* d_path, PathPoint* aln_path, 
						 SW_Parameters* swp, OutputStore* result, const int right_extend);

int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        int* U, int* V, Alignment* align, DPathData2* d_path, 
        PathPoint* aln_path, OutputStore* result, SW_Parameters* swp,
	    double error_rate, const int min_aln_size);


inline int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        DiffRunningData &drd, double error_rate, const int min_aln_size) {

    return dw(query, query_size, query_start, target, target_size, target_start, 
			  drd.DynQ, drd.DynT, drd.align, drd.d_path, drd.aln_path, drd.result, &drd.swp, 
              error_rate, min_aln_size);

}

} // 

#endif  // DW_H
