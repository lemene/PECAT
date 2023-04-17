#include "diff.hpp"
#include <memory.h>
#include <cassert>
#include <iostream>
#include "utility.hpp"


namespace test {

#define GAP_ALN 4


DiffRunningData::DiffRunningData(int segsize)
 : segment_size_(segsize)
 , U(row_size, 0)
 , V(column_size, 0)
 , align(segment_aln_size)
 , result(max_aln_size)
 , d_path(d_path_size)
 , aln_path(aln_path_size) {
}


DiffRunningData::DPathData2* DiffRunningData::GetDPathIdx(const int d, const int k, const unsigned int max_idx, DPathData2* base)
{
    DPathData2 target;
    target.d = d;
    target.k = k;

    DPathData2* ret = (DPathData2*)bsearch(&target, base, max_idx, sizeof(DPathData2), [](const void* a, const void* b) {
        const DPathData2* d1 = (const DPathData2*)a;
        const DPathData2* d2 = (const DPathData2*)b;
        return (d1->d == d2->d) ? (d1->k - d2->k) : (d1->d - d2->d);
    });
    return ret;
}

bool DiffRunningData::Align(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, const int right_extend, double error_rate)
{
    
    align.init();
    int new_min_k, new_max_k, pre_k;
    int x, y;
    int ck, cd, cx, cy, nx, ny; 
    
    int aln_path_idx, aln_pos, aligned = 0;
    
    int max_d = (int)(2.0 * error_rate * (q_len + t_len));
    int k_offset = max_d;
    int band_size = band_tolerance * 2;
    int best_m = -1;
    int min_k = 0;
    int max_k = 0;
    unsigned long d_path_idx = 0;
    unsigned long max_idx = 0;
    
    for (int d = 0; d < max_d; ++d) {

        if (max_k - min_k > band_size) break;
        
        int  k;
        for (k = min_k; k <= max_k; k += 2) {
            if( k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset]) ) { 
                pre_k = k + 1; x = V[k + 1 + k_offset]; 
            } else { 
                pre_k = k - 1; x = V[k - 1 + k_offset] + 1; 
            }
            y = x - k;
            d_path[d_path_idx].d = d;
            d_path[d_path_idx].k = k;
            d_path[d_path_idx].x1 = x;
            d_path[d_path_idx].y1 = y;
			
			if (right_extend)
				while( x < q_len && y < t_len && query[x] == target[y]) { ++x; ++y; }
			else
				while( x < q_len && y < t_len && query[-x] == target[-y]) { ++x; ++y; }

            d_path[d_path_idx].x2 = x;
            d_path[d_path_idx].y2 = y;
            d_path[d_path_idx].pre_k = pre_k;
            ++d_path_idx;

            V[k + k_offset] = x;
            U[k + k_offset] = x + y;
            best_m = std::max(best_m, x + y);
            if (x >= q_len || y >= t_len)
            { aligned = 1; max_idx = d_path_idx; break; }
        }

        // for banding
        new_min_k = max_k;
        new_max_k = min_k;
        for (int k2 = min_k; k2 <= max_k; k2 += 2) {
            if (U[k2 + k_offset] >= best_m - band_tolerance) { 
                new_min_k = std::min(new_min_k, k2); new_max_k = std::max(new_max_k, k2); 
            }
        }
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;

        if (aligned)
        {
            align.aln_q_e = x;
            align.aln_t_e = y;
            align.dist = d;
            align.aln_str_size = (x + y + d) / 2;
            align.aln_q_s = 0;
            align.aln_t_s = 0;

            if (get_aln_str)
            {
                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
                {
                    DPathData2* d_path_aux = GetDPathIdx(cd, ck, max_idx, &d_path[0]);
                    aln_path[aln_path_idx].x = d_path_aux->x2;
                    aln_path[aln_path_idx].y = d_path_aux->y2;
                    ++aln_path_idx;
                    aln_path[aln_path_idx].x = d_path_aux->x1;
                    aln_path[aln_path_idx].y = d_path_aux->y1;
                    ++aln_path_idx;
                    ck = d_path_aux->pre_k;
                    cd -= 1;
                }
                --aln_path_idx;
                cx = aln_path[aln_path_idx].x;
                cy = aln_path[aln_path_idx].y;
                align.aln_q_s = cx;
                align.aln_t_s = cy;
                aln_pos = 0;
                while (aln_path_idx > 0)
                {
                    --aln_path_idx;
                    nx = aln_path[aln_path_idx].x;
                    ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny) continue;
                    if (cx == nx && cy != ny)
                    {
						if (right_extend)
						{
							for (int i = 0; i < ny - cy; ++i) align.q_aln_str[aln_pos + i] = GAP_ALN;
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[cy + i];
						}
						else
						{
							for (int i = 0; i < ny - cy; ++i) align.q_aln_str[aln_pos + i] = GAP_ALN;
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[-(cy + i)];
						}
                        aln_pos += ny - cy;
                    }
                    else if (cx != nx && cy == ny)
                    {
						if (right_extend)
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[cx + i];
							for (int i = 0; i < nx - cx; ++i) align.t_aln_str[aln_pos + i] = GAP_ALN;
						}
						else
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[-(cx + i)];
							for (int i = 0; i < nx - cx; ++i) align.t_aln_str[aln_pos + i] = GAP_ALN;
						}
                        aln_pos += nx - cx;
                    }
                    else
                    {
						if (right_extend)
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[cx + i];
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[cy + i];
						}
						else
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[-(cx + i)];
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[-(cy + i)];
						}
                        aln_pos += ny - cy;
                    }
                    cx = nx;
                    cy = ny;
                }
                align.aln_str_size = aln_pos;
            }
            break;
        }
    }

    return align.aln_q_e == q_len || align.aln_t_e == t_len;
}

bool DiffRunningData::Align_old(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, const int right_extend, double error_rate)
{
    
    align.init();
    int new_min_k, new_max_k, pre_k;
    int x, y;
    int ck, cd, cx, cy, nx, ny; 
    
    int aln_path_idx, aln_pos, aligned = 0;
    
    int max_d = (int)(2.0 * error_rate * (q_len + t_len));
    int k_offset = max_d;
    int band_size = band_tolerance * 2;
    int best_m = -1;
    int min_k = 0;
    int max_k = 0;
    unsigned long d_path_idx = 0;
    unsigned long max_idx = 0;
    
    for (int d = 0; d < max_d; ++d) {
        
        if (max_k - min_k > band_size) break;
        
        int  k;
        for (k = min_k; k <= max_k; k += 2) {
            if( k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset]) ) { 
                pre_k = k + 1; x = V[k + 1 + k_offset]; 
            } else { 
                pre_k = k - 1; x = V[k - 1 + k_offset] + 1; 
            }
            y = x - k;
            d_path[d_path_idx].d = d;
            d_path[d_path_idx].k = k;
            d_path[d_path_idx].x1 = x;
            d_path[d_path_idx].y1 = y;
			
			if (right_extend)
				while( x < q_len && y < t_len && query[x] == target[y]) { ++x; ++y; }
			else
				while( x < q_len && y < t_len && query[-x] == target[-y]) { ++x; ++y; }

            d_path[d_path_idx].x2 = x;
            d_path[d_path_idx].y2 = y;
            d_path[d_path_idx].pre_k = pre_k;
            ++d_path_idx;

            V[k + k_offset] = x;
            U[k + k_offset] = x + y;
            best_m = std::max(best_m, x + y);
            if (x >= q_len || y >= t_len)
            { aligned = 1; max_idx = d_path_idx; break; }
        }

        // for banding
        new_min_k = max_k;
        new_max_k = min_k;
        for (int k2 = min_k; k2 <= max_k; k2 += 2) {
            if (U[k2 + k_offset] >= best_m - band_tolerance) { 
                new_min_k = std::min(new_min_k, k2); new_max_k = std::max(new_max_k, k2); 
            }
        }
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;

        if (aligned)
        {
            align.aln_q_e = x;
            align.aln_t_e = y;
            align.dist = d;
            align.aln_str_size = (x + y + d) / 2;
            align.aln_q_s = 0;
            align.aln_t_s = 0;

            if (get_aln_str)
            {
                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
                {
                    DPathData2* d_path_aux = GetDPathIdx(cd, ck, max_idx, &d_path[0]);
                    aln_path[aln_path_idx].x = d_path_aux->x2;
                    aln_path[aln_path_idx].y = d_path_aux->y2;
                    ++aln_path_idx;
                    aln_path[aln_path_idx].x = d_path_aux->x1;
                    aln_path[aln_path_idx].y = d_path_aux->y1;
                    ++aln_path_idx;
                    ck = d_path_aux->pre_k;
                    cd -= 1;
                }
                --aln_path_idx;
                cx = aln_path[aln_path_idx].x;
                cy = aln_path[aln_path_idx].y;
                align.aln_q_s = cx;
                align.aln_t_s = cy;
                aln_pos = 0;
                while (aln_path_idx > 0)
                {
                    --aln_path_idx;
                    nx = aln_path[aln_path_idx].x;
                    ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny) continue;
                    if (cx == nx && cy != ny)
                    {
						if (right_extend)
						{
							for (int i = 0; i < ny - cy; ++i) align.q_aln_str[aln_pos + i] = GAP_ALN;
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[cy + i];
						}
						else
						{
							for (int i = 0; i < ny - cy; ++i) align.q_aln_str[aln_pos + i] = GAP_ALN;
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[-(cy + i)];
						}
                        aln_pos += ny - cy;
                    }
                    else if (cx != nx && cy == ny)
                    {
						if (right_extend)
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[cx + i];
							for (int i = 0; i < nx - cx; ++i) align.t_aln_str[aln_pos + i] = GAP_ALN;
						}
						else
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[-(cx + i)];
							for (int i = 0; i < nx - cx; ++i) align.t_aln_str[aln_pos + i] = GAP_ALN;
						}
                        aln_pos += nx - cx;
                    }
                    else
                    {
						if (right_extend)
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[cx + i];
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[cy + i];
						}
						else
						{
							for (int i = 0; i < nx - cx; ++i) align.q_aln_str[aln_pos + i] = query[-(cx + i)];
							for (int i = 0; i < ny - cy; ++i) align.t_aln_str[aln_pos + i] = target[-(cy + i)];
						}
                        aln_pos += ny - cy;
                    }
                    cx = nx;
                    cy = ny;
                }
                align.aln_str_size = aln_pos;
            }
            break;
        }
    }

    return align.aln_q_e == q_len || align.aln_t_e == t_len;
}

void DiffRunningData::dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size,
						const int right_extend, double error_rate)
{
	const idx_t ALN_SIZE = segment_size_;
    int extend1 = 0, extend2 = 0;
    int flag_end = 1;
    
    int extend_size = std::min(query_size, target_size);
    while (flag_end)
    {
        
        int seg_size = 0;
        if (extend_size > (ALN_SIZE + 100)) { 
            seg_size = ALN_SIZE; 
        } else { 
            seg_size = extend_size; flag_end = 0; 
        }
        U.assign(U.size(), 0);
        V.assign(V.size(), 0);

        const char* seq1 = right_extend ? query + extend1 : query - extend1;
        const char* seq2 = right_extend ? target + extend2 : target - extend2;

        bool align_flag = Align(seq1, seg_size, seq2, seg_size, 0.3 * seg_size, 400, right_extend, error_rate);

        if (align_flag)
        {
            int i, j, k, num_matches;
            for (k = align.aln_str_size - 1, i = 0, j = 0, num_matches = 0; k > -1 && num_matches < 4; --k) {
                if (align.q_aln_str[k] != GAP_ALN) ++i;
                if (align.t_aln_str[k] != GAP_ALN) ++j;
                if (align.q_aln_str[k] == align.t_aln_str[k]) ++num_matches;
                else num_matches = 0;
            }

            if (flag_end)
            {
                i = ALN_SIZE - align.aln_q_e + i;
                j = ALN_SIZE - align.aln_t_e + j;
                if (i == ALN_SIZE) align_flag = 0;
				extend1 = extend1 + ALN_SIZE - i; extend2 = extend2 + ALN_SIZE - j;
            }
            else
            {
                i = extend_size - align.aln_q_e;
                j = extend_size - align.aln_t_e;
                if (i == extend_size) align_flag = 0;
				extend1 += (extend_size - i); extend2 += (extend_size - j);
                k = align.aln_str_size - 1;
            }
            if (align_flag)
            {
                if (right_extend) {
                    result.right_store1.insert(result.right_store1.end(), align.q_aln_str.begin(), align.q_aln_str.begin()+k+1);
                    result.right_store2.insert(result.right_store2.end(), align.t_aln_str.begin(), align.t_aln_str.begin()+k+1);
                } else {
                    result.left_store1.insert(result.left_store1.end(), align.q_aln_str.begin(), align.q_aln_str.begin()+k+1);
                    result.left_store2.insert(result.left_store2.end(), align.t_aln_str.begin(), align.t_aln_str.begin()+k+1);
                }
                extend_size = std::min(query_size - extend1, target_size - extend2);
            }
        }
        if (!align_flag) break;
    }
}

fsa::TimeCounter tc("align");
int  DiffRunningData::dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        double error_rate)
{
    fsa::TimeCounter::Mark m(tc);
    result.init();
    
    // left extend
    dw_in_one_direction(query + query_start - 1, query_start,
						target + target_start - 1, target_start,
						0, error_rate);
    // right extend
    dw_in_one_direction(query + query_start, query_size - query_start,
						target + target_start, target_size - target_start,
						1, error_rate);

    // merge the results
    int i, j, k, idx = 0;
    const char* encode2char = "ACGT-";
    for (k = (int)result.left_store1.size() - 1, i = 0, j = 0; k > - 1; --k, ++idx)
    {
		unsigned char ch = result.left_store1[k];
		ch = encode2char[ch];
		result.out_store1[idx] = ch;
		if (ch != '-') ++i;
		
		ch = result.left_store2[k];
		ch = encode2char[ch];
		result.out_store2[idx] = ch;
		if (ch != '-')++j;
    }
    result.query_start = query_start - i;
    result.target_start = target_start - j;

    for (k = 0, i = 0, j = 0; k < (int)result.right_store1.size(); ++k, ++idx)
    {
        char ch = encode2char[(int)result.right_store1[k]];
		result.out_store1[idx] = ch;
		if (ch != '-') ++i;
		
		ch = result.right_store2[k];
		ch = encode2char[(int)ch];
		result.out_store2[idx] = ch;
		if (ch != '-') ++j;
    }
    result.out_store_size = idx;
    result.out_store1[result.out_store_size] = '\0';
    result.out_store2[result.out_store_size] = '\0';
    result.query_end = query_start + i;
    result.target_end = target_start + j;

    StatAlignment();
    tc.Add(result.out_store_size);
    return 1;
}

void DiffRunningData::StatAlignment() {
    for (int j = 0; j < result.out_store_size; ++j)
    {
        if (result.out_store1[j] == result.out_store2[j])
        {
            ++result.mat;
        }
        else if (result.out_store1[j] == '-')
        {
            ++result.ins;
        }
        else if (result.out_store2[j] == '-')
        {
            ++result.del;
        }
        else
        {
            ++result.mis;
        }
    }
}

} // e
