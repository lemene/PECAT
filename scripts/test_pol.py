#!/usr/bin/env python3

"""暂时存放测试polish工程的脚本"""

import sys
sys.path.insert(0, "/home/niefan/work/mbio")
import mbio
import mbio.ftype.table as tb
import mbio.utils.utils as utils

def ts_show_win_cov(argv):

    idx = int(argv[0])

    covs = []
    for line in open(argv[1]):
        its = line.split()
        covs.append(int(its[1+idx]))

    diff = [abs(a-b) for a, b in zip(covs[0:-1],covs[1:])]
    rate = [abs(a-b)*2/(a+b) if a+b > 0 else 0 for a, b in zip(covs[0:-1],covs[1:])]
    #tb.show_hist(diff)
    #tb.show_hist(covs)
    tb.show_hist(rate)

def calculate_similarity(cigar):
    ref_length = 0
    matches = 0
    i = 0
    while i < len(cigar):
        # Extract number
        num_str = ''
        while i < len(cigar) and cigar[i].isdigit():
            num_str += cigar[i]
            i += 1
        if i >= len(cigar):
            break
        # Extract operation
        op = cigar[i]
        num = int(num_str) if num_str else 0
        # Update reference length and matches
        if op in ('=', 'X', 'D', 'I'):
            ref_length += num
        if op == '=':
            matches += num
        i += 1
    return matches / ref_length if ref_length > 0 else 0.0


def calc_local_similarity(cigar, winsize):
    pairs = []
    i = 0
    while i < len(cigar):
        # Extract number
        num_str = ''
        while i < len(cigar) and cigar[i].isdigit():
            num_str += cigar[i]
            i += 1
        if i >= len(cigar):
            break
        # Extract operation
        op = cigar[i]
        num = int(num_str) if num_str else 0
        # Update reference length and matches
        if op in ('=', 'M'):
            for j in range(num):
                pairs.append((1, 1))
        elif op in ('X', 'D'):
            for j in range(num):
                pairs.append((0, 1))
        elif op == 'I':
            pairs.append((0, num))
        i += 1

    local = []
    for i in range(len(pairs) - winsize ):
        m = sum([x[0] for x in pairs[i:i+winsize]])
        r = sum([x[1] for x in pairs[i:i+winsize]])
        local.append(m/r if r > 0 else 0.0)

    return local


def ts_cigar_identity(argv):

    similarity = calculate_similarity(argv[0])
    print("Similarity: ",  similarity)

    if len(argv) >= 2:
        print(argv[1])
        winsize = int(argv[1])
        local_similarity = calc_local_similarity(argv[0], winsize)
        for i, v in enumerate(local_similarity):
            print(f"Window {i}: {v:.4f}")
        print("Local Similarity: ", local_similarity)

def ts_bed_range(argv):
    fname = argv[0]
    range = 0
    for line in open(fname):
        its = line.split()
        s, e = int(its[1]), int(its[2])
        range += e - s
    print("Range:", range)
_local_func = locals()
def main():
    utils.script_entry(sys.argv, _local_func, "ts_")

if __name__ == '__main__':
    main()
