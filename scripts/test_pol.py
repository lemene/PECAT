#!/usr/bin/env python3

"""暂时存放测试polish工程的脚本"""

import sys
sys.path.insert(0, "/home/niefan/work/mbio")
import mbio
import mbio.ftype.table as tb
import mbio.utils.utils as utils
import argparse
import matplotlib.pyplot as plt

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

def load_from_gaep(fname):
    result = []
    for i, line in enumerate(open(fname)):
        if i % 4 == 3:
            its = line.split()
            if its[3] == 'inv':
                chr, start, end, = its[2], its[6], its[7]
            else:
                chr, start, end, = its[2], its[5], its[6]
            result.append((chr, int(start), int(end), line))
            if result[-1][1] >= result[-1][2]:
                print(result[-1])
                assert 0
    result.sort()
    return result
def load_from_bed1(fname):
    beds = []
    for line in open(fname):
        its = line.split()
        beds.append((its[0], int(its[1]), int(its[2])))
    return beds


def ts_save_coverage_graph(argv):
    import matplotlib.pyplot as plt
    from collections import defaultdict

    bed_fname = argv[0]
    cov_fname = argv[1]

    line = open(bed_fname, 'r').readline()
    if len(line.split()) == 1:
        # If the first line has only one column, it's a GAEP file
        beds = load_from_gaep(bed_fname)
    else:
        # Otherwise, it's a standard BED file  
        beds = load_from_bed1(bed_fname)
        

    covs = defaultdict(list)
    for line in open(cov_fname):
        its = line.split()
        assert len(its) >= 3
        bed = its[0].split(':')
        
        cov = [float(i) for i in its[1:4]]
        covs[bed[0]].append(cov)

    for bed in beds:
        assert bed[0] in covs, f"Bed {bed[0]} not found in coverage data"
        s = bed[1] // 200
        e = (bed[2] + 199) // 200
        s = max(0, s - 50)
        e = min(e + 50, len(covs[bed[0]]))

        cov_data = covs[bed[0]]

        x = [i*200 for i in range(s, e)]
        for i in range(3):
            y = [cov[i] for cov in cov_data[s:e]]
            plt.plot(x, y, label=f"Coverage {i+1}")
        plt.savefig(f"{bed[0]}_{bed[1]}_{bed[2]}.png")
        plt.close()
    

def ts_show_cov_graph(argv):
    import matplotlib.pyplot as plt
    from collections import defaultdict

    bed_str = argv[0]
    cov_fname = argv[1]

    ctg, start_end = bed_str.split(':')
    start, end = map(int, start_end.split('-')) 

    covs = defaultdict(list)
    for line in open(cov_fname):
        its = line.split()
        assert len(its) >= 3
        bed = its[0].split(':')
        
        cov = [float(i) for i in its[1:4]]
        covs[bed[0]].append(cov)

    assert ctg in covs, f"Bed {ctg} not found in coverage data"
    s = start // 200
    e = (end + 199) // 200
    s = max(0, s - 50)
    e = min(e + 50, len(covs[bed[0]]))

    cov_data = covs[bed[0]]

    x = list(range(s, e))
    for i in range(3):
        y = [cov[i] for cov in cov_data[s:e]]
        plt.plot(x, y, label=f"Coverage {i+1}")
    plt.show()
    
def ts_kmer_to_int(argv):
    kmer = argv[0]
    k = len(kmer)
    val = 0
    for i in range(k):
        val <<= 2
        if kmer[i] == 'A':
            val |= 0
        elif kmer[i] == 'C':
            val |= 1
        elif kmer[i] == 'G':
            val |= 2
        elif kmer[i] == 'T':
            val |= 3
        else:
            raise ValueError(f"Invalid character in kmer: {kmer[i]}")
    print(val)
    print(f"{val:0{k*2}b}")
    print(f"{val:0{k*2}x}")
    print(f"{val:0{k*2}X}") 

def ts_show_win_cov(argv):
    """显示窗口覆盖度变化率曲线图"""
    parser = argparse.ArgumentParser(ts_show_win_cov.__doc__)
    parser.add_argument("fname", type=str, help="覆盖度文件名")
    parser.add_argument("--range", type=lambda s: tuple(int(n) for n in s.split(',')), default="2,100")

    try:
        args = parser.parse_args(argv)
        x, ys = [], []
        colnames = []
        for lineno, line in enumerate(open(args.fname)):
            if lineno == 0:
                colnames = line.strip().split()[1:]
                ys = [[] for _ in colnames]
                continue

            its = line.split()
            p = int(its[0])
            if p < args.range[0]:
                continue
            if p > args.range[1]:
                break
            x.append(int(its[0]))
            for i, y in enumerate(its[1:]):
                ys[i].append(float(y))
        
        plt.figure()
        for n, y in zip(colnames, ys):
            plt.plot(x, y, label=n)
        plt.legend()
        plt.title(f"Coverage")
        plt.xlabel("Position")
        plt.show()  

    except SystemExit:
        print("参数解析错误！")
        return


_local_func = locals()
def main():
    utils.script_entry(sys.argv, _local_func, "ts_")

if __name__ == '__main__':
    main()
