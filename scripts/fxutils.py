#!/usr/bin/env python3
import sys, os
import traceback
import gzip
import argparse
from collections import defaultdict
import subprocess
import re

mydir = os.path.split(__file__)
sys.path.append(mydir)
from misc import *
import misc
import prjfile as prj



def fx_fa2fq(argv):
    parser = argparse.ArgumentParser("Convert fasta to fastq")
    parser.add_argument("fasta", type=str)
    parser.add_argument("fastq", type=str)
    parser.add_argument("--quality", type=str, default="g")
    try:
        args = parser.parse_args(argv)
        from Bio import SeqIO
        with open(args.fastq, "w") as fastq:
            for i, rec in enumerate(SeqIO.parse(args.fasta, "fasta")):
                fastq.write("@%s\n%s\n+\n%s\n" % (rec.name, rec.seq, "g"*len(rec.seq)))
    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_get_phased_reads(argv):
    parser = argparse.ArgumentParser("从bam文件中得到phased reads")
    parser.add_argument("sam", type=str)
    parser.add_argument("phased", type=str)

    phased = defaultdict(lambda : [set(),set()])
    try:
        args = parser.parse_args(argv)

        for line in open_file(args.sam):
            if not line.startswith("@"):
                its = line.split()
                ps, hp = "", 0
                for i in its:
                    if i.startswith("HP:i:"):
                        hp = int(i[5:])
                    if i.startswith("PS:i:"):
                        ps = i[5:]

                if hp != 0:
                    phased[its[2]+"-"+ps][hp-1].add(its[0])
                #else:
                #    untagged.append(its[0])

        with open(args.phased, "w") as ofile:
            for k, v in phased.items():
                ofile.write("%s:%s\n" % (",".join(v[0]), ",".join(v[1])))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_purge_paf(argv):
    
    parser = argparse.ArgumentParser("Purge paf files")
    parser.add_argument("paf", type=str)
    parser.add_argument("tile", type=str)
    parser.add_argument("id2name", type=str)
    parser.add_argument("phased", type=str)
    parser.add_argument("--count", type=int, default=2)


    try:
        args = parser.parse_args(argv)
        from collections import defaultdict
        groups = defaultdict(set)

        logger.info("Loading phased reads")
        for line in open(args.phased):
            its = line.strip().split(":")
            assert len(its) == 2
            its[0] = its[0].split(',')
            its[1] = its[1].split(',')
            for i in its[0]:
                groups[i].update(its[1])
            for i in its[1]:
                groups[i].update(its[0])

        logger.info("Loading reads in contigs")
        ctgs = PrjFile.get_contig_reads(args.tile, args.id2name)

        logger.info("Collecting reads in differenct haplotype")
        ctgs_v = {}
        for k, v in ctgs.items():
            its = defaultdict(int)
            for iv in v:
                for g in groups[iv]:
                    its[g] += 1
            ctgs_v[k] = its

        logger.info("Filtering overlaps, count=%d" % (args.count))
        for line in open(args.paf):
            its = line.split()
            c = ctgs_v[its[5]].get(its[0], 0)
            if c < args.count:
                print(line, end="")


    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_purge_paf2(argv):
    
    parser = argparse.ArgumentParser("Purge paf files")
    parser.add_argument("paf", type=str)
    parser.add_argument("tile", type=str)
    parser.add_argument("readsnps", type=str)


    try:
        args = parser.parse_args(argv)

        filter = PurgeInconsistentInPaf(args.tile, args.readsnps)
        
        logger.info("Filtering overlaps" )
        for line in open(args.paf):
            if not filter(line):
                print(line, end="")

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_purge_overlaps(argv):
    
    parser = argparse.ArgumentParser("Purge paf files")
    parser.add_argument("overlaps", type=str)
    parser.add_argument("tile", type=str)
    parser.add_argument("readsnps", type=str)


    try:
        args = parser.parse_args(argv)

        if args.overlaps.endswith(".paf"):
            filter = PurgeInconsistentInPaf(args.tile, args.readsnps)
        elif args.overlaps.endswith(".sam"):
            filter = PurgeInconsistentInSam(args.tile, args.readsnps)
        else:
            filter = None
        logger.info("Filtering overlaps" )
        for line in open(args.overlaps):
            if not filter(line):
                print(line, end="")

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_purge_rd2ctg(argv):
    
    parser = argparse.ArgumentParser("Purge rd2ctg files")
    parser.add_argument("ifname", type=str)
    parser.add_argument("ofname", type=str)
    parser.add_argument("tile", type=str)
    parser.add_argument("readsnps", type=str)
    parser.add_argument("id2name", type=str)
    parser.add_argument("--threads", type=int, default=10)


    try:
        args = parser.parse_args(argv)
        block_size = 5000

        filter = PurgeInconsistentInSam(args.tile, args.readsnps, args.id2name)
        
        threads = (args.threads / 2, 1 ,args.threads / 2)

        cmd0 = "samtools view -h -@ %d %s" % (threads[0], args.ifname)
        cmd1 = "samtools view -b -@ %d -o %s" % (threads[2], args.ofname)
        ip = subprocess.Popen(cmd0, shell=True, stdout=subprocess.PIPE, cwd=os.getcwd(),encoding="utf-8")
        op = subprocess.Popen(cmd1, shell=True, stdin=subprocess.PIPE, cwd=os.getcwd(),encoding="utf-8")
        count = [0, 0]
        for line in ip.stdout:
            count[0] += 1
            if line.startswith("@"):
                print(line, file=op.stdin, end="")
                count[1] += 1
            else:
                if not filter(line):
                    print(line, file=op.stdin, end="")
                    count[1] += 1

            if count[0] % block_size == 0:
                logger.info("done: %d/%d" % (count[1], count[0]))

        logger.info("reserves: %d/%d" % (count[1], count[0]))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()



def fx_stat_runtime(argv):
    parser = argparse.ArgumentParser("stat plgd running time")
    parser.add_argument("logdir", type=str)
    parser.add_argument("--threads", type=int, default=1)

    try:
        args = parser.parse_args(argv)

        import re, glob, datetime
        
        stats = defaultdict(datetime.timedelta)
        start_prefix = "Plgd script start:"
        end_prefix = "Plgd script end:"
        for f in glob.glob(os.path.join(args.logdir, "*.sh.log")):
            _, n = os.path.split(f)
            t = n.split('_')[0]

            s, e = None, None
            for line in open(f):
                if re.match(start_prefix, line):
                    s = datetime.datetime.strptime(line[len(start_prefix):].strip(), "%Y-%m-%d %H:%M:%S")

                if re.match(end_prefix, line):
                    e = datetime.datetime.strptime(line[len(end_prefix):].strip(), "%Y-%m-%d %H:%M:%S")

            assert s != None and e != None and e >= s
            stats[t] += e - s

        total = datetime.timedelta()
        for k, v in stats.items():
            total += v
            print("%s:\t%s, %s" % (k, v, v*args.threads))

        print("Total:\t%s, %s" % (total, total*args.threads))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()



def fx_kmer(argv):
    parser = argparse.ArgumentParser("s")
    parser.add_argument("kmers", type=str)
    parser.add_argument("reads", type=str)

    try:
        args = parser.parse_args(argv)

        logger.info("load kmer")
        import re, glob, datetime
        from Bio import SeqIO
        kmers = set()
        for line in open(args.kmers):
            kmers.add(line.split()[0])
        logger.info("load kmers: %d" % len(kmers))

        for i in kmers: k = len(i); break
        logger.info("kmer len: %d" % k)

        def get_kmers(s, k):
            return set([s[i:i+k] for i in range(0, len(s)-k+1)])
        
        ifile, iftype = open_seq_file(args.reads)
        for i, rec in enumerate(SeqIO.parse(ifile, iftype)):
            ks = get_kmers(rec.seq, k)
            count = len(ks & kmers)
            if count > 0:
                print(rec.name, count)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_filter_paf_with_bininfos(argv):
    '''使用bin信息过滤paf文件'''
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("ifname", type=str)
    parser.add_argument("ofname", type=str)
    parser.add_argument("bininfos", type=str)

    try:
        args = parser.parse_args(argv)

        logger.info("加载bininfos（fsa_kmer_bin产生）信息")

        classified = [set(), set(), set()]
        for line in open(args.bininfos):
            its = line.split()      # name ... type[1|0|-1]
            classified[int(its[-1])+1].add(its[0])

        logger.info("过滤")
        with open(args.ofname, "w") as ofile:
            for i, line in enumerate(open(args.ifname)):
                its = line.split()
                s0, s1 = its[0], its[5]
                if s0 in classified[0] and s1 in classified[2] or s0 in classified[2] and s1 in classified[0]:
                    pass
                else:
                    ofile.write(line)

                if i % 1000000 == 0:
                    print(i)
    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def detect_overlap_name_postion(fname):
    if fname.endswith(".paf") or fname.endswith(".paf.gz"):
        return (0, 5)
    elif fname.endswith(".m4") or fname.endswith(".m4.gz"):
        return (0, 1)
    elif fname.endswith(".m4a") or fname.endswith(".m4a.gz"):
        return (0, 1)
    else:
        assert 0, "Failed to recognize overlap format"

def fx_eval_overlap_haplotype(argv):
    parser = argparse.ArgumentParser("使用reads的分类信息（fsa_kmer_bin生成），评估overlaps中有多少单倍型混合")
    parser.add_argument("overlaps", type=str)
    parser.add_argument("binfos", type=str)


    try:
        args = parser.parse_args(argv)

        logger.info("加载binfos（fsa_kmer_bin产生）信息")
        binfos = prj.BinInfos(args.binfos)

        logger.info("检查 overlaps")
        pos = detect_overlap_name_postion(args.overlaps)
        count = [0, 0, 0]
        for line in open_file(args.overlaps):
            its = line.split()
            r0, r1 = its[pos[0]], its[pos[1]]
            count[1+binfos.test(r0, r1)] += 1

        print(count, 1 - count[0]/sum(count))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_eval_phased(argv):
    parser = argparse.ArgumentParser("使用reads的分类信息（fsa_kmer_bin生成），检查phased reads的正确率")
    parser.add_argument("binfos", type=str)
    parser.add_argument("phased", type=str)

    try:
        args = parser.parse_args(argv)

        logger.info("加载binfos（fsa_kmer_bin产生）信息")

        
        binfo = prj.Binfo(args.binfos)

        done = dict()
        logger.info("加载phased信息")
        for line in open(args.phased):
            its = line.strip().split(':')
            set0 = its[0].split(',')
            set1 = [i.split("|")[0] for i in its[1].split(',')[:-1]]
            # assert len(set0) == 1
            for s0 in set0:
                for s1 in set1:
                    if (s0,s1) in done or (s1, s0) in done: continue
                    done[(s0, s1)] = binfo.cmp(s0, s1)

        count = sum([1 for _, v in done.items() if v == -1])

        print("%d / %d = %.02f" % (count, len(done), count/len(done) * 100))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_eval_phased3(argv):
    parser = argparse.ArgumentParser("使用reads的分类信息（fsa_kmer_bin生成），检查phased reads的正确率")
    parser.add_argument("binfos", type=str)
    parser.add_argument("overlaps", type=str)
    parser.add_argument("filtered", type=str)
    parser.add_argument("--detail", action="store_true")

    try:
        args = parser.parse_args(argv)

        logger.info("加载binfos（fsa_kmer_bin产生）信息")
        binfo = prj.Binfo(args.binfos)
        logger.info("加载overlaps信息")

        ols = prj.load_overlap_ids(args.overlaps)
        flt = prj.load_overlap_ids(args.filtered)
        # ols, flt are sorted

        TP, TN, FP, FN = 0, 0, 0, 0
        for r0, r1 in ols:
            s = binfo.cmp(r0, r1)
            if s == -1:
                if (r0, r1) not in flt:
                    TP += 1
                else:
                    FN += 1
            else:
                if (r0, r1) not in flt:
                    if args.detail: print(r0, r1)
                    FP += 1
                else:
                    TN += 1

        print("TP = %d, TN = %d, FP = %d, FN = %d" % (TP, TN, FP, FN), TP, TN, FP, FN)
        precision = TP/(TP+FP)
        recall = TP/(TP+FN)
        f1_score = 2* (precision*recall) / (precision + recall) 
        print("precision   = %.02f" % precision)
        print("recall      = %.02f" % recall)
        print("F1-score    = %.02f" % f1_score)


    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_compare_binfos(argv):
    parser = argparse.ArgumentParser("比较两个BIN的信息")
    parser.add_argument("binfos0", type=str)
    parser.add_argument("binfos1", type=str)
    parser.add_argument("--method", type=str, default="")

    try:
        args = parser.parse_args(argv)
        binfos0 = prj.Binfo(args.binfos0)
        binfos1 = prj.Binfo(args.binfos1)

        if args.method == "":
            count = [0, 0]
            for n, v in binfos1.infos.items():
                #if (v * binfos0.infos[n] == -1):
                if (v != binfos0.infos[n]):

                    count[1] += 1
                count[0] += 1

            print(count, 1-count[1]/count[0])
        elif args.method == "read":
            for n, v0 in binfos0.kmers.items():
                if n in binfos1.kmers:
                    v1 = binfos1.kmers[n]
                    x0 = v0[0] / sum(v0) if sum(v0) > 0 else 0
                    y0 = v0[1] / sum(v0) if sum(v0) > 0 else 0
                    x1 = v1[0] / sum(v1) if sum(v1) > 0 else 0
                    y1 = v1[1] / sum(v1) if sum(v1) > 0 else 0
                    print(n, x0, y0, x1, y1)
        elif args.method == "read1":
            for n, v0 in binfos0.kmers.items():
                if n in binfos1.kmers:
                    v1 = binfos1.kmers[n]
                    r0 = min(v0) / sum(v0) if sum(v0) > 0 else 0
                    r1 = min(v1) / sum(v1) if sum(v1) > 0 else 0
                    print(n, r0, r1)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_eval_inconsistent_overlaps(argv):
    parser = argparse.ArgumentParser("比较两个BIN的信息")
    parser.add_argument("binfos", type=str)
    parser.add_argument("--detail", action="store_true")
    parser.add_argument("overlaps", type=str)

    try:
        args = parser.parse_args(argv)
        binfos = prj.BinInfos(args.binfos)

        pos = [0, 5] if args.overlaps.endswith(".paf") else [0, 1]

        count = [0, 0, 0]
        for line in open(args.overlaps):
            its = line.split()
            q, t = its[pos[0]], its[pos[1]]
            c = binfos.infos[q] * binfos.infos[t]
            count [c+1] += 1
            if c == -1: print(q, t)

        print(count, count[0]/(count[0] + count[2]), count[2]/(count[0] + count[2]))
        

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_eval_inconsistent(argv):
    parser = argparse.ArgumentParser("evaluate inconsistent overlaps/edges/phased")
    parser.add_argument("binfos", type=str)
    parser.add_argument("--detail", action="store_true")
    parser.add_argument("fname", type=str)

    try:
        args = parser.parse_args(argv)
        binfos = prj.BinInfos(args.binfos)

        if args.fname.endswith(".paf"):
            def parse_line(line):
                its = line.split()
                return ((its[0], its[5]),)
        elif args.fname.endswith(".m4a"):
            def parse_line(line):
                its = line.split()
                return ((its[0], its[1]),)
        elif args.fname.endswith("graph_edges.gz"):
            def parse_line(line):
                its = line.split()
                if its[-1] == "active":
                    return ((its[0][:-2],its[1][:-2]),)
                else:
                    return ()
        elif args.fname.endswith("inconsistent"):
            def parse_line(line):
                pairs = []
                its = line.strip().split(':')
                set0 = its[0].split(',')
                set1 = [i.split("|")[0] for i in its[1].split(',')[:-1]]
                # assert len(set0) == 1
                for s0 in set0:
                    for s1 in set1:
                        pairs.append((s0, s1))
                return pairs

        else:
            assert 0

        count = [0, 0, 0]
        done = set()
        import glob
        for fn in glob.glob(args.fname):
        
            for line in misc.open_file(fn):
                its = parse_line(line)
                for t, q in its:
                    if (t,q) in done or (q,t) in done: continue

                    done.add((t,q))
                    c = binfos.infos[q] * binfos.infos[t]
                    count [c+1] += 1
                    if c == 1 and args.detail: print(q, t)

        print(count, count[0]/sum(count), count[2]/sum(count))
        

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_analyze_snps(argv):
    parser = argparse.ArgumentParser("分析SNP")
    parser.add_argument("variants", type=str)

    try:
        args = parser.parse_args(argv)

        for line in open(args.variants):
            its = line.split()
            name, pos, *vs = its
            vs = [int(v) for v in vs]
            s = sum(vs)
            count = sum([1 for v in vs if v >= 0.2*s])
            if count >= 2:
                print(pos)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_compare_variants(argv):
    parser = argparse.ArgumentParser("比较两个Variants文件")
    parser.add_argument("vars0", type=str)
    parser.add_argument("vars1", type=str)
    parser.add_argument("--valid0", type=str)
    parser.add_argument("--valid1", type=str)
    parser.add_argument("--common", type=str)

    try:
        args = parser.parse_args(argv)

        def load_snps(fname):
            snps = {}
            for line in open(fname):
                its = line.split()
                if its[-1] == '1':
                    snps[(its[0], int(its[1]))] = [int(i) for i in its[2:-1]]
            return snps

        snps0 = load_snps(args.vars0)
        snps1 = load_snps(args.vars1)

        comm = set(snps0.keys()).intersection(set(snps1.keys()))

        logger.info("snps0 size: %d" % len(snps0))
        logger.info("snps1 size: %d" % len(snps1))
        logger.info("common size: %d" % len(comm))

        def calc_freq(v):
            sv = list(sorted(v))
            return sv[-3] / sv[-1] # secondary / cov

        if args.valid0:
            logger.info("output: %s" % args.valid0)
            ofile = open(args.valid0, "w")
            for k in sorted(snps0.keys()):
                print("_".join([str(ik) for ik in k]), calc_freq(snps0[k]), file=ofile)

        if args.valid1:
            logger.info("output: %s" % args.valid1)
            ofile = open(args.valid1, "w")
            for k in sorted(snps1.keys()):
                print("_".join([str(ik) for ik in k]), calc_freq(snps1[k]), file=ofile)

        if args.common:
            logger.info("output: %s" % args.common)
            ofile = open(args.common, "w")
            for k in comm:
                print("_".join([str(ik) for ik in k]), calc_freq(snps0[k]), calc_freq(snps1[k]), file=ofile)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_compare_snps_in_reads(argv):
    parser = argparse.ArgumentParser("比较snps在不同reads的分布")
    parser.add_argument("vars0", type=str)
    parser.add_argument("vars1", type=str)
    parser.add_argument("--snps_in_contigs", type=str)

    try:
        args = parser.parse_args(argv)

        def load_snps_in_reads(fname, comm):
            snps = {}
            for line in open(fname):
                its = line.split()
                count = 0
                for i in its[4:]:
                    print(i)
                    p, *_ = i.split("-")
                    if comm == None or (its[1], int(p)) in comm:
                        count += 1
                snps[(its[1], its[3])] = count
            return snps


        def load_comm_snps(fname):
            comm = set()
            for line in open(fname):
                c, p = line.split()[0].split("_")

                comm.add((c, int(p)))
            return comm

        comm_snps = load_comm_snps(args.snps_in_contigs) if args.snps_in_contigs else None

        snps0 = load_snps_in_reads(args.vars0, comm_snps)
        snps1 = load_snps_in_reads(args.vars1, comm_snps)



        comm = set(snps0.keys()).intersection(set(snps1.keys()))

        logger.info("snps0 size: %d" % len(snps0))
        logger.info("snps1 size: %d" % len(snps1))
        logger.info("common size: %d" % len(comm))

        for k in comm:
            print("_".join([str(ik) for ik in k]), snps0[k], snps1[k])


    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()



def fx_max_haperr(argv):
    parser = argparse.ArgumentParser("比较两个Variants文件")
    parser.add_argument("binfos", type=str)

    try:
        args = parser.parse_args(argv)
        infos = []
        for line in open(args.binfos):
            its = line.split()
            pkmer, mkmer = int(its[2]), int(its[3])
            infos.append((line, min(pkmer, mkmer)))

        infos.sort(key = lambda x:  -x[1])

        for i, it in enumerate(infos):
            if i > 10: break
            print(it[0], end="")

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()



def fx_snp_diff(argv):
    parser = argparse.ArgumentParser("两条reads的SNPs的异同。")
    parser.add_argument("rinfos", type=str)
    parser.add_argument("read0", type=str)
    parser.add_argument("read1", type=str)
    parser.add_argument("--w", type=int, default=0)

    try:
        args = parser.parse_args(argv)
        rinfos = PrjFile.load_snp_in_read2(args.rinfos)

        n0, n1 = args.read0, args.read1
        w = args.w
        print(n0, n0 in rinfos)
        print(n1, n1 in rinfos)
        if n0 in rinfos and n1 in rinfos:
            r0 = rinfos[n0]
            r1 = rinfos[n1]

            for ir0 in r0:
                for ir1 in r1:
                    if ir0[0] != ir1[0]: continue
                    count = [0, 0]
                    detail = defaultdict(set)
                    for k, v in ir0[1].items():
                        if k in ir1[1]:
                            if v[w] == "-1" or ir1[1][k][w] == "-1": continue
                            if v[w] == ir1[1][k][w]:
                                count[0] += 1
                            else:
                                count[1] += 1
                            detail[k].add(v[w])
                            detail[k].add(ir1[1][k][w])

                    #print(detail)
                    print(count)
                    print(ir0[0],ir1[0])
                    print(ir0[0])
    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_falcon_format(argv):
    parser = argparse.ArgumentParser("将fsa生成的ctg名称改为falcon名称格式")
    parser.add_argument("ifname", type=str)
    parser.add_argument("ofname", type=str)
    parser.add_argument("--namemap", action="store_true")
    parser.add_argument("--nfname", type=str, default="")
   
    try:
        args = parser.parse_args(argv)
        

        from Bio import SeqIO
        
        ifile, iftype = open_seq_file(args.ifname)
 
        re_main = re.compile("ctg([0-9]+)")
        re_sub = re.compile("ctg([0-9]+)-([0-9]+)-([0-9]+)$")
        def split_name(rec):
            r = re_main.match(rec.id)
            assert r != None
            main, homo = r.group(1), ""

            r = re_sub.match(rec.id)
            if r:
                homo = main
            else:
                homo = rec.description.split()[1].split(':')[0][3:] 
            return main, homo            

        logger.info("load namemap")
        nfname = args.nfname if len(args.nfname) > 0 else args.ifname
        indexs = defaultdict(int)
        namemap = {}
        nfile, nftype = open_seq_file(nfname)
        for rec in SeqIO.parse(nfile, nftype):
            main, homo = split_name(rec)
            if homo != "":
                name = '0'*(6-len(homo)) + homo + 'F_' + str(indexs[homo]).rjust(3, "0")
                indexs[homo] += 1
            else:
                name = '0'*(6-len(main)) + main + 'F' 
            namemap[rec.name] = name

        if args.namemap:
            print("%s\t%s" % (name, rec.name))

        logger.info("rename")
        ofile, _ = open_seq_file(args.ofname, "w")
        ifile, iftype = open_seq_file(args.ifname)
        for rec in SeqIO.parse(ifile, iftype):
            ofile.write(">%s\n%s\n" % (namemap[rec.name], rec.seq))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_rename(argv):
    parser = argparse.ArgumentParser("修改序列名称")
    parser.add_argument("mapfile", type=str)
    parser.add_argument("ifname", type=str)
    parser.add_argument("ofname", type=str)
    parser.add_argument("--map", type=str, default="0->1")
   
    try:
        args = parser.parse_args(argv)
        
        indexs = defaultdict(int)

        from Bio import SeqIO
        
        ifile, iftype = open_seq_file(args.ifname)
        ofile, oftype = open_seq_file(args.ofname, "w")

        m = args.map.split("->")
        assert len(m) == 2
        m = (int(m[0]), int(m[1]))

        namemaps = {}
        for line in open(args.mapfile):
            its = line.split()
            namemaps[its[m[0]]] = its[m[1]] 

        for i, rec in enumerate(SeqIO.parse(ifile, iftype)):
            assert rec.id in namemaps       
            rec.id = namemaps[rec.id]
            rec.description = rec.id
            SeqIO.write(rec, ofile, oftype)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_stat_binfo(argv):
    parser = argparse.ArgumentParser("统计binfo信息")
    parser.add_argument("binfo", type=str)
   
    try:
        args = parser.parse_args(argv)
        
        binfo = prj.Binfo(args.binfo)
        stat = binfo.stat()
        
        print("All Kmer: %d(%.04f)" % (stat[0], 1))
        print("Err Kmer: %d(%.04f)" % (stat[1], stat[1]/stat[0]))
        print("Hap Kmer:  %d(%.04f)" % (stat[2], stat[2]/stat[0]))
        print("Hap Err:  %d(%.04f, %.04f)" % (stat[3], stat[3]/stat[2], 1 - stat[3]/stat[2]))
        print("Kmer Err:  %d(%.04f)" % ((stat[2]+stat[1]), (stat[2]+stat[1])/stat[0]))


    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()



def fx_eval_select(argv):
    parser = argparse.ArgumentParser("Evalute accuracy of selecting reads for correction")
    parser.add_argument("binfos", type=str)
    parser.add_argument("select", type=str)
    parser.add_argument("--weight", type=str, default="0.2:0.8")
    try:
        args = parser.parse_args(argv)

        logger.info("Loading bin infos from %s" % args.binfos)
        bininfos = {}
        for line in open(args.binfos):
            its = line.split()
            bininfos[int(its[0])] = int(its[-1])  # 1 paternal -1 maternal 0 ambiguous

        judge = lambda a, b: 0 if bininfos[a] == 0 or bininfos[b] == 0 else 1 if bininfos[a] == bininfos[b] else -1

        logger.info("Evaluting infos: %s" % args.select)
        scores = defaultdict(list)
        pair = [0, 0, 0]    # error, amb, ok
        TP, FP, FN, TN = 0, 0, 0, 0
        SEL = 0
        
        import glob
        
        for fsel in glob.glob(args.select):
            logger.info("Load infos: %s" % fsel)
            for line in open(fsel):
                its = line.split()
                a, b, score, select = int(its[0]), int(its[1]), float(its[2]), int(its[3])

                j = judge(a, b)
                pair[1+j] += 1
                if select == 1: 
                    scores[a].append([score, j])
                    SEL += 1

                if select == 1 and j == 1:
                    TP += 1
                elif select == 1 and j == -1:
                    FP += 1
                elif select == 0 and j == 1:
                    FN += 1
                elif select == 0 and j == -1:
                    TN += 1

        SUM = sum(pair)
        
        wts = [0, 0, 0]
        for k, v in scores.items():
            ss = [i[0] for i in v]
            maxs = max(ss)
            mins = min(ss)
            rs = [float(w) for w in args.weight.split(":")]
            linear = lambda x: (x - mins) * (rs[1]-rs[0]) / (maxs-mins) + rs[0] if maxs != mins else sum(rs)/len(rs)
            
            for i in v:
                wts[i[1]+1] += linear(i[0])


        print("Selected in all        : %.02f%%" % ((SEL / SUM) * 100))
        print("Error without selection: %.02f%%" % (pair[0] / SUM *100))
        print("Error with selection   : %.02f%%" % (FP/SEL *100))
        print("Error with weight      : %.02f%%" % (wts[0] / sum(wts)*100))

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_n50(argv):
    parser = argparse.ArgumentParser("Calculate N50")
    parser.add_argument("fname", type=str)
    parser.add_argument("--column", type=int, default=2)
    parser.add_argument("--genome_size", type=int, default=0)

    try:
        args = parser.parse_args(argv)
        lens = [int(line.split()[args.column]) for line in open(args.fname)]

        genome_size = args.genome_size if args.genome_size > 0 else sum(lens)

        lens.sort(key=lambda x: -x)

        acuu = 0
        for l in lens:
            acuu += l
            if acuu >= genome_size / 2:
                print("N50", l)
                break

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_split_polctgs(argv):
    parser = argparse.ArgumentParser("split polished contigs")
    parser.add_argument("polctgs", type=str)
    parser.add_argument("--tool", type=str, default="")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--prjdir", type=str, default="")
    parser.add_argument("--prictgs", type=str)
    parser.add_argument("--altctgs", type=str)

    try:
        from Bio import SeqIO
        
        args = parser.parse_args(argv)
        def load_names(fname):
            names = set()
            for line in open(fname):
                if line.startswith(">"):
                    names.add(line.split()[0][1:])
            return names

        def classify_gcpp(name, pri, alt):
            its = name.split("|")
            assert len(its) == 2
            if its[0] in pri:
                return "primary"
            elif its[0] in alt:
                return "alternate"
            else:
                return "none"

                
        def classify_pilon(name, pri, alt):
            while name.endswith("_pilon"):
                name = name[:-6]
            if name in pri:
                return "primary"
            elif name in alt:
                return "alternate"
            else:
                return "none"

        def detect_tool(fname):
            fpol, tpol = open_seq_file(fname)
            for rec in SeqIO.parse(fpol, tpol):
                if rec.name.endswith("_pilon"):
                    return "pilon"
                elif rec.name.endswith("|arrow"):
                    return "gcpp"
                else:
                    return ""
            else:
                return ""

        tool = args.tool if args.tool != "" else detect_tool(args.polctgs)
        logger.info("tool: %s" % tool)
        classify = classify_gcpp if tool == "gcpp" else classify_pilon if tool == "pilon"  else None

        if args.prictgs != "" and args.altctgs:
            prictgs = args.prictgs
            altctgs = args.altctgs
        elif args.prjdir != "":
            prictgs = args.prjdir + "/6-polish/primary.fasta"
            altctgs = args.prjdir + "/6-polish/alternate.fasta"
        else:
            p = prj.find_prjpath("6-polish", 5)
            assert p != "" and "find empty path"
            prictgs = p + "/primary.fasta"
            altctgs = p + "/alternate.fasta"

        logger.info("prictgs:%s" % prictgs)
        prinames = load_names(prictgs)
        logger.info("altctgs:%s" % altctgs)
        altnames = load_names(altctgs)

        pripol = "primary.fasta"
        altpol = "alternate.fasta"

        assert not os.path.exists(pripol) and not os.path.exists(altpol) or args.force

        
        fpol, tpol = open_seq_file(args.polctgs)
        fpri, tpri = open_seq_file(pripol, "w")
        falt, talt = open_seq_file(altpol, "w")

        for rec in SeqIO.parse(fpol, tpol):
            t = classify(rec.id, prinames, altnames)
            if t == "primary":
                SeqIO.write(rec, fpri, tpri)
            elif t == "alternate":
                SeqIO.write(rec, falt, talt)
            else:
                assert "Unknown"

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

def fx_split_mappings(argv):
    parser = argparse.ArgumentParser("split polished contigs")
    parser.add_argument("maps", type=str)
    parser.add_argument("ctgs", type=str, default="")
    parser.add_argument("subpaf", type=str, default="")
    parser.add_argument("--type", type=str, default="")

    try:
        
        args = parser.parse_args(argv)
        def load_names(fname):
            names = set()
            for line in open(fname):
                if line.startswith(">"):
                    names.add(line.split()[0][1:])
            return names

        ctgs = args.ctgs.split(",")
        subs = args.subpaf.split(",")
        assert len(ctgs) == len(subs)

        logger.info("loading read names in fasta")
        ofiles = [(load_names(c), open(s, "w"))  for c,s  in zip(ctgs,subs)]

        def parse_paf(line):
            its = line.split()
            return [0, its[5]]

        def parse_sam(line):
            if line.startswith("@"):
                return [1, None]
            else:
                its = line.split()
                return [0, its[2]]

        if args.type == "paf" or args.maps.endswith(".paf"):
            parse_mapping = parse_paf
        elif args.type == "sam" or args.maps.endswith(".sam"):
            parse_mapping = parse_sam
        else:
            parse_mapping = None         

        logger.info("starting split mappings")
        for line in sys.stdin if args.maps == "-" else open(args.maps):
            c, n = parse_mapping(line)
            for f in ofiles:
                if c == 1 or n in f[0]:
                    f[1].write(line)
            

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()
        exit(-1)

if __name__ == '__main__':
    script_entry(sys.argv, locals())
