#!/usr/bin/env python3
 
import sys,os
import traceback
import gzip
from collections import defaultdict
import argparse
import re
import math
import multiprocessing as mp

mydir = os.path.split(__file__)
sys.path.append(mydir)

import misc

logger = misc.logger

def fx_eval_heterozygous(argv):
    '''Evalute heterozygous'''
    parser = argparse.ArgumentParser("Evalute heterozygous")
    parser.add_argument("reads", type=str, nargs="+")
    parser.add_argument("--kmer", type=int, default=21)
    parser.add_argument("--out", type=str, default="out")
    parser.add_argument("--threads", type=int, default="24")
    parser.add_argument("--memory", type=int, default=4000000000)
    try:
        args = parser.parse_args(argv)

        readsjf = "%s/reads.jf" % (args.out)
        readshisto = "%s/reads.histo" % (args.out)
        result = "%s/result" % (args.out)

        cmd = "mkdir %s && jellyfish count -C -m %d -t %d -s %d %s -o %s" %     \
            (args.out, args.kmer, args.threads, args.memory, " ".join(args.reads), readsjf)
        misc.run_if_modified(args.reads, [readsjf], cmd)

        cmd = "jellyfish histo -t %d %s > %s" % (args.threads, readsjf, readshisto)
        misc.run_if_modified([readsjf], [readshisto], cmd)

        cmd = "genomescope.R -i %s -o %s -k %d > %s" % (readshisto, args.out, args.kmer, result)
        misc.run_if_modified([readshisto], [result], cmd)

        os.system("cat %s" % result)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def pilon_round(contigs, reads1, reads2, out, threads):

    ctgs = out + "/contigs.fasta"
    pols = out + "/polished"
    sam = out + "/reads.sam"
    bam = out + "/reads.bam"
    sorted = out + "/reads.sorted.bam"

    cmd = "ln -s -r %s %s" % (contigs, ctgs)
    misc.run_if_modified([contigs], [ctgs], cmd)

    cmd = "bwa index %s" % (ctgs)
    misc.run_if_modified([ctgs], [ctgs+".amb"], cmd)
 



    cmd = "bwa mem -t %d %s %s %s > %s" % (threads, ctgs, reads1, reads2, sam)
    misc.run_if_modified([ctgs], [sam], cmd)

    cmd = "samtools view -@%d %s -b -o %s" % (threads, sam, bam)
    misc.run_if_modified([sam], [bam], cmd)

    cmd = "samtools sort -@%d %s -o %s" % (threads, bam, sorted)
    misc.run_if_modified([bam], [sorted], cmd)

    cmd = "samtools index %s" % (sorted)
    misc.run_if_modified([sorted], [sorted+".bai"], cmd)

    cmd = "java -Xmx32G -jar ~/tool/pilon-1.24.jar --genome %s --fix snps,indels,gaps  --bam %s --threads %d --output %s" % (ctgs, sorted, threads, pols)
    misc.run_if_modified([ctgs], [pols+".fasta"], cmd)


def run_pilon(n, ctg, bam, wrkdir):
    subctg = wrkdir + ("/%s.fasta" % n)
    subbam = wrkdir + ("/%s.bam" % n)
    subpol = wrkdir + ("/pol_%s" % n)

    cmd = "samtools faidx %s %s > %s" % (ctg, n, subctg)
    misc.run_if_modified([ctg], [subctg], cmd)

    cmd = "samtools view -b %s %s > %s" % (bam, n, subbam)
    misc.run_if_modified([bam], [subbam], cmd)

    cmd = "samtools index %s" % (subbam)
    misc.run_if_modified([subbam], [subbam+".bai"], cmd)

    cmd = "java -Xmx32G -jar ~/tool/pilon-1.24.jar --genome %s --fix snps,indels,gaps  --bam %s --output %s" % (subctg, subbam, subpol)
    misc.run_if_modified([subctg, subbam], [subpol+".fasta"], cmd)


def pilon_round2(contigs, reads1, reads2, out, threads, clean):

    ctgs = out + "/contigs.fasta"
    pols = out + "/polished.fasta"
    sorted = out + "/reads.sorted.bam"

    subdir = out + "/sub"

    if misc.file_newer(ctgs, pols):
        misc.run_cmd("mkdir -p %s" % subdir)

        cmd = "ln -s -r %s %s" % (contigs, ctgs)
        misc.run_if_modified([contigs], [ctgs], cmd)

        cmd = "bwa index %s" % (ctgs)
        misc.run_if_modified([ctgs], [ctgs+".amb"], cmd)
    
        cmd = "samtools faidx %s" % (ctgs)
        misc.run_if_modified([ctgs], [ctgs+".fai"], cmd)


        cmd = "bwa mem -t %d %s %s %s | samtools sort -@%d --output-fmt BAM -o %s" % (threads, ctgs, reads1, reads2, threads, sorted)
        misc.run_if_modified([ctgs], [sorted], cmd)

        cmd = "samtools index %s" % (sorted)
        misc.run_if_modified([sorted], [sorted+".bai"], cmd)

        ctgnames = []
        for line in open(ctgs+".fai"):
            its = line.split()
            ctgnames.append((its[0], int(its[1])))

        ctgnames.sort(key=lambda x: x[1])

        pool = mp.Pool(threads)
        for n, _ in ctgnames:
            pool.apply_async(run_pilon, (n, ctgs, sorted, subdir))
        pool.close()
        pool.join()

        misc.run_cmd("cat %s/pol_*.fasta >> %s" % (subdir, pols))

        if clean:
            logger.info("remove intermediate files")
            logger.info("remove %s" % subdir)
            logger.info("remove %s*" % sorted)
            logger.info("remove %s.*" % ctgs)
            misc.run_cmd("rm -rf %s" % (subdir))
            misc.run_cmd("rm -rf %s*" % (sorted))
            misc.run_cmd("rm -rf %s.*" % (ctgs))

    return pols

def fx_pilon(argv):
    parser = argparse.ArgumentParser("run pilon")
    parser.add_argument("contigs", type=str)
    parser.add_argument("reads1", type=str)
    parser.add_argument("reads2", type=str)
    parser.add_argument("--iteration", type=int, default=1)
    parser.add_argument("--output", type=str, default="output")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--clean", type=int, default=1)
   
    try:
        args = parser.parse_args(argv)

        ctgs = args.contigs

        for i in range(args.iteration):
            pols = args.output + ("/%d/polished.fasta" % i)

            if misc.file_newer(ctgs, pols):

                wrkdir = args.output + ("/%d" % i)
                os.system("mkdir -p %s" % wrkdir)
                pilon_round2(ctgs, args.reads1, args.reads2, wrkdir, args.threads, args.clean)
                ctgs = pols
   
    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_gcpp(argv):
    parser = argparse.ArgumentParser("run pbmm2 and gcpp")
    parser.add_argument("contigs", type=str)
    parser.add_argument("reads", type=str)
    parser.add_argument("polished", type=str)
    parser.add_argument("--pbmm2_options", type=str, default="")
    parser.add_argument("--gcpp_options", type=str, default="--algorithm=arrow -x 5 -X 120 -q 0")
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--wrkdir", type=str, default="temp")

    try:
        args = parser.parse_args(argv)
        rd2ctg = "rd2ctg.bam"

        cmd = "pbmm2 align -j %d --sort %s %s %s %s" % (args.threads, args.pbmm2_options, args.contigs, args.reads, rd2ctg)
        misc.run_if_modified([args.contigs, args.reads], [rd2ctg], cmd)

        cmd = "gcpp -j %d %s -r %s %s -o %s " % (args.threads, args.gcpp_options, args.contigs, rd2ctg, args.polished)
        misc.run_if_modified([args.contigs, rd2ctg], [args.polished], cmd)

        
    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()


def fx_mummerplot(argv):
    parser = argparse.ArgumentParser("run mummerplot")
    parser.add_argument("ctgs", type=str)
    parser.add_argument("refs", type=str)
    try:
        args = parser.parse_args(argv)
        ctg = args.ctgs
        ref = args.refs

        cmd = "nucmer --mum -l 100 -c 1000 -d 10 --banded -D 5 %s %s" % (ref, ctg)
        os.system(cmd)
        cmd = "delta-filter -i 95 -o 95 out.delta > out.best.delta"
        os.system(cmd)
        cmd = "dnadiff -d out.best.delta"
        os.system(cmd)
        cmd = "mummerplot out.best.delta --fat -f --png"
        os.system(cmd)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()

if __name__ == '__main__':
    if len(sys.argv) > 1:
       locals()[sys.argv[1]](sys.argv[2:])
    else:
       for func in list(locals().keys()):
           if func.startswith("fx_"):
               print("%s: %s" % (func, locals()[func].__doc__.split("\n")[0]))