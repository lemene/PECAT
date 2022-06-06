#!/usr/bin/env python3
import sys, os
import re
import traceback
import argparse
import glob
import multiprocessing as mp
import functools

BIN_PATH, _ = os.path.split(__file__)

def split_vcf(vcf, wrkdir, prefix):
    ifile = open(vcf)

    os.system("mkdir -p %s" % wrkdir)
    
    def write_heads(ofile, heads, n):
        for line in heads:
            if not line.startswith("##contig=") or line.startswith("##contig=<ID="+n):
                ofile.write(line)

    heads = []
    curr, ofile = None, None
    for line in ifile:
        if line.startswith("#"):
            heads.append(line)
        else:
            its = line.split()
            if its[0] != curr:
                curr = its[0]
                ofile = open(os.path.join(wrkdir, curr + "_ctg"+prefix+".vcf"), "w")
                write_heads(ofile, heads, curr)

            ofile.write(line)

def phase_all(contigs, rd2ctg, wrkdir, threads, usehic):
    vcfs = list(glob.glob(os.path.join(wrkdir, "*_ctg.vcf")))

    func = functools.partial(phase_vcf, contigs, rd2ctg, usehic)
    with mp.Pool(threads) as p:
        p.map(func, vcfs)

def phase_vcf(contigs, rd2ctg, usehic, vcf):
    path, fname = os.path.split(vcf)
    ctg = fname[:-8]
    main, ext = os.path.splitext(vcf)
    phased = main  + ".phased" + ext
    haplotag = main  + ".haplotag.bam"
    phreads = main  + ".phased.reads"
    subr2c = main  + ".bam"
    hic = "" if not usehic else main + "_hic.vcf"

    cmd = "whatshap phase --reference %s %s %s %s -o %s" % (contigs, vcf, hic, rd2ctg, phased)
    os.system(cmd)
    cmd = "bgzip -c %s > %s.gz" % (phased, phased)
    os.system(cmd)
    cmd = "tabix -p vcf %s.gz" % phased
    os.system(cmd)

    cmd = "samtools view -b %s %s > %s" % (rd2ctg, ctg, subr2c)
    os.system(cmd)
    
    cmd = "samtools index %s" %  subr2c
    os.system(cmd)

    cmd = "whatshap haplotag --reference %s %s.gz %s -o %s" % (contigs, phased, subr2c, haplotag)
    os.system(cmd)

    cmd = "samtools view %s | %s/fxtools.py fx_get_phased_reads - %s" % (haplotag, BIN_PATH, phreads)
    os.system(cmd)

def collect_phased_reads(wrkdir, phreads):
    with open(phreads, "w") as ofile:
        for f in glob.glob(os.path.join(wrkdir, "*.phased.reads")):
            for line in open(f):
                ofile.write(line)


def main(argv):
    parser = argparse.ArgumentParser("Purge paf files")
    parser.add_argument("contigs", type=str)
    parser.add_argument("vcf", type=str)
    parser.add_argument("rd2ctg", type=str)
    parser.add_argument("phreads", type=str)
    parser.add_argument("--wrkdir", type=str, default="wrkdir")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--hic", type=str, default="")

    try:
        args = parser.parse_args(argv)
        split_vcf(args.vcf, args.wrkdir, "")
        if (len(args.hic) > 0):
            split_vcf(args.hic, args.wrkdir, "_hic")
        phase_all(args.contigs, args.rd2ctg, args.wrkdir, args.threads, len(args.hic) > 0)
        collect_phased_reads(args.wrkdir, args.phreads)
    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()



if __name__ == "__main__":
    main(sys.argv[1:])