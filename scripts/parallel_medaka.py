#!/usr/bin/env python3


import sys, os
import subprocess
import traceback
import argparse
import logging

# medaka_cmd = "singularity exec --containall -B `pwd -P`:`pwd -P` medaka_1.7.2--aa54076.sif medaka"

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt = '%Y-%m-%d %H:%M:%S'
                    )
logger = logging.getLogger("utils")

def required_succ(r):
    if r != 0:
        exit(-1)

def get_ref_name_len(bam):

    name_len = []

    cmd = "samtools view -H %s" % bam
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    for line in  proc.stdout.readlines():
        its = line.decode("utf8").split()
        if its[0] == "@SQ":
            assert len(its) == 3
            n = its[1].split(':')
            l = its[2].split(':')
            assert len(n) == 2 and n[0] == "SN" and len(l) == 2 and l[0] == "LN"

            name_len.append((n[1], int(l[1])))
    
    proc.wait()
    name_len.sort(key = lambda x: -x[1])
    logger.info(name_len)
    return name_len


def group_ref_name(name_len, group_size):
    groups = [list() for _ in range(group_size)]

    for i, nl in enumerate(name_len):
        a = i // group_size
        b = i % group_size

        if a % 2 == 0:
           groups[b].append(nl[0])
        else:
           groups[group_size-1-b].append(nl[0])
    return groups

def run_medaka(groups, args):
    procs = []
    for i, g in enumerate(groups):
        if len(g) > 0:
            cmd = "%s consensus %s --threads 2 %s  %s/hdf-%d  --regions %s 2>&1 | tee %s/g%d.log" % (
                args.medaka,
                args.options,
                args.bam,
                args.wrkdir, i ,
                " ".join(g),
                args.wrkdir, i)
            logger.info(cmd)
            proc = subprocess.Popen(cmd, shell=True)
            procs.append(proc)

    for p in procs:
        r = p.wait()
        required_succ(r)
        logger.info(r)

    cmd = "%s stitch --threads %d %s/hdf-* %s %s" % (
           args.medaka,
           args.threads,
           args.wrkdir,
           args.contigs,
           args.polished)
    logger.info(cmd)
    proc = subprocess.Popen(cmd, shell=True)
    required_succ(proc.wait())


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Run medaka in parallel")
    parser.add_argument("contigs", type=str)
    parser.add_argument("bam", type=str)
    parser.add_argument("polished", type=str)
    parser.add_argument("--threads", type=int, default="24")
    parser.add_argument("--medaka", type=str, default="medaka")
    parser.add_argument("--options", type=str, default="")
    parser.add_argument("--wrkdir", type=str, default=".")

    try:
        args = parser.parse_args(sys.argv[1:])
        groups = group_ref_name(get_ref_name_len(args.bam), args.threads//2)
        args.contigs = os.path.abspath(args.contigs)
        args.bam = os.path.abspath(args.bam)
        args.polished = os.path.abspath(args.polished)
        args.wrkdir = os.path.abspath(args.wrkdir)
        run_medaka(groups, args)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()
        exit(-1)
