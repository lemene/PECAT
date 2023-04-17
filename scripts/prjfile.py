import logging
import sys,os
from collections import defaultdict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("utils")

'''fsa工程文件的读取函数'''


class BinInfos(object):
    '''对于fsa_kmer_bin产生的文件
    id len patkmer matkmer offsetkmer errkmer type
'''

    def __init__(self, fname):
        self.load(fname)
        
    def load(self, fname):
        
        logger.info("加载bininfos（fsa_kmer_bin产生）信息")
        self.classified = [set(), set(), set()]
        self.infos = {}
        self.kmers = {}
        for line in open(fname):
            its = line.split()      # name ... type[1|0|-1]
            self.infos[its[0]] = int(its[-1])
            self.kmers[its[0]] = int(its[2]), int(its[3])
            self.classified[int(its[-1])+1].add(its[0])
    

    def test(self, n0, n1):
        return self.infos[n0] * self.infos[n1]



class Binfo(object):
    def __init__(self, fname):
        self.binfo = self.load(fname)
    

    def load(self, fname):
        '''name len pat_only_kmer mat_only_kmer ok_kmer err_kmer type'''
        binfo = {}
        for line in open(fname):
            its = line.split()
            binfo[its[0]] = (its[0], int(its[1]), int(its[2]), int(its[3]), int(its[4]), int(its[5]), int(its[6]) )

        return binfo

    def stat(self):
        kmers = [0, 0, 0, 0] # all, err, hap, haperr
        for its in self.binfo.values():
            kmers[0] += int(its[4]) + int(its[5])
            kmers[1] += int(its[5])
            kmers[2] += int(its[2]) + int(its[3])
            kmers[3] += min(int(its[2]), int(its[3]))
        return kmers

    def cmp(self, r0, r1):
        assert r0 in self.binfo and r1 in self.binfo
        if self.binfo[r0][6] == 0 or self.binfo[r1][6] == 0:
            return 0
        elif self.binfo[r0][6] == self.binfo[r1][6]:
            return 1
        else:
            return -1

def find_item_in_dir(item, dir):
    for i in os.listdir(dir):
        if item == i:
            return os.path.join(dir, item)
    return ""
    
def find_prjpath(item, depth=3):
    cwd = os.getcwd()

    for _ in range(depth):
        path = find_item_in_dir(item, cwd)
        if path != "":
            return path

        cwd = cwd + "/.."

    return ""


def load_overlap_ids_file(fname, ols, idx):
    for line in open(fname):
        its = line.split()
        m0, m1 = sorted([its[idx[0]], its[idx[1]]])
        ols.add((m0, m1))
        
def load_overlap_ids_txt(fname, ols):
    for line in open(fname):
        n = line.strip()
        if len(n) > 0:
            load_overlap_ids_xxx(n, ols)

def load_overlap_ids_xxx(fname, ols):
    logger.info("load overps ids:%s" % fname)
    if fname.endswith(".txt"):
        load_overlap_ids_txt(fname, ols)
    elif fname.endswith(".paf"):
        load_overlap_ids_file(fname, ols, (0,5))
    elif fname.endswith(".m4a"):
        load_overlap_ids_file(fname, ols, (0, 1))
    else:
        logger.error("不支持文件 %s" % fname)


def load_overlap_ids(fname):
    ols = set()
        
    import glob
    for fn in glob.glob(fname):
        load_overlap_ids_xxx(fn, ols)
    return ols


def detect_overlap_name_postion(fname):
    if fname.endswith(".paf") or fname.endswith(".paf.gz"):
        return (0, 5)
    elif fname.endswith(".m4") or fname.endswith(".m4.gz"):
        return (0, 1)
    elif fname.endswith(".m4a") or fname.endswith(".m4a.gz"):
        return (0, 1)
    else:
        assert 0, "Failed to recognize overlap format"

def load_snps_in_reads(fname):
    logger.info("Load snps of reads in file: %s" % fname)
    readsnps = defaultdict(list)
    item_start = 4
    for line in open(fname):
        its = line.strip().split()
        assert len(its) >= item_start
        ctg = its[0]
        snps = []
        for i in its[item_start:]:
            p, _, _, v, *_ = i.split("|")
            if v != '-1':              # skip -1
                snps.append((p, v))
        if len(snps) > 0:
            readsnps[its[1]].append((ctg, snps))

    return readsnps