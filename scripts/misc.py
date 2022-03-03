import sys, os


from collections import defaultdict
import logging
import gzip
import multiprocessing as mp;

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt = '%Y-%m-%d %H:%M:%S'
                    )
logger = logging.getLogger("utils")


def run_if_modified(ifiles, ofiles, cmd, force=False, reqsucc=True):
    '''如果ofiles修改时间落后ifiles，则运行func'''

    print(ifiles)
    print(ofiles)
    logger.info("Run: %s" % cmd)
    # 使支持单文件
    _ifiles = [ifiles] if type(ifiles) == str else ifiles
    _ofiles = [ofiles] if type(ofiles) == str else ofiles

    itimes =[ os.path.getmtime(f) for f in _ifiles]
    otimes = [os.path.getmtime(f) if os.path.lexists(f) else 0 for f in _ofiles]

    if len(otimes) == 0 or max(itimes) > min(otimes) or force:
        r = os.system(cmd)
        if reqsucc: assert r == 0

def run_cmd(cmd):
    logger.info("Run: %s" % cmd)
    r = os.system(cmd)
    assert r == 0

def file_newer(f0, f1):
    return not os.path.lexists(f1) or os.path.getmtime(f0) > os.path.getmtime(f1)

def open_file(fname):
    if fname == "-":
        return sys.stdin
    elif fname.endswith(".gz"):
        return gzip.open(fname, "rt") 
    else:
        return open(fname, "r")


def open_seq_file(fname, mode="r"):
    if fname[-6:] == ".fasta":
        return open(fname, mode), "fasta"
    if fname[-4:] == ".fna":
        return open(fname, mode), "fasta"
    elif fname[-6:] == ".fastq":
        return open(fname, mode), "fastq"
    elif fname[-9:] == ".fasta.gz":
        return gzip.open(fname, mode+"t"), "fasta"
    elif fname[-9:] == ".fastq.gz":
        return gzip.open(fname, mode+"t"), "fastq"
    else:
        return None, ""


def sample_seqs(ifname, ofname, base_size):
    from Bio import SeqIO
    import random    

    ifile, iftype = open_seq_file(ifname)

    total_base_size = 0
    for i, rec in enumerate(SeqIO.parse(ifile, iftype)):
        total_base_size += len(rec.seq)

    ratio = (base_size*1.2) / total_base_size   # 1.2 多采样0.5
    
    ifile.seek(0)
    allrec = [] 
    for i, rec in enumerate(SeqIO.parse(ifile, iftype)):
        if random.random() <= ratio:
            allrec.append(rec)

    random.shuffle(allrec)

    accu = 0
    ofile, oftype = open_seq_file(ofname, "w")
    for r in allrec:
        SeqIO.write(r, ofile, oftype)
        accu += len(r.seq)
        if accu >= base_size: break


def script_entry(argv, locals):
    
    if len(argv) > 1:
        locals[argv[1]](argv[2:])
    else:
        maxlen = max([len(func) for func in list(locals.keys())])
        for func in list(locals.keys()):
            if func.startswith("fx_"):
                doc = locals[func].__doc__
                print("%s%s: %s" % (func, " "*(maxlen-len(func)), doc.split("\n")[0] if doc != None else ""))

class PrjFile(object):
    @staticmethod
    def load_id2name(fname, reverse=False):
        id_name = {}
        idx = (0, 1) if not reverse else (1, 0)
        for line in open_file(fname):
            its = line.split()
            id_name[its[idx[0]]] = its[idx[1]]

        return id_name

    @staticmethod
    def get_contig_reads(tile):

        ctgs = defaultdict(list)

        for line in open_file(tile):
            its = line.split()
            s = ctgs[its[0]]
            r = its[1].split("=")[1].split("~")
            a = str(int(r[0].split(":")[0]))
            if len(s) == 0: s.append(a)
            a = str(int(r[1].split(":")[0]))
            s.append(a)

        return ctgs

    @staticmethod
    def load_phased(fname):
        phased = defaultdict(set)
        for line in open(fname):
            its = line.strip().split(':')
            for s in its[1].split(','):
                if len(s) > 0:
                    phased[its[0]].add(s.split('|')[0])

        return phased

    @staticmethod
    def load_snp_in_read(fname):

        rinfos = defaultdict(list)
        for line in open(fname):
            its = line.split()

            ss = {}
            for i in its[4:]:
                p, b, *_ = i.split("-")

                ss[(its[0], p)] = b

            rinfos[its[1]].append(ss)
        return rinfos

    @staticmethod
    def process_line_snp_in_read(line):
            its = line.split() 
            ss = {}
            for i in its[4:]:
                p, b, p1, b1 = i.split("|")

                ss[(its[0], p)] = (b,b1)
            return (its[1], its[0]), ss

    @staticmethod
    def load_snp_in_read2(fname):

        rinfos = defaultdict(list)

        buffer = []
        pool = mp.Pool(processes=24)


        for line in open(fname):
            its = line.split()
            buffer.append(line)
            # if len(buffer) == 24:
            #     result = pool.map(PrjFile.process_line_snp_in_read, buffer)
            #     buffer = []

            #     for r in result:
            #         rinfos[r[0]].append(r[1])

        if buffer:
            result = pool.map(PrjFile.process_line_snp_in_read, buffer)


            for r in result:
                rinfos[r[0][0]].append((r[0][1], r[1]))

        return rinfos

class PurgeInconsistent(object):
    def __init__(self, tile, readsnps):
        self.readsnps = defaultdict(list)

        item_start = 4
        
        logger.info("Loading reads snps")
        for line in open(readsnps):
            its = line.strip().split()
            assert len(its) >= item_start
            ctg = its[0]
            snps = []
            for i in its[item_start:]:
                p, v, *_ = i.split("|")
                snps.append((p, v))
            self.readsnps[its[1]].append((ctg, snps))
        logger.info("Loading reads in contigs")
        ctgs = PrjFile.get_contig_reads(tile)

        logger.info("Collecting snps in contigs")
        self.ctgsnps = {}
        for c, r in ctgs.items():
            infos = defaultdict(lambda :defaultdict(lambda: defaultdict(int)))
            for ir in set(r):
                if ir in self.readsnps:
                    for rctg, rsnps in self.readsnps[ir]:

                        for p, v in rsnps:
                            infos[rctg][p][v] += 1
            self.ctgsnps[c] = infos

    def __call__(self, rd, ctg):
        if rd in self.readsnps and ctg in self.ctgsnps:
            csnps = self.ctgsnps[ctg]
            count = [0, 0]
            for rctg, rsnps in self.readsnps[rd]:
                if rctg in csnps:
                    for p, v in rsnps:
                        if p in csnps[rctg]:
                            if v in csnps[rctg][p]:
                                
                                if csnps[rctg][p][v] > sum(csnps[rctg][p].values()) / 2:
                                    count[0] += 1
                                elif csnps[rctg][p][v] < sum(csnps[rctg][p].values()) / 2:
                                    count[1] += 1
                            else:
                                count[1] += 1
            #print(count, its[0], its[2])
            return count[0] < count[1]
              
        else:
            return False


class PurgeInconsistentInSam(object):
    #def __init__(self, tile, readsnps, id2name):
    def __init__(self, tile, readsnps):
        self.purge = PurgeInconsistent(tile, readsnps)
        #self.n2id = PrjFile.load_id2name(id2name, True)

    def __call__(self, line):
        its = line.split()
        return self.purge(its[0], its[2])
        #if its[0] in self.n2id:
            #return self.purge(self.n2id[its[0]], its[2])
            
        #else:
        #    return False

class PurgeInconsistentInPaf(object):
    def __init__(self, tile, readsnps):
        self.purge = PurgeInconsistent(tile, readsnps)

    def __call__(self, line):
        its = line.split()
        return self.purge(its[0], its[5])



if __name__ == "__main__":
    import sys
    locals()[sys.argv[1]](*sys.argv[2:]) 
