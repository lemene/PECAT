#!/usr/bin/env python3
import sys

from Bio import SeqIO


def load_infos(fname):
    infos = set()

    for line in open(fname):
        its = line.split()
        r0, r1, s = int(its[0]), int(its[1]), int(its[-1])
        if s == 1:
            if r0 < r1:
                infos.add((r0, r1))
            else:
                infos.add((r1, r0))

    return infos

def load_read(fname):
    rd = {}
    for rec in SeqIO.parse(fname, "fasta"):
        s, e = [int(i) for i in rec.description.split()[1].split("=")[1].split("-") ]
        rd[int(rec.id)] = (s, e, len(rec.seq))

    return rd
    
def adjust_filtered_overlaps(flt_ols_fname, adj_ols_fname, rd, infos):

    def map_position(s0, e0, s1, e1, p0):
        return int((e1 - s1) / (e0 - s0) * (p0 - s0) + s1)
    
    with open(adj_ols_fname, "w") as f:
        done = set()
        for line in open(flt_ols_fname):
            #print(line, end="")
            its = line.split()
            qid, qlen, qstart, qend = [int(i) for i in its[0:4]]
            tid, tlen, tstart, tend = [int(i) for i in its[5:9]]
            d = its[4]
            qual = int(its[9]) / int(its[10])

            
            tqid = (qid, tid) if qid < tid else (tid, qid)

            if tqid in infos and qid in rd and tid in rd:
                if tqid not in done:
                    done.add(tqid)
    
                    qs, qe, ql = rd[qid]
                    ts, te, tl = rd[tid]

                    #print("rd:", qs, qe, ql, ts, te, tl)

                    if d == '+':
                        qs_t = map_position(qstart, qend, tstart, tend, qs)
                        qe_t = map_position(qstart, qend, tstart, tend, qe)

                        adj_tstart = max(qs_t, ts)
                        adj_tend = min(qe_t, te)

                        ts_q = map_position(tstart, tend, qstart, qend, ts)
                        te_q = map_position(tstart, tend, qstart, qend, te)
                        adj_qstart = max(ts_q, qs)
                        adj_qend = min(te_q, qe)

                        #print("adj:", adj_qstart, adj_qend, adj_tstart, adj_tend)

                        if adj_qend > adj_qstart + 1000 and adj_tend > adj_tstart + 1000:

                            qrate = ql / (qe - qs)
                            adj_qlen = ql
                            adj_qstart = int((adj_qstart - qs) * qrate)
                            adj_qend = int((adj_qend - qs) * qrate)

                            trate = tl / (te - ts)
                            adj_tlen = tl
                            adj_tstart = int((adj_tstart - ts) * trate)
                            adj_tend = int((adj_tend - ts) * trate)

                            ol_len = (adj_qend - adj_qstart + adj_tend - adj_tstart) // 2
                            match = int(ol_len * qual)

                            f.write(f"{qid}\t{adj_qlen}\t{adj_qstart}\t{adj_qend}\t{d}\t{tid}\t{adj_tlen}\t{adj_tstart}\t{adj_tend}\t{ol_len}\t{match}\t{its[-1]}\n")

                    else:
                        assert d == '-'                    
                        qs_t = map_position(qstart, qend, tend, tstart, qs)
                        qe_t = map_position(qstart, qend, tend, tstart, qe)

                        adj_tstart = max(qe_t, ts)
                        adj_tend = min(qs_t, te)

                        ts_q = map_position(tstart, tend, qend, qstart, ts)
                        te_q = map_position(tstart, tend, qend, qstart, te)
                        adj_qstart = max(te_q, qs)
                        adj_qend = min(ts_q, qe)

                        #print("adj:", adj_qstart, adj_qend, adj_tstart, adj_tend)

                        if adj_qend > adj_qstart + 1000 and adj_tend > adj_tstart + 1000:

                            qrate = ql / (qe - qs)
                            adj_qlen = ql
                            adj_qstart = int((adj_qstart - qs) * qrate)
                            adj_qend = int((adj_qend - qs) * qrate)

                            trate = tl / (te - ts)
                            adj_tlen = tl
                            adj_tstart = int((adj_tstart - ts) * trate)
                            adj_tend = int((adj_tend - ts) * trate)

                            ol_len = (adj_qend - adj_qstart + adj_tend - adj_tstart) // 2
                            match = int(ol_len * qual)

                            f.write(f"{qid}\t{adj_qlen}\t{adj_qstart}\t{adj_qend}\t{d}\t{tid}\t{adj_tlen}\t{adj_tstart}\t{adj_tend}\t{match}\t{ol_len}\t{its[-1]}\n")
                else:
                    done.remove(tqid)



if __name__ == "__main__":
    flt_ols_fname = sys.argv[1]
    adj_ols_fname = sys.argv[2]
    rd_fname = sys.argv[3]
    infos_fname = sys.argv[4]

    infos = load_infos(infos_fname)
    rd = load_read(rd_fname)

    adjust_filtered_overlaps(flt_ols_fname, adj_ols_fname, rd, infos)
