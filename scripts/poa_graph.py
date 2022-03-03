import sys
import traceback
from functools import reduce



class POAGraph:
    class Link(object):
        def __init__(self, curr, prev, count, seqs):
            self.curr = curr
            self.prev = prev
            self.count = count
            self.seqs = set([i for i, b in enumerate(seqs) if b == '1'])

        def __str__(self):
            return str(self.prev) + "->" + str(self.curr)+"|"+str(len(self.seqs))

        def __repr__(self):
            return self.__str__()

    class Node(object):
        def __init__(self, b, col, row):
            self.links = []
            self.base = b
            self.col = col
            self.row = row

        def add_link(self, link):
            self.links.append(link)

        def link_size(self):
            return len(self.links)

        def get_link(self, i):
            return self.links[i]

        def __str__(self):
            return "%d_%d_%s" % (self.col, self.row, self.base)

        def __repr__(self):
            return self.__str__()

    
    class NodeGroup(object):
        def __init__(self, col, row):
            self.col = col
            self.row = row
            self.bases = [POAGraph.Node('A', col, row), POAGraph.Node('C', col, row), POAGraph.Node('G', col, row), POAGraph.Node('T', col, row), POAGraph.Node('-', col, row)]

        def get_base(self, i):
            if type(i) == int:
                assert i >= 0 and i < 5
                return self.bases[i]
            elif type(i) == str:
                b2i = {'A':0, 'C':1, 'G':2, 'T':3, '-':4};
                return self.bases[b2i[i]]

        def get_links(self):
            links = []
            for b in self.bases:
                links.extend(b.links)
            return links

        

    class Column(object):
        def __init__(self, col):
            self.col = col
            self.rows = [POAGraph.NodeGroup(col, 0)]

        def get_row(self, i):
            while len(self.rows) < i+1:
                self.rows.append(POAGraph.NodeGroup(self.col, len(self.rows)));
            return self.rows[i]

        def get_links(self):
            return self.rows[0].get_links()


    def __init__(self):
        self.cols = []

    def load(self, fname):
        for i, line in enumerate(open(fname)):
            if i == 0: continue

            its = line.split(",")
            col, row, base = its[0].split('_')
            pcol, prow, pbase = its[1].split('_')
            count = its[2]
            seqs = its[3]

            self.add_link((int(col), int(row), base), (int(pcol), int(prow), pbase), int(count), seqs)

    def add_link(self, curr, prev, count, seqs):
        curr = self.get_col(curr[0]).get_row(curr[1]).get_base(curr[2])
        prev = self.get_col(prev[0]).get_row(prev[1]).get_base(".ACGT-".index(prev[2])-1) if prev[0] > 0 else None

        curr.add_link(POAGraph.Link(curr, prev, count, seqs))

    def get_col(self, i):
        while len(self.cols) < i+1:
            self.cols.append(POAGraph.Column(len(self.cols)));
        return self.cols[i]

    def cross(self):
        cross_count = 0
        for i, c in enumerate(self.cols):
            links = c.get_links()
            min_branch_coverage = 42*2/3
            total = reduce(lambda x,y:x+y.count, links, 0)
            if len(links) >= 2 and total >= min_branch_coverage * 2:
                links.sort(key=lambda x: -x.count)

                sz = 0
                while sz < len(links):
                    if links[sz].count < min_branch_coverage or  links[sz].count / total < 0.2:
                        break

                    sz += 1

                if sz >= 2:
                    print(i, total, links[0].count, links[1].count)
                    cross_count += 1

        print(cross_count)

    def dump_segments(self, start, sb, end, eb):
        paths = []

        curr_node = self.get_col(end).get_row(0).get_base(eb)
        curr = [curr_node, 0]

        stack = [curr]
        path = []
        while len(stack) > 0:
            #print(stack, stack[-1][1] , stack[-1][0].link_size())
            if stack[-1][1] < stack[-1][0].link_size():
                path.append(stack[-1][0].get_link(stack[-1][1]))
                stack[-1][1] += 1
                next = [path[-1].prev, 0]
                if next[0].col > start or next[0].col == start and next[0].row > 0:
                    stack.append(next)
                else:
                    if next[0].col == start and next[0].row == 0 and next[0].base == sb:
                        paths.append(path.copy())
                    path.pop()
            else:
                stack.pop()
                if len(path) > 0: path.pop()

        for p in paths:
            total = reduce(lambda x,y:x.intersection(y.seqs), p, p[0].seqs)
            print("ooo", len(total))

            for s in p:
                print(sorted(s.seqs))

       
            print("xxx", p)
        



def poa_segments(argv):
    try:
        fname = argv[0]
        start = int(argv[1])
        sbase = argv[2]
        end   = int(argv[3])
        ebase = argv[4]

        poa = POAGraph()
        poa.load(fname)

        poa.dump_segments(start, sbase, end, ebase)
    except:
        traceback.print_exc()


def poa_togv(argv):
    
    try:
        toACGT = ['A','C','G','T','N']

        toACGTN = {'A':'A', 'C':'C', 'G':'G', 'T':'T', '-':'N' }
        ifname = argv[0]
        ofname = argv[1]
        start = int(argv[2])
        end = int(argv[3])
        threshold = int(argv[4])

        ofile = open(ofname, "w")
        ofile.write("""digraph st2 {
    fontname = "Verdana";
    fontsize = 10;
    rankdir=TB;
    splines="curved"
    node [fontname = "Verdana", fontsize = 10, color="skyblue", shape="record", width="0.5"]; 
    edge [fontname = "Verdana", fontsize = 10, color="crimson", style="solid",splines="curved" ];
    """)
    
        done = set()
        for i, line in enumerate(open(ifname)):
            if line.startswith("Target"): continue

            items = [i.strip() for i in line.split(",")]
            curr = items[0].split("_")
            curr[0] = int(curr[0])
            curr[1] = int(curr[1])
            if curr[2] == '-': curr[2] = 'N' 
            prev = items[1].split("_")
            prev[0] = int(prev[0])
            prev[1] = int(prev[1])
            if prev[2] == '-': prev[2] = 'N' 

            weight = items[2]
            if int(weight) <= threshold: continue

            if start <= prev[0]  and curr[0] < end:
                print(done)
                if (curr[0], curr[1]) not in done:
                    ofile.write("""    p%s [label="{<A>A|<C>C|<G>G|<T>T|<N>-|%s}" pos="%d,%d!"]; \n""" % (items[0][:-2], items[0][:-2], (int(curr[0]))*3, int(curr[1])*3))
                    done.add((curr[0], curr[1]))

                if (prev[0], prev[1]) not in done:
                    ofile.write("""    p%s [label="{<A>A|<C>C|<G>G|<T>T|<N>-|%s}" pos="%d,%d!"]; \n""" % (items[1][:-2], items[1][:-2], (int(prev[0]))*3, int(prev[1])*3))
                    done.add((prev[0], prev[1]))

                
                ofile.write("""    p%s:%s -> p%s:%s [label=%s];\n""" %(items[0][:-2], curr[2], items[1][:-2], prev[2], weight))

    
        ofile.write("""}""")

    except:
        traceback.print_exc()

    


def poa_test(argv):
    try:
        fname = argv[0]
        ag = POAGraph()
        ag.load(fname);
        ag.cross()

    except:
        traceback.print_exc()



if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1].startswith("poa_"):
       locals()[sys.argv[1]](sys.argv[2:])
    else:
       for func in list(locals().keys()):
           if func.startswith("poa_"):
               print(func)
