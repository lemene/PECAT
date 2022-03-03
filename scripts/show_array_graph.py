import sys


toACGT = ['A','C','G','T','N']

toACGTN = {'A':'A', 'C':'C', 'G':'G', 'T':'T', '-':'N' }

def to_gv(msafile, gvfile, start, end):
    gv = open(gvfile, "w")
    gv.write("""digraph st2 {
    fontname = "Verdana";
    fontsize = 10;
    rankdir=TB;
    splines="curved"
    node [fontname = "Verdana", fontsize = 10, color="skyblue", shape="record", width="0.5"]; 
    edge [fontname = "Verdana", fontsize = 10, color="crimson", style="solid",splines="curved" ]; 

    """)
    
    nodes = set()
    pos = None
    for i, line in enumerate(open(msafile)):
        
        items = line.split()
        if start <= int(items[0]) and int(items[0]) < end:
    
            if pos == None or pos != (items[0], items[1]):
                node = "_".join(["p"]+items[:2])
                gv.write("""    %s [label="{<A>A|<C>C|<G>G|<T>T|<N>-|%s}" pos="%d,%d!"]; \n""" % (node, node, (int(items[0]))*3, int(items[1])*3))
                pos = (items[0], items[1])
                nodes.add(pos)

            if items[5] != '.':
            
                tid = "_".join(["p"]+items[:2])
                sid = "_".join(["p"]+items[3:5])
                if (items[0], items[1]) in nodes and (items[3], items[4]) in nodes:
                    gv.write("""    %s:%s -> %s:%s [label=%s];\n""" %(tid, toACGT[int(items[2])], sid, toACGTN[items[5]], items[6]))

    
    gv.write("""}""")
if __name__ == "__main__":
    to_gv(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
