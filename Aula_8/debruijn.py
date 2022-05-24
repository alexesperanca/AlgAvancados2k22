# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class DeBruijnGraph (MyGraph):
    
    def __init__(self, frags:list):
        MyGraph.__init__(self, {})
        self.create_deBruijn_graph(frags)

    def add_edge(self, o, d):
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        self.graph[o][d] = None

    def in_degree(self, v):
        res = 0
        for k in self.graph.keys():
            if v in self.graph[k]:
                res += list(self.graph[k].keys()).count(v)
        return res

    def create_deBruijn_graph(self, frags:list):        
        for seq in frags:
            pref = prefix(seq)
            self.add_vertex(pref)
            suf = suffix(seq)
            self.add_vertex(suf)
            self.add_edge(pref, suf)

    def seq_from_path(self, path:list) -> str:
        seq = path[0]
        for i in range(1,len(path)):
            nxt = path[i]
            seq += nxt[-1]
        return seq 
    
def suffix(seq:str) -> str: 
    return seq[1:]
    
def prefix(seq:str) -> str:
    return seq[:-1]

def composition(k:int, seq:str) -> list:
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
    res.sort()
    return res

def test1():
    print("* Test 1 *\n")
    frags = ["ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    
def test2():
    print("\n* Test 2 *\n")
    frags = ["ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    print (dbgr.check_nearly_balanced_graph())
    print (dbgr.eulerian_path())

def test3():
    print("\n* Test 3 *\n")
    orig_sequence = "ATGCAATGGTCTG"
    frags = composition(3, orig_sequence)
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    print (dbgr.check_nearly_balanced_graph())
    p = dbgr.eulerian_path()
    print (p)
    print (dbgr.seq_from_path(p))

def test4():
    gr = MyGraph( {1:[2], 2:[3,1], 3:[4], 4:[2,5],
    5:[6], 6:[4]} )
    gr.print_graph()
    print (gr.check_balanced_graph() )
    print (gr.eulerian_cycle() )

test1()
test2()
test4()
test3()