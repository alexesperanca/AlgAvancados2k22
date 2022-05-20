# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class DeBruijnGraph (MyGraph):
    
    def __init__(self, frags):
        MyGraph.__init__(self, {})
        self.create_deBruijn_graph(frags) # Maybe senão for em frags e sim uma string, chamar composition para originar os frags

    def add_edge(self, o, d):
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        self.graph[o].append(d)

    def in_degree(self, v):
        res = 0
        for k in self.graph.keys(): 
            if v in self.graph[k]: 
                res += self.graph[k].count(v)
        return res

    def create_deBruijn_graph(self, frags):        
        pass

    def seq_from_path(self, path):
        seq = path[0]
        for i in range(1,len(path)):
            nxt = path[i]
            seq += nxt[-1]
        return seq 
    
def suffix (seq): 
    return seq[1:]
    
def prefix(seq):
    return seq[:-1]

def composition(k, seq):
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
    res.sort()
    return res


def test1():
    frags = ["ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    
def test2():
frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
dbgr = DeBruijnGraph(frags)
dbgr.print_graph()
print (dbgr.check_nearly_balanced_graph())
print (dbgr.eulerian_path())


def test3():
    orig_sequence = "ATGCAATGGTCTG"
    frags = composition(3, orig_sequence)
    # ... completar


test1()
print()
#test2()
#print()
#test3()
    
