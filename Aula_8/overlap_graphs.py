# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class OverlapGraph(MyGraph):
    
#    def __init__(self, frags):
#        MyGraph.__init__(self, {})
#        self.create_overlap_graph(frags)

    def __init__(self, frags:list, reps = False):
        MyGraph.__init__(self, {})
        if reps: self.create_overlap_graph_with_reps(frags)
        else: self.create_overlap_graph(frags)
        self.reps = reps
        
    
    ## create overlap graph from list of sequences (fragments)
    def create_overlap_graph(self, frags:list):
        for s in frags:
            self.add_vertex(s)
        
        ## add edges
        for seq in frags:
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    self.add_edge(seq, seq2)
        
    def create_overlap_graph_with_reps(self, frags:list):  # caso de replicas de fragmentos
        n = 1
        for s in frags:
            self.add_vertex(f"{s}-{n}")
            n += 1
        n = 1
        for seq in frags:           # Verifica cada elemento do frag com os outros se têm ligação e faz os nodos entre eles basicamente
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    for x in self.get_instances(seq2):
                        self.add_edge(seq+ "-" + str(n), x)
            n += 1
    
    def get_instances(self, seq:str) -> list:
        res = []
        for k in self.graph.keys():
            if seq in k: res.append(k)
        return res
    
    def get_seq(self, node):
        if node not in self.graph.keys(): return None
        if self.reps: return node.split("-")[0]
        else: return node
    
    def seq_from_path(self, path):
        # ...
        return seq    
   
                    
# auxiliary
def composition(k, seq):
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
    res.sort()
    return res
    
def suffix (seq): 
    return seq[1:]
    
def prefix(seq):
    return seq[:-1]

  
# testing / mains
def test1():
    print("* Test 1 *\n")
    seq = "CAATCATGATG"
    k = 3
    print (composition(k, seq))
   
def test2():
    print("\n* Test 2 *\n")
    frags = ["ACC", "ATA", "CAT", "CCA", "TAA"]
    ovgr = OverlapGraph(frags, False)
    ovgr.print_graph()

def test3():
    print("\n* Test 3 *\n")
    frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    ovgr = OverlapGraph(frags, True)
    ovgr.print_graph()

def test4():
    print("\n* Test 4 *\n")
    frags = ["ATA",  "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA" , "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    ovgr = OverlapGraph(frags, True)
    path = ['ACC−2', 'CCA−8', 'CAT−5', 'ATG−3']
    print (ovgr.check_if_valid_path(path))
    print (ovgr.check_if_hamiltonian_path(path))
    path2 = ['ACC−2', 'CCA−8', 'CAT−5', 'ATG−3', 'TGG−13', 'GGC−10', 'GCA−9', 'CAT−6', 'ATT−4', 'TTT−15', 'TTC−14', 'TCA−12', 'CAT−7', 'ATA−1', 'TAA−11']
    print (ovgr.check_if_valid_path(path2))
    print (ovgr.check_if_hamiltonian_path(path2))
    print (ovgr.seq_from_path(path2))

def test5():
    print("\n* Test 5 *\n")
    frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    ovgr = OverlapGraph(frags, True)

    path = ovgr.search_hamiltonian_path()
    print(path)
    print (ovgr.check_if_hamiltonian_path(path))
    print (ovgr.seq_from_path(path))

def test6():
    print("\n* Test 6 *\n")
    orig_sequence = "CAATCATGATGATGATC"
    frags = composition(3, orig_sequence)
    print (frags)
    ovgr = OverlapGraph(frags, True)
    ovgr.print_graph()
    path = ovgr.search_hamiltonian_path()
    print (path)
    print (ovgr.seq_from_path(path))
   
test1()
test2()
test3()
test4()
#test5()
#test6()
