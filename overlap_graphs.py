# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`OverlapGraph` class.
This class includes diverse strategies that uses Hamiltonian methods for faster sequence analysis, sequence obtainment from graphs, and suffix/prefix iteration such as:
    - Creation of an overlap graph that essentially constructs a graph based on the portion of sequences that overlap each other (prefix and suffix)
    - Get original sequence from a path of fragmented sequences by considering their overlap 
"""

from MyGraph import MyGraph

class OverlapGraph(MyGraph):
    '''DeBruijnGraph class depends on its parent class "MyGraph" and executes a graph construction with edges between overlap of prefixes and suffixes of fragments. Hereon, sequence reconstruction is possible with Hamiltonian method passing through the nodes of the graph once assuming their overlap

    Parameters
    ----------
    MyGraph : Class
        Essential class to construct the graph and access to the hamiltonian methods
    '''
    def __init__(self, frags:list, reps: bool = False):
        '''Construction of a graph from the elements of the list "frags" provided. Graph can accept repeats or not, depending on the user preference

        Parameters
        ----------
        frags : list
            List of elements to create the graph
        reps : bool, optional
            If True the graph will add to each element an identifier to avoid repititions. If False repititions are not considered. By default False
        '''
        MyGraph.__init__(self, {})
        if reps: self.create_overlap_graph_with_reps(frags)
        else: self.create_overlap_graph(frags)
        self.reps = reps

    
    ## create overlap graph from list of sequences (fragments)
    def create_overlap_graph(self, frags: list):
        '''Graph creation without repititions consideration. Relation between nodes is created if overlap of prefix and suffix exists

        Parameters
        ----------
        frags : list
            List of elements to add to the graph
        '''
        for s in frags:
            self.add_vertex(s)
        
        ## add edges
        for seq in frags:
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    self.add_edge(seq, seq2)
        
    def create_overlap_graph_with_reps(self, frags: list):
        '''Graph creation considering repititions. Here each element possesses an identifier as "-" + number. Relation between nodes is created if overlap of prefix and suffix exists

        Parameters
        ----------
        frags : list
            List of elements to add to the graph
        '''
        n = 1
        for s in frags:
            self.add_vertex(f"{s}-{n}")
            n += 1
        n = 1
        for seq in frags:           # Verifica cada elemento do frag com os outros se têm overlaps (prefixo com sufixo) e faz os nodos entre eles basicamente
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    for x in self.get_instances(seq2):
                        self.add_edge(seq+ "-" + str(n), x)
            n += 1
    
    def get_instances(self, seq:str) -> list:
        '''Auxiliary function of "create_overlap_graph" and "create_overlap_graph_with_reps" methods where obtains the elements of the graph where "seq" occurs in the graph

        Parameters
        ----------
        seq : str
            _description_

        Returns
        -------
        list
            _description_
        '''
        res = []
        for k in self.graph.keys():
            if seq in k: res.append(k)
        return res
    
    def get_seq(self, node: str) -> str:
        '''Auxiliary method of "seq_from_path" that obtains the sequence separated from the number identifier. Example: "ACC-1" -> "ACC"

        Parameters
        ----------
        node : str
            Node to obtain the only the sequence

        Returns
        -------
        str
            Final sequence without the identifier
        '''
        if node not in self.graph.keys(): return None
        if self.reps: return node.split("-")[0]
        else: return node
    
    def seq_from_path(self, path: list) -> str:
        '''Obtain sequence from a give path. Essentially obtains the first element of the list and adds to the string the final character of the rest of the elements 

        Parameters
        ----------
        path : list
            List elements of a path

        Returns
        -------
        str
            String with the final sequence obtained from the path
        '''
        if not self.check_if_hamiltonian_path(path): return None
        seq = self.get_seq(path[0])
        for i in range(1, len(path)):
            n = self.get_seq(path[i])
            seq += n[-1]
        return seq
   
                    
# auxiliary
def composition(k: int, seq: str) -> list:
    '''Method to obtain each element of "k" size possible from the sequence "seq"

    Parameters
    ----------
    k : int
        Size of each element to retrieve
    seq : str
        Original sequence

    Returns
    -------
    list
        List of elements from the "seq" with "k" size
    '''
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
    res.sort()
    return res
    
def suffix(seq: str) -> str:
    '''Method to get the suffix of a given sequence

    Parameters
    ----------
    seq : str
        Original sequence

    Returns
    -------
    str
        Suffix from the given sequence
    ''' 
    return seq[1:]
    
def prefix(seq):
    '''Method to get the prefix of a given sequence

    Parameters
    ----------
    seq : str
        Original sequence

    Returns
    -------
    str
        Prefix from the given sequence
    ''' 
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
    ovgr = OverlapGraph(frags, False)
    path = ['acC', 'CCA', 'CAT', 'ATG']
    ovgr.print_graph()
    print (ovgr.check_if_valid_path(path))
    print (ovgr.check_if_hamiltonian_path(path))
    path2 = ['ACC', 'CCA', 'CAT', 'ATG', 'TGG', 'GGC', 'GCA', 'CAT', 'ATT', 'TTT', 'TTC', 'TCA', 'CAT', 'ATA', 'TAA']
    print (ovgr.check_if_valid_path(path2))
    print (ovgr.check_if_hamiltonian_path(path2))
    print (ovgr.seq_from_path(path2))
    print(ovgr.search_hamiltonian_path())


def test5():
    print("\n* Test 5 *\n")
    frags = ["ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    ovgr = OverlapGraph(frags, True)
    ovgr.print_graph()
    path = ovgr.search_hamiltonian_path()
    print (path)
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
    print (ovgr.check_if_hamiltonian_path(path))
    print (ovgr.seq_from_path(path))

if __name__ == "__main__":
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
