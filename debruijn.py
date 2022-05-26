# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`DeBruijnGraph` class.
This class includes diverse strategies for deBrujin graph construction to more easily obtain and deconstruct sequences, such as:
    - Graph construction of edges between prefixes and suffixes
    - Division of a sequence into fragmented element of a given size
    - Sequence obtainment from an eulerian path acquired from "MyGraph" class

"""

from MyGraph import MyGraph

class DeBruijnGraph (MyGraph):
    '''DeBruijnGraph class depends on its parent class "MyGraph" and executes a graph construction with edges between prefixes and suffixes of fragments. Hereon, if the graph obtained is balanced, eulerian cycles (cycle that passes throw every node and returns to the start) are obtained to generate an eulerian path

    Parameters
    ----------
    MyGraph : class
        Essential class to construct the graph and access to the eulerian methods
    '''
    def __init__(self, frags:list):
        '''DeBruijn graph constructor initialization. Runs the method "create_deBruijn_graph" to obtain the prefix and suffix edges

        Parameters
        ----------
        frags : list
            Elements to form the graph
        '''
        MyGraph.__init__(self, {})
        self.create_deBruijn_graph(frags)

    def add_edge(self, o: str, d: str):
        '''Creation of edges between the two strings/nodes provided

        Parameters
        ----------
        o : str
            Parent node of "d"
        d : str
            Child node of "o"
        '''
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        self.graph[o][d] = None

    def in_degree(self, v: str) -> int:
        '''Method that returns the entry degree of the given node

        Parameters
        ----------
        v : str
            Node to calculate the entry degree

        Returns
        -------
        int
            Value of the entry degree
        '''
        res = 0
        for k in self.graph.keys():
            if v in self.graph[k]:
                res += list(self.graph[k].keys()).count(v)
        return res

    def create_deBruijn_graph(self, frags:list):
        '''Creation of a deBruijn graph where the prefixes of each element of the "frags" has an edge with its suffix

        Parameters
        ----------
        frags : list
            List of fragments of sequence to add to the graph
        '''
        for seq in frags:
            pref = prefix(seq)
            self.add_vertex(pref)
            suf = suffix(seq)
            self.add_vertex(suf)
            self.add_edge(pref, suf)

    def seq_from_path(self, path:list) -> str:
        '''Method that builds the original sequence from the eulerian path provided

        Parameters
        ----------
        path : list
            Eulerian path

        Returns
        -------
        str
            Sequence obtained from the "path"
        '''
        try:
            seq = path[0]
        except:
            return None
        for i in range(1,len(path)):
            nxt = path[i]
            seq += nxt[-1]
        return seq 
    
def suffix(seq:str) -> str: 
    '''Method that obtains the suffix of the "seq" provided

    Parameters
    ----------
    seq : str
        Sequence to obtain the suffix

    Returns
    -------
    str
        Suffix of "seq"
    '''
    return seq[1:]
    
def prefix(seq:str) -> str:
    '''Method that obtains the prefix of the "seq" provided

    Parameters
    ----------
    seq : str
        Sequence to obtain the prefix

    Returns
    -------
    str
        Prefix of "seq"
    '''
    return seq[:-1]

def composition(k:int, seq:str) -> list:
    '''Method that retrieves the sequence "seq" divided into fragmented sequences of "k" size

    Parameters
    ----------
    k : int
        Size of the sequence fragments
    seq : str
        Original sequence

    Returns
    -------
    list
        List of the sequence fragments of "k" size
    '''
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
    frags = composition(4, orig_sequence)
    print("Frags", frags)
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    print (dbgr.check_nearly_balanced_graph())
    p = dbgr.eulerian_path()
    print (p)
    print (dbgr.seq_from_path(p))

def test4():
    print("\n* Test 4 *\n")
    gr = MyGraph( {1:[2], 2:[3,1], 3:[4], 4:[2,5],
    5:[6], 6:[4]} )
    gr.print_graph()
    print (gr.check_balanced_graph())
    print (gr.eulerian_cycle())

if __name__ == "__main__":
    test1()
    test2()
    test3()
    test4()