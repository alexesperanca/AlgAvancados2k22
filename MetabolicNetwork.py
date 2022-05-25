# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`MetabolicNetwork` class, for construction of graphs, given a dictionary of values or strings.
This class includes diverse strategies for (...), such as:
"""

from sympy import EX
from MyGraph import MyGraph

class MetabolicNetwork (MyGraph):
    '''_summary_

    Parameters
    ----------
    MyGraph : _type_
        _description_
    '''
    
    def __init__(self, network_type: str = "metabolite-reaction", split_rev: bool = False):
        '''Construction of a graph with metabolites and/or reactions provided to the class

        Parameters
        ----------
        network_type : str, optional
            Informs the type of graph network the algorithm is working with, by default "metabolite-reaction"
            "metabolite-reaction": Graph with metabolites and reactions occuring 
            "metabolite-metabolite": Graph with the metabolite that origin others
            "reaction-reaction": Graph with all the reactions that induce others to happen
        split_rev : bool, optional
            Indicates if the reactions are considered reversible or not, by default False
        '''
        MyGraph.__init__(self, {})
        self.net_type = network_type
        self.node_types = {}
        if network_type == "metabolite-reaction":
            self.node_types["metabolite"] = []
            self.node_types["reaction"] = []
        self.split_rev = split_rev
    
    def add_vertex_type(self, v: str, nodetype: str):
        '''Addition of vertex to the graph constructed

        Parameters
        ----------
        v : str
            Vertex to be added
        nodetype : str
            Vertex type, only "metabolite" or "reaction" supported
        '''
        assert nodetype in ["metabolite", "reaction"], "Node type given not supported"
        self.add_vertex(v)
        self.node_types[nodetype].append(v)
    
    def get_nodes_type(self, nodetype: str) -> list:
        '''Obtain the nodes of a given node type

        Parameters
        ----------
        nodetype : str
            Node type, only "metabolite" or "reaction" supported

        Returns
        -------
        list
            List of nodes in the type provided
        '''
        assert nodetype in ["metabolite", "reaction"], "Nodetype given not supported"
        if nodetype in self.nodetypes:
            return self.nodetypes[nodetype]
        else: return None
    
    def load_from_file(self, filename: str):
        '''Loading of a file content to graph

        Parameters
        ----------
        filename : str
            File name with the content desired to be loaded to a graph. Restrictive content exhibition is required

        Raises
        ------
        Exception
            Invalid line exhibited in the file loaded
        Exception
            Invalid line exhibited in the file loaded
        '''
        rf = open(filename)
        gmr = MetabolicNetwork("metabolite-reaction")
        for line in rf:
            if ":" in line:
                tokens = line.split(":")
                reac_id = tokens[0].strip()
                gmr.add_vertex_type(reac_id, "reaction")
                rline = tokens[1]
            else: raise Exception("Invalid line:")                
            if "<=>" in rline:
                left, right = rline.split("<=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_vertex_type(reac_id+"_b", "reaction")
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id+"_b", met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_edge(met_id, reac_id+"_b")
                        gmr.add_edge(reac_id, met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
            elif "=>" in line:
                left, right = rline.split("=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(met_id, reac_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(reac_id, met_id)
            else: raise Exception("Invalid line:")    

        
        if self.net_type == "metabolite-reaction": 
            self.graph = gmr.graph
            self.node_types = gmr.node_types
        elif self.net_type == "metabolite-metabolite":
            self.convert_metabolite_net(gmr)
        elif self.net_type == "reaction-reaction": 
            self.convert_reaction_graph(gmr)
        else: self.graph = {}
        
        
    def convert_metabolite_net(self, gmr: dict):
        '''Addition to the main graph the metabolites that possess an relation, information extracted from the loaded file

        Parameters
        ----------
        gmr : dict
            Current graph
        '''
        for m in gmr.node_types["metabolite"]: # Acede a cada metabolito obtido do ficheiro
            self.add_vertex(m)
            sucs = gmr.get_successors(m)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if m != s2: self.add_edge(m, s2)
        
    def convert_reaction_graph(self, gmr: dict):
        '''Addition to the main graph the reactions that possess an relation, information extracted from the loaded file


        Parameters
        ----------
        gmr : dict
            Current graph
        '''
        for r in gmr.node_types["reaction"]:
            self.add_vertex(r)
            sucs = gmr.get_successors(r)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if r != s2: self.add_edge(r, s2)

    def degrees_centrality(self) -> dict:
        res = {}
        deg = self.all_degrees()
        for k, v in deg.items():
            res[k] = v/(len(deg)-1)
        return res

    def betweenness_centrality(self, s, t):
        '''
        Calcula a centralidade de todos os nodos restantes do caminho entre 2 nodos
        É suposto aqui haver vários caminhos, no entanto, apenas retornamos na outra função o melhor.. Como fazer?
        Será que verifico apenas se passa?
        '''
        path = self.get_path(s, t)
        if path == None: return "No path found"
        res = {}
        for elem in self.graph:
            if elem in path and elem not in (s, t):
                res[elem] = 1
            else: res[elem] = 0
        return res

    def closeness_centrality(self):
        '''Daria o mesmo que o degrees_centrality caso consideremos grafo não-direcionado, o que não acontece aqui'''
        res = {}
        deg = self.all_degrees()
        for k, v in deg.items():
            res[k] = (len(self.graph[k]) / (len(deg) - 1)) * (len(self.graph[k]) / v) # Aqui deveria ser a soma das distâncias entre os nodos em vez de v, mas como estamos num grafo sem pesos assumi que tudo é =1
        return res

# Exs Portfólio do PP

    def _prod(self, init, f):
        '''Função auxiliar para as 2 primeiras do porfólio já que todas seguem o mesmo padrão (acho eu)'''
        if self.net_type != "metabolite-reaction": raise Exception("No Metabolites or Reactions represented in the system")
        n = init.copy()
        visited = []
        res = []
        while len(n) != 0:
            m = n.pop(0)
            for j in self.graph[m]:
                if f in j and j not in res: 
                    res.append(j)
                if j not in visited:
                    visited.append(j)
                    n.append(j)
        return res

    def active_reactions(self, met: list) -> list:
        '''_summary_

        Parameters
        ----------
        met : list
            _description_

        Returns
        -------
        list
            _description_
        '''
        return self._prod(met, "R")

    def met_prod(self, reac):
        '''Metabolitos que poderão ser produzidos dada uma lista de reações'''
        return self._prod(reac, "M")

    def final_met(self, init_met):
        '''Metabolitos finais produzidos, tendo em conta uma lista inicial de metabolitos'''
        if self.net_type != "metabolite-reaction": raise Exception("No Metabolites or Reactions represented in the system")
        visited = []
        it = init_met.copy()
        res = []
        while len(it) != 0:
            x = it.pop(0)
            reac = self.active_reactions([x])
            met = self.met_prod(reac)
            for m in met:
                if m not in visited:
                    visited.append(m)
                    it.append(m)
                    res.append(m)
        return res

def test1():
    print("* Test 1 *\n")
    m = MetabolicNetwork("metabolite-reaction")
    m.add_vertex_type("R1","reaction")
    m.add_vertex_type("R2","reaction")
    m.add_vertex_type("R3","reaction")
    m.add_vertex_type("M1","metabolite")
    m.add_vertex_type("M2","metabolite")
    m.add_vertex_type("M3","metabolite")
    m.add_vertex_type("M4","metabolite")
    m.add_vertex_type("M5","metabolite")
    m.add_vertex_type("M6","metabolite")
    m.add_edge("M1","R1")
    m.add_edge("M2","R1")
    m.add_edge("R1","M3")
    m.add_edge("R1","M4")
    m.add_edge("M4","R2")
    m.add_edge("M6","R2")
    m.add_edge("R2","M3")
    m.add_edge("M4","R3")
    m.add_edge("M5","R3")
    m.add_edge("R3","M6")
    m.add_edge("R3","M4")
    m.add_edge("R3","M5")
    m.add_edge("M6","R3")
    m.print_graph()
    print("Reactions: ", m.get_nodes_type("reaction") )
    print("Metabolites: ", m.get_nodes_type("metabolite") )

        
def test2():
    print("\n* Test 2 *\n")
    print("metabolite-reaction network:")
    mrn = MetabolicNetwork("metabolite-reaction")
    mrn.load_from_file("example-net.txt")
    mrn.print_graph()
    print("Reactions: ", mrn.get_nodes_type("reaction") )
    print("Metabolites: ", mrn.get_nodes_type("metabolite") )
    print()
    
    print("metabolite-metabolite network:")
    mmn = MetabolicNetwork("metabolite-metabolite")
    mmn.load_from_file("example-net.txt")
    mmn.print_graph()
    print()
    
    print("reaction-reaction network:")
    rrn = MetabolicNetwork("reaction-reaction")
    rrn.load_from_file("example-net.txt")
    rrn.print_graph()
    print()
    
    print("metabolite-reaction network (splitting reversible):")
    mrsn = MetabolicNetwork("metabolite-reaction", True)
    mrsn.load_from_file("example-net.txt")
    mrsn.print_graph()
    print()
    
    print("reaction-reaction network (splitting reversible):")
    rrsn = MetabolicNetwork("reaction-reaction", True)
    rrsn.load_from_file("example-net.txt")
    rrsn.print_graph()
    print()

    print("* TESTES *")
    mmn.print_graph()
    print(mmn.degrees_centrality())
    print(mmn.closeness_centrality())
    print(mmn.betweenness_centrality('M6', 'M2'))

    print("\n* Portfólio *\n")
    mrn.print_graph()
    print(mrn.active_reactions(["M3", "M2"]))
    # print(rrsn.active_reactions(["M1", "M2"])) -> Funciona bem, levanta exceção como suposto
    # print(mmn.active_reactions(["M1", "M2"])) -> Funciona bem, levanta exceção como suposto
    print(mrn.met_prod(["M4"]))
    print(mrn.final_met(["M5", "M2"])) # Not sure se isto está correto pq estou a considerar ciclos (ou seja, um metabolito produz outro e este porduz o anterior), mas pode n ser o caso

if __name__ == "__main__":
    test1()
    test2()