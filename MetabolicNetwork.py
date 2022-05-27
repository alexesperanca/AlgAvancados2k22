# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`MetabolicNetwork` class.
This class includes diverse strategies for the construction of graphs based in metabolites and reactions and information retained, such as:
    - Loading of external files
    - Addition of vertex to the graph
    - Calculate centrality values
    - Obtain metabolites produced
    - Obtain active reactions
"""

from MyGraph import MyGraph

class MetabolicNetwork (MyGraph):
    '''Subclass of MyGraph that creates a graph with metabolites and reactions indicating the type of network construction

    Parameters
    ----------
    MyGraph : class
        Parent class of MetabolicNetwork
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
        if nodetype in self.node_types:
            return sorted(self.node_types[nodetype])
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

        rf.close()
        
        
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
        '''Method calculates the degree centrality of all the nodes in the graph. Node centrality represents the fraction of nodes that is connected to

        Returns
        -------
        dict
            Dictionary of nodes and respective degree centrality
        '''
        res = {}
        deg = self.all_degrees()
        for k, v in deg.items():
            res[k] = v/(len(deg)-1)
        return dict(sorted(res.items(), key = lambda x:x[0]))

    def betweenness_centrality(self, s: str, t: str) -> dict:
        '''Method computes the betweenness centrality of all the nodes of a path between 2 given nodes. Betweenness centrality represents the centrality value of each node in the path of the 2 points

        Parameters
        ----------
        s : str
            Node path starts
        t : str
            Node path ends

        Returns
        -------
        dict
            Dictionary with all the nodes and respective betweenness centrality value 
        '''
        path = self.get_path(s, t)                  # Only the best path is returned, should we calculate to all the paths?
        if path == None: return "No path found"
        res = {}
        for elem in self.graph:
            if elem in path and elem not in (s, t):
                res[elem] = 1
            else: res[elem] = 0
        return dict(sorted(res.items(), key = lambda x:x[0]))

    def closeness_centrality(self) -> dict:
        '''Method calculates the closeness centrality of all the nodes in the graph. Bigger value indicates that a node is more central in the graph

        Returns
        -------
        dict
            Dictionary with all the nodes and respective closeness centrality value
        '''
        res = {}
        deg = self.all_degrees()
        for k, v in deg.items():
            res[k] = (len(self.graph[k]) / (len(deg) - 1)) * (len(self.graph[k]) / v) # Aqui deveria ser a soma das distâncias entre os nodos em vez de v, mas como estamos num grafo sem pesos assumi que tudo é =1
        return dict(sorted(res.items(), key = lambda x:x[0]))

# Portfólio

    def _prod(self, init: list, f: str) -> list:
        '''Auxiliary method of "active_reactions" and "met_prod" functions. Verifies for the given parameter "f" which reactions or metabolites are possibly produced in the system (graph)

        Parameters
        ----------
        init : list
            List of reaction/metabolites initialy active/present.
        f : str
            Indicates if the algorithm should filter reactions or metabolites

        Returns
        -------
        list
            List of the active reactions/metabolites
        '''
        assert self.net_type == "metabolite-reaction", "No Metabolites or Reactions represented in the system"
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
        '''Method that depends on "_prod" function to obtain the active reactions for a given list of initial metabolites

        Parameters
        ----------
        met : list
            List of metabolites

        Returns
        -------
        list
            List of active reactions
        '''
        return self._prod(met, "R")

    def met_prod(self, reac: list) -> list:
        '''Method that depends on "_prod" function to obtain the possibly produced metabolites for a given list of active reactions

        Parameters
        ----------
        reac : list
            List of reactions

        Returns
        -------
        list
            List of metabolites possibly produced
        '''
        return self._prod(reac, "M")

    def final_met(self, init_met: list) -> list:
        '''Method that obtains the possibly produced metabolites for a given list of initial metabolites

        Parameters
        ----------
        init_met : list
            List of metabolites for

        Returns
        -------
        list
            List of final possible produced metabolites
        '''
        assert self.net_type == "metabolite-reaction", "No Metabolites or Reactions represented in the system"
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
