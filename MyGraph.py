# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`MyGraph` class, for construction of graphs, given a dictionary of values or strings.
This class includes diverse strategies for graph iteration, path discovery, etc, such as:
    - Incremention of nodes and edges construction
    - Obtain successors, predecessors and adjacent nodes
    - Degrees calculation
    - Reachables nodes obtainment
    - Path discovery between differente locations (nodes)
    - Calculation of clustering coefficients
    - Hamiltonian path discovery
    - Eulerian path discovery
"""

from typing import Union

class MyGraph:
    '''Class MyGraph is a parent class of several different ones such as: MetabolicNetwork, OverlapGraph and DeBruijnGraph. This class essentially builds graphs with a given dictionary of values or strings.
    '''
    def __init__(self, g: dict = {}):
        '''Constrution of the initial graph represented by a dictionary

        Parameters
        ----------
        g : dict, optional
            Dictionary given to represent the graph, by default {}
        '''
        self.graph = g    

    def print_graph(self):
        '''Prints the content of the graph as adjacency list 
        '''
        for v in self.graph.keys():
            print (v, " -> ", self.graph[v])

## Basic info

    def get_nodes(self) -> list:
        '''Obtain list of nodes in the graph

        Returns
        -------
        list
            Nodes of the graph
        '''
        return list(self.graph.keys())
        
    def get_edges(self) -> list: 
        '''Method that gets the edges in the graph as a list of tuples (origin, destination)

        Returns
        -------
        list
            List of tuples with information about the origin and destination
        '''
        edges = []
        for v in self.graph.keys():
            for d in self.graph[v]:
                edges.append((v,d))
        return edges
      
    def size(self) -> tuple:
        '''Method that obtains the size of the graph: number of nodes and number of edges

        Returns
        -------
        tuple
            Respectively, number of nodes and number of edges
        '''
        return len(self.get_nodes()), len(self.get_edges())
      
    ## add nodes and edges    
    
    def add_vertex(self, v: Union[str, int, float]):
        '''Add a vertex to the graph and tests if vertex exists or not, adding if True

        Parameters
        ----------
        v : Union[str, int, float]
            Element to add in the graph
        '''
        if v not in self.graph:
            self.graph[v] = {}
        
    def add_edge(self, o: Union[str, int, float], d: Union[str, int, float]):
        '''Add edge to the graph; if vertices do not exist, they are added to the graph

        Parameters
        ----------
        o : Union[str, int, float]
            Element to add in the beginning of the graph and edge with "d"
        d: Union[str, int, float]
            Element to edge with "o"
        ''' 
        if o not in self.graph:
            self.graph[o] = {d: None}
        else:
            g = self.graph[o]
            g[d] = None

    ## successors, predecessors, adjacent nodes
        
    def get_successors(self, v: Union[str, int, float]) -> list:
        '''Obtain successors of the given element

        Parameters
        ----------
        v : Union[str, int, float]
            Element from which successors are returned

        Returns
        -------
        list
            Successors of the given element
        '''
        return list(self.graph[v].keys())
             
    def get_predecessors(self, v: Union[str, int, float]) -> list:
        '''Obtain predecessors of the given element

        Parameters
        ----------
        v : Union[str, int, float]
            Element from which predecessors are returned

        Returns
        -------
        list
            Predecessors of the given element
        '''
        return [k for k in self.graph.keys() if v in self.graph[k]]
    
    def get_adjacents(self, v: Union[str, int, float]) -> list:
        '''Obtain adjacents of the given element

        Parameters
        ----------
        v : Union[str, int, float]
            Element from which adjacents are obtained

        Returns
        -------
        list
            Adjacents of the given element
        '''
        pred = self.get_predecessors(v)
        suc = self.get_successors(v)
        for i in suc:
            if i not in pred: pred.append(i)
        return pred
        
    ## degrees    
    
    def out_degree(self, v: Union[str, int, float]) -> int:
        '''Obtain the number of out degree. Represents the number of successors/ramifications this node of the graph possesses 

        Parameters
        ----------
        v : Union[str, int, float]
            Node to obtain number of out degree

        Returns
        -------
        int
            Number of successors/ramifications from the given node "v"
        '''
        return len(self.graph[v])
    
    def in_degree(self, v: Union[str, int, float]) -> int:
        '''Obtain number of in degree. Represents the number of predecessors this node of the graph possesses

        Parameters
        ----------
        v : Union[str, int, float]
            Node to obtain number of in degree

        Returns
        -------
        int
            Number of predecessors from the given node "v"
        '''
        return len(self.get_predecessors(v))
        
    def degree(self, v: Union[str, int, float]) -> int:
        '''Obtain the number of degree. Represents the number adjacentes nodes of the given one

        Parameters
        ----------
        v : Union[str, int, float]
            Node to obtain number of degree

        Returns
        -------
        int
            Number of adjacent nodes of the given one "v"
        '''
        return len(self.get_adjacents(v))
        
    def all_degrees(self, deg_type: str = "inout") -> dict:
        '''Computes the degree (of a given type) for all nodes. "Deg_type" can be "in", "out", or "inout"

        Parameters
        ----------
        deg_type : str, optional
            Informs which type of degree to compute and return, by default "inout".
            "inout": Count of the in degree and out degree
            "in": Count of the in degree
            "out": Count of the out degree

        Returns
        -------
        dict
            Dictionary with all the nodes of the graph and respective degree value
        '''
        degs = {}
        for v in self.graph.keys():
            if deg_type == "out" or deg_type == "inout":
                degs[v] = len(self.graph[v])
            else: degs[v] = 0
        if deg_type == "in" or deg_type == "inout":
            for v in self.graph.keys():
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degs[d] = degs[d] + 1
        return degs

    ## BFS and DFS searches    
    
    def reachable_bfs(self, v: Union[str, int, float]) -> list:
        '''Method that obtains the reachable nodes from a given one using the BFS algorithm. BFS algorithm searches each level of the graph before passing to the following one

        Parameters
        ----------
        v : Union[str, int, float]
            Node to search for the reachables ones

        Returns
        -------
        list
            List of reachables nodes
        '''
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res
        
    def reachable_dfs(self, v: Union[str, int, float]) -> list:
        '''Method that obtains the reachable nodes from a given one using the DFS algorithm. DFS algorithm searches all successors (all levels) of an node before passing to the following one

        Parameters
        ----------
        v : Union[str, int, float]
            Node to search for the reachables ones

        Returns
        -------
        list
            List of reachables nodes
        '''
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            s = 0
            for elem in self.graph[node]:
                if elem not in res and elem not in l:
                    l.insert(s, elem)
                    s += 1
        return res    

    def _path(self, graph: dict, s: Union[str, int, float], d: Union[str, int, float]) -> list:
        '''Auxiliary function that obtains the path performed to reach from node "s" to "d"

        Parameters
        ----------
        graph : dict
            Path to analyze
        s : Union[str, int, float]
            Node where the path begins
        d : Union[str, int, float]
            Node where the path ends

        Returns
        -------
        list
            Path performed to reach from "s" to "d"
        '''
        visited = [s]
        path = {s: []}
        tries = 0
        try:
            while tries <= len(graph) and len(visited) != len(graph):
                tries += 1
                node = list(path.keys())[0]
                p = path.pop(node)
                for elem in graph[node]:
                    if elem == d:
                        return p + [node, elem]
                    elif elem not in visited:
                        path[elem] = p + [node]
                        visited.append(elem)
            return None
        except:
            return None

    def get_path(self, s: Union[str, int, float], d: Union[str, int, float]) -> list:
        '''Decision of the best path to perform. Methos executes the path for the normal graph and reverted graph.

        Parameters
        ----------
        s : Union[str, int, float]
            Node where the path begins
        d : Union[str, int, float]
            Node where the path ends

        Returns
        -------
        list
            Path performed to reach from "s" to "d"
        '''
        frontpath = self._path(self.graph, s, d)
        rev_graph = self.revert_graph()
        backpath = self._path(rev_graph, s, d)
        if frontpath == None and backpath == None: return None 
        elif backpath == None and frontpath != None: return frontpath
        elif frontpath == None and backpath != None: return backpath
        elif len(frontpath) <= len(backpath): return frontpath
        else: return backpath

    def revert_graph(self) -> dict:
        '''Reversion of the original graph

        Returns
        -------
        dict
            Reverted graph
        '''
        def add_edge(graph: dict, o: Union[str, int, float], d: Union[str, int, float]) -> dict:
            '''Edge addition to the new graph. Auxillary function to obtain the reverted graph.

            Parameters
            ----------
            graph : dict
                Actual graph to be modified
            o : Union[str, int, float]
                Element to add in the beginning of the graph and edge with "d"
            d : Union[str, int, float]
                Element to edge with "o"

            Returns
            -------
            dict
                Actual graph
            '''
            if o not in graph:
                graph[o] = {d: None}
            else:
                g = graph[o]
                g[d] = None
            return graph

        novo = {}
        for U in self.graph:
            for V in self.graph[U]:
                novo = add_edge(novo, V, U)
        return novo        

    def distance(self, s: Union[str, int, float], d: Union[str, int, float]) -> int:
        '''Method that calculates the distance of the path performed to reach from the node "s" to "d"

        Parameters
        ----------
        s : Union[str, int, float]
            Node where the path begins
        d : Union[str, int, float]
            Node where the path ends

        Returns
        -------
        int
            Distance of the path
        '''
        if s == d: return 0
        path = self.get_path(s, d)
        return len(path) - 1
        
    def shortest_path(self, s: Union[str, int, float], d: Union[str, int, float]) -> str:
        '''Print of the shortest path obtained

        Parameters
        ----------
        s : Union[str, int, float]
            Node where the path begins
        d : Union[str, int, float]
            Node where the path ends

        Returns
        -------
        str
            Shortest path performed
        '''
        if s == d: return [s,d]
        path = self.get_path(s, d)
        nodes = " -> ".join(f"{i}" for i in path)
        return nodes
        
    def reachable_with_dist(self, s: Union[str, int, float]) -> list:
        '''Method that returns an list of the reachable nodes from the given "s" and respective distance needed

        Parameters
        ----------
        s : Union[str, int, float]
            Node to get the reachable ones

        Returns
        -------
        list
            Tuples of the reachable nodes and distance: (Node, Distance)
        '''
        res = []
        l = [(s,0)]
        while len(l) > 0:
            node, dist = l.pop(0)
            if node != s: res.append((node,dist))
            for elem in self.graph[node]:
                if not is_in_tuple_list(l,elem) and not is_in_tuple_list(res,elem): 
                    l.append((elem,dist+1))
        return res

    ## mean distances ignoring unreachable nodes
    def mean_distances(self) -> tuple:
        '''Method to calculate the mean distance of the totality of the reachable nodes

        Returns
        -------
        tuple
            Mean distance of the reachable nodes: (Mean distance, Reachable nodes)
        '''
        tot = 0
        num_reachable = 0
        for k in self.graph.keys(): 
            distsk = self.reachable_with_dist(k)
            for _, dist in distsk:
                tot += dist
            num_reachable += len(distsk)
        meandist = float(tot) / num_reachable
        n = len(self.get_nodes())
        return meandist, float(num_reachable)/((n-1)*n) 

## cycles
    def node_has_cycle (self, v: Union[str, int, float]) -> bool:
        '''Method that verifies if a node has a cycle

        Parameters
        ----------
        v : Union[str, int, float]
            Node from which it is checked if a cycle is detected

        Returns
        -------
        bool
            True if a cycle is detected, False otherwise
        '''
        l = [v]
        res = False
        visited = [v]
        while len(l) > 0:
            node = l.pop(0)
            for elem in self.graph[node]:
                if elem == v: return True
                elif elem not in visited:
                    l.append(elem)
                    visited.append(elem)
        return res

    def has_cycle(self) -> bool:
        '''Method that verifies if the graph has any cycle

        Returns
        -------
        bool
            True if a cycle exists, False otherwise
        '''
        res = False
        for v in self.graph.keys():
            if self.node_has_cycle(v): return True
        return res

## topological metrics over degrees

    def mean_degree(self, deg_type: str = "inout") -> float:
        '''Computes the mean value of the degrees (of a given type) for all nodes. "Deg_type" can be "in", "out", or "inout"

        Parameters
        ----------
        deg_type : str, optional
            Informs which type of degree to compute and return, by default "inout".
            "inout": Count of the in degree and out degree
            "in": Count of the in degree
            "out": Count of the out degree

        Returns
        -------
        float
            Mean value of all the degrees
        '''
        degs = self.all_degrees(deg_type)
        return sum(degs.values()) / float(len(degs))
        
    def prob_degree(self, deg_type: str = "inout") -> dict:
        '''Computes the probable value of each node in the degree confered considering the totality of the degrees. "Deg_type" can be "in", "out", or "inout"

        Parameters
        ----------
        deg_type : str, optional
            Informs which type of degree to compute and return, by default "inout".
            "inout": Count of the in degree and out degree
            "in": Count of the in degree
            "out": Count of the out degree
        Returns
        -------
        dict
            Returns a dictionary with all the nodes and respective probable degree value
        '''
        degs = self.all_degrees(deg_type)
        res = {}
        for k in degs.keys():
            if degs[k] in res.keys():
                res[degs[k]] += 1
            else:
                res[degs[k]] = 1
        for k in res.keys():
            res[k] /= float(len(degs))
        return res   

## clustering
        
    def clustering_coef(self, v: Union[str, int, float]) -> float:
        '''Method calculates the clustering coefficient for a given node. Clustering coefficient represents the tendency that a node has to form a cluster with the adjacents

        Parameters
        ----------
        v : Union[str, int, float]
            Node to calculate the clustering coefficient

        Returns
        -------
        float
            Value of the clustering coefficient
        '''
        adjs = self.get_adjacents(v)
        if len(adjs) <=1: return 0.0
        ligs = 0
        for i in adjs:
            for j in adjs:
                if i != j:
                    if j in self.graph[i] or i in self.graph[j]: 
                        ligs = ligs + 1
        return float(ligs)/(len(adjs)*(len(adjs)-1))
        
    def all_clustering_coefs(self) -> dict:
        '''Method that collects all clustering coefficients of the nodes in the graph. Clustering coefficient represents the tendency that a node has to form a cluster with the adjacents

        Returns
        -------
        dict
            Dictionary of all nodes and respective clustering coefficients
        '''
        ccs = {}
        for k in self.graph.keys():
            ccs[k] = self.clustering_coef(k)
        return ccs
        
    def mean_clustering_coef(self) -> float:
        '''Method that calculates the mean value of the clustering coefficients

        Returns
        -------
        float
            Mean value of the clustering coefficients
        '''
        ccs = self.all_clustering_coefs()
        return sum(ccs.values()) / float(len(ccs))
            
    def mean_clustering_perdegree(self, deg_type: str = "inout") -> dict:
        '''Method that calculates the mean clustering value per degree of the graph

        Parameters
        ----------
        deg_type : str, optional
            Informs which type of degree to compute and return, by default "inout".
            "inout": Count of the in degree and out degree
            "in": Count of the in degree
            "out": Count of the out degree
        Returns
        -------
        dict
            Dictionary with each degree and respective clustering value
        '''
        degs = self.all_degrees(deg_type)
        ccs = self.all_clustering_coefs()
        degs_k = {}
        for k in degs.keys():
            if degs[k] in degs_k.keys(): degs_k[degs[k]].append(k)
            else: degs_k[degs[k]] = [k]
        ck = {}
        for k in degs_k.keys():
            tot = 0
            for v in degs_k[k]: tot += ccs[v]
            ck[k] = float(tot) / len(degs_k[k])
        return ck

## Caminhos Hamiltonianos

    def check_if_valid_path(self, p:list) -> bool:
        '''Method to verify if the given path is valid in the graph

        Parameters
        ----------
        p : list
            Path to take in the graph

        Returns
        -------
        bool
            Verification if the given path is valid
        '''
        if p[0] not in self.graph.keys(): return False
        for i in range(1,len(p)):
            if p[i] not in self.graph.keys() or p[i] not in self.graph[p[i-1]]:
                return False
        return True

    def check_if_hamiltonian_path (self, p:list) -> bool:
        '''Verify if the path provided is an hamiltonian path. Hamiltonian path is a path that passes through all the nodes exactly once

        Parameters
        ----------
        p : list
            Path given to verify if is an hamiltonian path

        Returns
        -------
        bool
            Verification if the path provided is an hamiltonian path
        '''
        if not self.check_if_valid_path(p): return False
        visited = {n: 0 for n in self.graph.keys()}
        pos = 1
        cur = p[0]
        visited[cur] += 1
        while pos < len(p):
            if p[pos] in self.graph[cur]:
                cur = p[pos]
                pos += 1
                visited[cur] += 1
            else: return False
        if all(v == 1 for v in visited.values()): return True
        else: return False

    def search_hamiltonian_path(self) -> list:
        '''Method that depends on the "search_hamiltonian_path_from_node" that verifies in every node if exists an hamiltonian path

        Returns
        -------
        list
            Hamiltonian path obtained
        '''
        for ke in self.graph.keys():
            p = self.search_hamiltonian_path_from_node(ke)
            if p != None: return p
        return None
    
    def _possible_steps(self, v: Union[str, int, float], path: list, visited: list) -> list:
        '''Auxiliary function of "search_hamiltonian_path_from_node" to obtain the possibilities of steps of path to take in the graph. Returns the nodes that are not already present in the path and were not visited by the node before 

        Parameters
        ----------
        v : Union[str, int, float]
            Current node to obtain the possibilities
        path : list
            Current path built
        visited : list
            Previously visited nodes

        Returns
        -------
        list
            Possibilities of path steps to take
        '''
        poss = []
        for p in self.graph[v]:
            if p in path or p in visited:
                pass
            else:
                poss.append(p)
        if poss == []: return None
        else: return poss
  
    def search_hamiltonian_path_from_node(self, start: str) -> list:
        '''Searches an hamiltonian path from a given node. Hamiltonian path is a path that passes through all the nodes exactly once

        Parameters
        ----------
        start : str
            Node to start the hamiltonian path

        Returns
        -------
        list
            Hamiltonian path obtained. Returns None if no hamiltonian path is found
        '''
        cur = start
        visited = {n: [] for n in self.graph.keys()}
        path = [cur]
        while len(path) < len(self.get_nodes()):
            poss = self._possible_steps(cur, path, visited[cur])
            if poss != None:
                for n in poss:
                    visited[cur].append(n)
                    cur = n
                    path.append(cur)
                    break
            elif len(path) > 1:
                rmv_node = path.pop()
                visited[cur] = []
                cur = path[-1]
            else: return None
        return path

## Ciclos Eulerianos

    def check_balanced_node(self, node: Union[str, int, float]) -> bool:
        '''Fuction that verifies if a given node is balanced. It is balanced if the degree of entry of the node is equal to its out degree 

        Parameters
        ----------
        node : Union[str, int, float]
            Node to verify balancing

        Returns
        -------
        bool
            If balanced return True, otherwise False
        '''
        return self.in_degree(node) == self.out_degree(node)
    
    def check_balanced_graph(self) -> bool:
        '''Fuction that verifies if the graph is balanced

        Returns
        -------
        bool
            If balanced return True, otherwise False
        '''
        for n in self.graph.keys():
            if not self.check_balanced_node(n): return False
        return True

    def eulerian_cycle(self) -> list:
        '''Method that obtains an eulerian cycle from the graph

        Returns
        -------
        list
            List of the nodes that form the path cycle
        '''
        # Execução: Obtém todos os edges de visita a partir daí tenta criar logo um ciclo pelo 1º, retirando sistematicamente os caminhos realizados dos "edges_visit"
        #           Caso ainda n esteja vazio o "edges_visit", volta a fazer um ciclo para com os edges não visitados ainda
        #           No fnal, junta o novo ciclo/ciclos ao realizado(s) anteriormente a partir do lugar do nodo inicial
        if not self.check_balanced_graph(): return None
        edges_visit = list(self.get_edges())
        res = []
        while edges_visit:
            if res:
                for i in range(len(edges_visit)):
                    if edges_visit[i][0] in res:
                        cycle = list(edges_visit.pop(i))
                        break
            else: cycle = list(edges_visit.pop(0))
            nxt = cycle[1]
            while nxt != cycle[0]:
                for new in self.graph[nxt]:
                    if (nxt, new) in edges_visit:
                        edges_visit.remove((nxt, new))
                        nxt = new
                        cycle.append(nxt)
            if not res: res = cycle
            else:
                pos = res.index(cycle[0])
                for i in range(len(cycle) - 1): res.insert(pos + i + 1, cycle[i + 1])
                        
        return res

    def check_nearly_balanced_graph(self) -> tuple:
        '''Method that verifies if the graph possesses 2 nodes semi-balanced. They are respectivly the start and finishing elements of an eulerian path

        Returns
        -------
        tuple
            (Start, Finish) are respectivly the positions of a possible eulerian path
        '''
        res = None, None
        for n in self.graph.keys():
            indeg = self.in_degree(n)
            outdeg = self.out_degree(n)
            if indeg - outdeg == 1 and res[1] is None: res = res[0], n      # + entrada do nodo que saída (para terminar o path)
            elif indeg - outdeg == -1 and res[0] is None: res = n, res[1]   # + saída do nodo que entrada (para começar o path)
            elif indeg == outdeg: pass
            else: return None, None
        return res

    def eulerian_path(self) -> list:
        '''Method that builds the eulerian path. Uses the "check_nearly_balanced_graph" function to obtain the first and last member of the path and the "eulerian_cycle" function to build the cycle essential to build it 

        Returns
        -------
        list
            List of nodes that form the eulerian path
        '''
        unb = self.check_nearly_balanced_graph()
        if unb[0] is None or unb[1] is None: return None
        self.graph[unb[1]][unb[0]] = None
        cycle = self.eulerian_cycle()
        for i in range(len(cycle) - 1):
            if cycle[i] == unb[1] and cycle[i + 1] == unb[0]: break
        path = cycle[i+1:] + cycle[1:i+1]                               # 2º "i+1" é porque no "cycle[i+1:]" já incluí o 1º e último membro do ciclo 
        return path

def is_in_tuple_list(tl: Union[list, tuple], val: Union[str, int, float]) -> bool:
    '''Method that verifies if an element "val" exists in the list/tuple (tl)

    Parameters
    ----------
    tl : Union[list, tuple]
        List/tuple to check if element "val" exists
    val : Union[str, int, float]
        Element to check if present in "tl"

    Returns
    -------
    bool
        If present, returns True. Otherwise, returns False
    '''
    res = False
    for (x,y) in tl:
        if val == x: return True
    return res


def test1():
    print("* Teste 1 *\n")
    gr = MyGraph({1:{2: None}, 2:{3: None}, 3:{2: None, 4: None}, 4:{2: None}})
    gr.print_graph()
    print (gr.get_nodes())
    print (gr.get_edges())
    

def test2():
    print("\n* Teste 2 *\n")
    gr2 = MyGraph()
    gr2.add_vertex(1)
    gr2.add_vertex(2)
    gr2.add_vertex(3)
    gr2.add_vertex(4)
    
    gr2.add_edge(1,2)
    gr2.add_edge(2,3)
    gr2.add_edge(3,2)
    gr2.add_edge(3,4)
    gr2.add_edge(4,2)
    
    gr2.print_graph()
    print (gr2.get_nodes())
    print (gr2.get_edges())
  
def test3():
    print("\n* Teste 3 *\n")
    gr = MyGraph({1:{2: None}, 2:{3: None}, 3:{2: None, 4: None}, 4:{2: None}})
    gr.print_graph()

    print (gr.get_successors(2))
    print (gr.get_predecessors(2))
    print (gr.get_adjacents(2))
    print (gr.in_degree(2))
    print (gr.out_degree(2))
    print (gr.degree(2))

def test4():
    print("\n* Teste 4 *\n")
    gr = MyGraph({1:{2: None}, 2:{3: None}, 3:{2: None, 4: None}, 4:{2: None}})
    print ("1º Graph:")
    print (gr.distance(1,4))
    print (gr.distance(4,3))
    
    print (gr.shortest_path(1,4))
    print (gr.shortest_path(4,3))

    print (gr.reachable_with_dist(1))
    print (gr.reachable_with_dist(3))
    
    gr2 = MyGraph({1:{2: None, 3: None}, 2:{4: None}, 3:{5: None}, 4:{}, 5:{}})
    print ("\n2º Graph:")
    print (gr2.distance(2,1))
    print (gr2.distance(1,5))
    
    print (gr2.shortest_path(1,5))
    print (gr2.shortest_path(2,1))

    print("\n* Reachable *\n")
    print (gr2.reachable_bfs(1))
    print (gr2.reachable_dfs(1))

    print (gr2.reachable_with_dist(1))
    print (gr2.reachable_with_dist(5))
    
def test5():
    print("\n* Teste 5 *\n")
    gr = MyGraph({1:{2: None}, 2:{3: None}, 3:{2: None, 4: None}, 4:{2: None}})
    print (gr.node_has_cycle(2))
    print (gr.node_has_cycle(1))
    print (gr.has_cycle())
    
    gr2 = MyGraph({1:{2: None, 3: None}, 2:{4: None}, 3:{5: None}, 4:{}, 5:{}})
    print (gr2.node_has_cycle(1))
    print (gr2.has_cycle())
    print("\n* Test Degrees *\n")
    print(gr2.all_degrees())
    print(gr2.mean_degree())
    print(gr2.prob_degree())
    print("\n* Test Distances *\n")
    print(gr2.mean_distances())    
    print("\n* Test Clustering *\n")
    print(gr.clustering_coef(3))
    print(gr.all_clustering_coefs())
    print(gr.mean_clustering_coef())
    print(gr.mean_clustering_perdegree())

if __name__ == "__main__":
    test1()
    test2()
    test3()
    test4()
    test5()
