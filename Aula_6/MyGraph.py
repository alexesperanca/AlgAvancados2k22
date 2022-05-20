# -*- coding: utf-8 -*-

## Graph represented as adjacency list using a dictionary
## keys are vertices
## values of the dictionary represent the list of adjacent vertices of the key node

class MyGraph:
    
    def __init__(self, g = {}):
        ''' Constructor - takes dictionary to fill the graph as input; default is empty dictionary '''
        self.graph = g    

    def print_graph(self):
        ''' Prints the content of the graph as adjacency list '''
        for v in self.graph.keys():
            print (v, " -> ", self.graph[v])

    ## get basic info

    def get_nodes(self):
        ''' Returns list of nodes in the graph '''
        return list(self.graph.keys())
        
    def get_edges(self): 
        ''' Returns edges in the graph as a list of tuples (origin, destination) '''
        edges = []
        for v in self.graph.keys():
            for d in self.graph[v]:
                edges.append((v,d))
        return edges
      
    def size(self):
        ''' Returns size of the graph : number of nodes, number of edges '''
        return len(self.get_nodes()), len(self.get_edges())
      
    ## add nodes and edges    
    
    def add_vertex(self, v):
        ''' Add a vertex to the graph; tests if vertex exists not adding if it does '''
        if v not in self.graph:
            self.graph[v] = {}
        
    def add_edge(self, o, d):
        ''' Add edge to the graph; if vertices do not exist, they are added to the graph ''' 
        if o not in self.graph:
            self.graph[o] = {d: None}
        else:
            g = self.graph[o]
            g[d] = None

    ## successors, predecessors, adjacent nodes
        
    def get_successors(self, v):
        return list(self.graph[v].keys())     # needed to avoid list being overwritten of result of the function is used
             
    def get_predecessors(self, v):
        return [k for k in self.graph.keys() if v in self.graph[k]]
    
    def get_adjacents(self, v):
        pred = self.get_predecessors(v)
        suc = self.get_successors(v)
        for i in suc:
            if i not in pred: pred.append(i)
        return pred
        
    ## degrees    
    
    def out_degree(self, v):
        return len(self.graph[v])
    
    def in_degree(self, v):
        return len(self.get_predecessors(v))
        
    def degree(self, v):
        return len(self.get_adjacents(v))
        
    def all_degrees(self, deg_type = "inout"):
        ''' Computes the degree (of a given type) for all nodes.
        deg_type can be "in", "out", or "inout" '''
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
    
    def reachable_bfs(self, v):
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res
        
    def reachable_dfs(self, v):
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

    def _path(self, graph, s, d):
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

    def get_path(self, s, d):
        '''Return the best path'''
        frontpath = self._path(self.graph, s, d)
        rev_graph = self.revert_graph()
        backpath = self._path(rev_graph, s, d)
        if frontpath == None and backpath == None: return None 
        elif backpath == None and frontpath != None: return frontpath
        elif frontpath == None and backpath != None: return backpath
        elif len(frontpath) <= len(backpath): return frontpath
        else: return backpath

    def revert_graph(self):
        '''Reverse of the original graph'''

        def add_edge(graph, o, d):
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

    def distance(self, s, d):
        if s == d: return 0
        path = self.get_path(s, d)
        return len(path) - 1
        
    def shortest_path(self, s, d):
        if s == d: return [s,d]
        path = self.get_path(s, d)
        nodes = " -> ".join(f"{i}" for i in path)
        return nodes
        
    def reachable_with_dist(self, s):
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
    def mean_distances(self):
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
    def node_has_cycle (self, v):
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

    def has_cycle(self):
        res = False
        for v in self.graph.keys():
            if self.node_has_cycle(v): return True
        return res

## topological metrics over degrees

    def mean_degree(self, deg_type = "inout"):
        degs = self.all_degrees(deg_type)
        return sum(degs.values()) / float(len(degs))
        
    def prob_degree(self, deg_type = "inout"):
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
        
    def clustering_coef(self, v):
        adjs = self.get_adjacents(v)
        if len(adjs) <=1: return 0.0
        ligs = 0
        for i in adjs:
            for j in adjs:
                if i != j:
                    if j in self.graph[i] or i in self.graph[j]: 
                        ligs = ligs + 1
        return float(ligs)/(len(adjs)*(len(adjs)-1))
        
    def all_clustering_coefs(self):
        ccs = {}
        for k in self.graph.keys():
            ccs[k] = self.clustering_coef(k)
        return ccs
        
    def mean_clustering_coef(self):
        ccs = self.all_clustering_coefs()
        return sum(ccs.values()) / float(len(ccs))
            
    def mean_clustering_perdegree(self, deg_type = "inout"):
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

    def check_if_valid_path(self, p ):
        if p[0] not in self.graph.keys(): return False
        for i in range(1,len(p)):
            if p[i] not in self.graph.keys() or p[i] not in self.graph[p[i-1]]:
                return False
        return True

    def check_if_hamiltonian_path (self, p):
        if not self.check_if_valid_path(p): return False
        #... 

def is_in_tuple_list (tl, val):
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
