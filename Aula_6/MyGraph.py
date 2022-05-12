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

    def _path(self, s, d):
        '''
        Código testa um caminho em frente nas possibilidades, caso chegue a um ponto sem fim, retorna a 2 anteriores (porque a anterior dará o mesmo resultado)
        A partir daí volta a testar, caso não haja anteriores (como no teste 4) verificamos os antecessores e corremos caminhos "em frente" por aí
        Caso não funcione, dará o mesmo erro e portanto devemos testar novamente os antecessores e assim sucessivamente...
        
        Disclaimer: Acho que poderá haver casos onde poderia haver um erro a ter o caminho caso seja muito grande o problema ;_;
                    Mas aqui passa todos os testes

        Nova forma: Fazer 1 caminho para a frente e outro para trás, verificar qual mais pequeno e fica esse -> Grafo e o seu contrário fazer caminho normal (com retornos caso sem fim)
        '''
        visited = [s]
        path = {s: 0}
        current = s
        poss = self.graph[current]
        while current != d:
            value = path[current]
            for n in poss:
                previous = current
                if n in visited: pass
                elif self.graph[n] == {} and n != d: 
                    visited.append(n)
                    pass
                else:
                    current = n
                    path[current] = value + 1
                    visited.append(n)
                    break
            
            if previous == current:
                if current != s:
                    del path[current]
                    current = list(path.keys())[-1]
                    poss = self.graph[current]
                else:
                    try:
                        current = list(path.keys())[-2] # Caso esteja no início dos caminhos para dar erro e entrar nos antecessores
                        poss = self.graph[current]
                    except:
                        poss = self.get_predecessors(current) # Ir através dos antecessores
            else:
                poss = self.graph[current]
        '''Mekieeee'''
        return path, current

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
        path, last = self._path(s, d)
        return path[last]
        
    def shortest_path(self, s, d):
        if s == d: return [s,d]
        path, l = self._path(s, d)
        nodes = "".join(f"{i} -> " for i in path)
        return nodes[:-4]
        
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
    

if __name__ == "__main__":
    test1()
    test2()
    test3()
    test4()
    test5()
