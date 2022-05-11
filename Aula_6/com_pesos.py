# -*- coding: utf-8 -*-

import io
import pprint

# Para desenhar o grafo: Graphviz
# com pesos

def create_graph():
    return {}

def add_node(graph, node):
    if node not in graph:
        graph[node] = {}
    
def add_edge(graph, u, v, w = None):
    add_node(graph, u)
    add_node(graph, v)
    graph[u][v] = w
    
def to_graphviz(g):
    with io.StringIO() as F:
        print("""digraph{
            node[shape = "circle", style = "filled"];
            """, file = F)
        for u in g:
            for v in g[u]:
                if g[u][v] is not None:
                    print(f'{u} -> {v}[label = "{g[u][v]}"];', file = F)
                else:
                    print(f"{u} -> {v};", file = F)
        print("}", file = F)
        return F.getvalue()
        
g = create_graph()
add_node(g, 7)
add_edge(g, 1, 2, 7)
add_edge(g, 2, 3)
add_edge(g, 2, 4)

pprint.pprint(g)
pprint.pprint(to_graphviz(g))