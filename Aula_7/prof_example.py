# -*- coding: utf-8 -*-

import io
import pprint
import re
import numpy as np

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
            print(f"{u} -> {v};", file = F)
            for v in g[u]:
                if g[u][v] is not None:
                    print(f'{u} -> {v}[label = "{g[u][v]}"];', file = F)
                else:
                    print(f"{u} -> {v};", file = F)
        print("}", file = F)
        return F.getvalue()
        
g = create_graph()
#add_node(g, 7)
#add_edge(g, 1, 2, 7)
#add_edge(g, 2, 3)
#add_edge(g, 2, 4)

#pprint.pprint(g)

'''
Aula seguinte (7) de Reações!!
Aqui vamos adicionar ao grafo g as reações, metabolitos e respetivas ligações, sendo cada reação a conversão de reagentes a produtos
'''

def separar_metabolitos(exp):
    return [M.strip() for M in re.split(r'\s*\+\s*', exp)]

def create_met_react_graph(ficheiro):
    g = create_graph()
    for reacao in ficheiro:
        R, Eq = reacao.split(':')
        bidirec = '<=>' in Eq
        if bidirec:
            esq, dir = [separar_metabolitos(X) for X in Eq.split('<=>')]
        else:
            esq, dir = [separar_metabolitos(X) for X in Eq.split('=>')]
        #print(R, bidirec, esq, dir)

        for M in esq:
            add_edge(g, M, R)
            if bidirec: add_edge(g, R, M)

        for M in dir:
            add_edge(g, R, M)
            if bidirec: add_edge(g, M, R)
    return g

def is_met(M): return M[:2] == "M_"
def is_react(M): return M[:2] == "R_"

def converte_graph(criterio, g):
    '''
    Associa seja metabolitos ou reações, aos que é possível chegar a partir de cada
    '''
    novo = create_graph()
    for M1 in g:
        if criterio(M1):
            for R in g[M1]:
                for M2 in g[R]:
                    add_edge(novo, M1, M2)
    return novo

def outdeg(g, node):
    '''
    Conta os graus associados a cada nodo -> Ou seja, length do dicionário associado
    '''
    return len(g[node])

def inverse(g):
    '''
    Permite verificar quais reações dão origem aos metabolitos e, por outro lado, quais metabolitos dão origem às reações
    '''
    novo = create_graph()
    for U in g:
        for V in g[U]:
            add_edge(novo, V, U)
    return novo

def graus(g):
    invg = inverse(g)
    return [outdeg(g, M) + outdeg(invg, M) for M in g]

def alcancaveis(g, M):
    '''
    A partir de um metabolito, a quais somos capazes de alcançar através deste
    '''
    lst = [(M, 0)]
    alc = {}
    while lst:
        U, dist = lst[0]
        lst = lst[1:]
        alc[U] = dist
        for V in g[U]:
            if V not in alc:
                lst.append((V, dist + 1))
    return alc

def distancias(g):
    return {U: alcancaveis(g, U) for U in g}

def dist_med(dists):
    return np.mean([dists[U][V] for U in dists for V in dists[U]])

ficheiro = io.StringIO("""R_1: M_1 + M_2 => M_3 + M_4
R_2: M_4 + M_6 => M_3
R_3: M_4 + M_5 <=> M_6""")

g = create_met_react_graph(ficheiro)
gmet = converte_graph(is_met, g)
greac = converte_graph(is_react, g)

'''
with open("ecoli.txt") as ficheiro:
    g = create_met_react_graph(ficheiro)
'''

pprint.pprint(g)
print()
#pprint.pprint(gmet)
#print()
#pprint.pprint(greac) # Dá o próprio metabolito associado porque a reação é reversível
#print()
#pprint.pprint(inverse(g))
#print(outdeg(g, "M_4"))
pprint.pprint(distancias(g))
print(dist_med(distancias(g)))