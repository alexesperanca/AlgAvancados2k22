# -*- coding: utf-8 -*-
"""
@author: Alexandre Esperança
"""

# from MySeq import MySeq
# from MyMotifs import MyMotifs

from Sequence import Sequence
from Motifs import Motifs

class MotifFinding:
    
    def __init__(self, size = 8, seqs = None):
        self.motifSize = size
        if (seqs != None):
            self.seqs = seqs
            self.alphabet = 'DNA'
        else:
            self.seqs = []

    def __len__ (self):
        return len(self.seqs)
    
    def __getitem__(self, n):
        return self.seqs[n]
    
    def seqSize (self, i):
        return len(self.seqs[i])
    
    def readFile(self, fic):
        for s in open(fic, 'r'):
            self.seqs.append(s.strip().upper())
        # self.alphabet = self.seqs[0].alphabet()
        
    def createMotifFromIndexes(self, indexes):
        pseqs = []
        for i,ind in enumerate(indexes):
            sequencia = self.seqs[i][ind:ind+self.motifSize]
            # print(sequencia)
            pseqs.append(sequencia)
            # pseqs.append(Sequence(self.seqs[i].seq[ind:(ind+self.motifSize)]))
        return Motifs(pseqs)
        
        
    # SCORES
        
    def score(self, s):
        score = 0
        motif = self.createMotifFromIndexes(s)
        mat = motif.__create_prof__()
        for dic in mat:
            m = max(*dic.values())
            score += m
        return score
    
    
    def scoreMult(self, s):
        score = 1.0
        motif = self.createMotifFromIndexes(s)
        mat = motif.__create_prof__()
        for dic in mat:
            m = max(*dic.values())
            score *= m
        return score
 
       
    # EXHAUSTIVE SEARCH_type_
       
    def nextSol (self, s):
        nextS = [0]*len(s)
        pos = len(s) - 1     
        while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if (pos < 0): 
            nextS = None
        else:
            for i in range(pos): 
                nextS[i] = s[i]
            nextS[pos] = s[pos]+1;
            for i in range(pos+1, len(s)):
                nextS[i] = 0
        return nextS
        
    def exhaustiveSearch(self, seqs):
        melhorScore = -1
        res = []
        s = [0]* len(seqs)
        while (s!= None):
            sc = self.score(s)
            if (sc > melhorScore):
                melhorScore = sc
                res = s
            s = self.nextSol(s)
        return res
     
    # BRANCH AND BOUND     
     
    def nextVertex (self, s):
        res =  []
        if len(s) < len(self.seqs): # internal node -> down one level
            for i in range(len(s)): 
                res.append(s[i])
            res.append(0)
        else: # bypass
            pos = len(s)-1 
            while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
                pos -= 1
            if pos < 0: res = None # last solution
            else:
                for i in range(pos): res.append(s[i])
                res.append(s[pos]+1)
        return res
    
    
    def bypass (self, s):
        res =  []
        pos = len(s) -1
        while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if pos < 0: res = None 
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos]+1)
        return res
        
    def branchAndBound (self):
        melhorScore = -1
        melhorMotif = None
        size = len(self.seqs)
        s = [0]*size
        while s != None:
            if len(s) < size:
                optimScore = self.score(s) + (size-len(s)) * self.motifSize
                if optimScore < melhorScore: s = self.bypass(s)
                else: s = self.nextVertex(s)
            else:
                sc = self.score(s)
                if sc > melhorScore:
                    melhorScore = sc
                    melhorMotif  = s
                s = self.nextVertex(s)
        return melhorMotif

    # Consensus (heuristic)
  
    def heuristicConsensus(self) -> list:
        '''Algoritmo heuristico (tipo consensus) com método exaustivo
        
        Returns
        -------
        list
            list of index for each sequence
        '''
        motif_finding = MotifFinding(self.motifSize, self.seqs) ## primeiras duas sequencias fixas
        search = motif_finding.exhaustiveSearch(self.seqs[:2]) # melhores posições para primeiras duas seq
        for i, seq in enumerate(self.seqs[2:]): # para as restantes sequencias // como tem enumerate, é necessário depois somar 2 ao 'i':
            search.append(0) # adiciona-se um zero como valor inicial para cada seq adicional (index = 0, é o primeiro a ser testado)
            best_score = -1
            best_ind = 0 # index onde ocorre melhor score
            for ind in range((len(seq))-self.motifSize+1): # testar cada index (com limite do tamanho da sequencia pelo motif)
                search[i+2] = ind 
                s_current = self.score(search)
                if s_current > best_score:
                    best_score = s_current #melhor score que o atual na lista
                    best_ind = ind 
            search[i+2] = best_ind
        return search

    # Consensus (Stochastic)

    def heuristicStochastic (self):
        from random import randint
        # ...
        return None

    # Gibbs sampling 

    def gibbs (self):
        from random import randint
        # ...
        return None

    def roulette(self, f):
        from random import random
        tot = 0.0
        for x in f: tot += (0.01+x)
        val = random()* tot
        acum = 0.0
        ind = 0
        while acum < val:
            acum += (f[ind] + 0.01)
            ind += 1
        return ind-1

# tests

def test1(): 
    print ("Test 1:")
    sm = MotifFinding()
    sm.readFile("exemploMotifs.txt")
    sol = [25,20,2,55,59]
    sa = sm.score(sol)
    print(sa)
    scm = sm.scoreMult(sol)
    print(scm)

def test2():
    print ("\nTest 2 - exhaustive:")
    seq1 = "ATAGAGCTGA"
    seq2 = "ACGTAGATGA"
    seq3 = "AAGATAGGGG"
    seq1 = "ATAGAGCTGA"
    seq2 = "ACGTAGATGA"
    seq3 = "AAGATAGGGG"
    mf = MotifFinding(3, [seq1,seq2,seq3])
    sol = mf.exhaustiveSearch([seq1,seq2,seq3])
    print ("Solution", sol)
    print ("Score: ", mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

    print ("\nTest 2 - Branch and Bound:")
    sol2 = mf.branchAndBound()
    print ("Solution: " , sol2)
    print ("Score:" , mf.score(sol2))
    print("Consensus:", mf.createMotifFromIndexes(sol2).consensus())
    
    print ("\nTest 2 - Heuristic consensus: ")
    sol1 = mf.heuristicConsensus()
    print ("Solution: " , sol1)
    print ("Score:" , mf.score(sol1))

def test3():
    print ("\nTest 3:")
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt")
    print ("Branch and Bound:")
    sol = mf.branchAndBound()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

def test4():
    print ("\nTest 4:")
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt")
    print("Heuristic stochastic")
    sol = mf.heuristicStochastic()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print ("Score mult:" , mf.scoreMult(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())
    
    sol2 = mf.gibbs(1000)
    print ("Score:" , mf.score(sol2))
    print ("Score mult:" , mf.scoreMult(sol2))


if __name__ == '__main__':
    
    test1()
    test2()
    test3()
    test4()

