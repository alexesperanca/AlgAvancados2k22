# -*- coding: utf-8 -*-
"""
@authors:
- Alexandre Esperança
- Cristiana Martins
- Mónica Leiras
- Tomás Sá
"""

from Motifs import Motifs
from random import randint


class MotifFinding:
    
    def __init__(self, size = 8, seqs = None):
        self.motifSize = size
        if (seqs != None):
            self.seqs = seqs
            self.alphabet = self.alphabet()
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
        self.alphabet = self.alphabet()
    
    def alphabet(self) -> str:
        if all (i in 'ACGT' for i in self.seqs[0]) is True:
            return 'ACGT'
        elif all (i in 'ACGU' for i in self.seqs[0]) is True:
            return 'ACGU'
        elif all (i in 'FLIMVSPTAY_HQNKDECWRG' for i in self.seqs[0]) is True:
            return 'FLIMVSPTAY_HQNKDECWRG'
    
    def createMotifFromIndexes(self, indexes, pseudocontagem = 0):
        pseqs = []
        for i,ind in enumerate(indexes):
            sequencia = self.seqs[i][ind:ind+self.motifSize]
            pseqs.append(sequencia)
        return Motifs(pseqs, pseudo=pseudocontagem)
        
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
        
    def exhaustiveSearch(self):
        melhorScore = -1
        res = []
        s = [0]* len(self.seqs)
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
        motif_finding = MotifFinding(self.motifSize, self.seqs[:2]) ## primeiras duas sequencias fixas
        search = motif_finding.exhaustiveSearch() # melhores posições para primeiras duas seq
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
        search = [randint(0, self.seqSize(i)-self.motifSize) for i in range(len(self.seqs))]
        best_score = self.score(search)
        best = False
        while best is False:
            motif = self.createMotifFromIndexes(search)
            for i in range(len(self.seqs)):
                seq_prob, ind = motif.seq_most_probable(self.seqs[i])
                search[i] = ind
            score = self.score(search)
            if score > best_score:
                best_score = score
                best = True
        return search
        
        
    # Gibbs sampling 

    def gibbs (self, iterations):
        search = [randint(0, self.seqSize(i)-self.motifSize) for i in range(len(self.seqs))] # lista de index inicial
        best_search = list(search)
        best_score = self.score(best_search)
        for ite in range(iterations):
            ind_s1 = randint(0, len(self.seqs)-1)
            s1 = self.seqs[ind_s1]
            self.seqs.remove(s1)
            search.pop(ind_s1)
            motif = self.createMotifFromIndexes(search, pseudocontagem = 1)
            list_prob = []
            p = 0
            while len(s1[p:p+self.motifSize]) == len(motif.profile):
                prob = motif.prob_seq(s1[p:p+self.motifSize])
                list_prob.append(prob)
                p +=1
            position = self.roulette(list_prob)
            search.insert(ind_s1, position)
            self.seqs.insert(ind_s1, s1)
            score = self.score(search)
            if score > best_score:
                best_score = score
                best_search = list(search)
        return best_search

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
    sm.readFile("AlgAvancados2k22/exemploMotifs.txt")
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
    sol = mf.exhaustiveSearch()
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
    mf.readFile("AlgAvancados2k22/exemploMotifs.txt")
    print ("Branch and Bound:")
    sol = mf.branchAndBound()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

def test4():
    print ("\nTest 4:")
    mf = MotifFinding()
    mf.readFile("AlgAvancados2k22/exemploMotifs.txt")
    print("Heuristic stochastic")
    sol = mf.heuristicStochastic()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print ("Score mult:" , mf.scoreMult(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())
    
    sol2 = mf.gibbs(1000)
    print ("Solution 2: " , sol2)
    print ("Score:" , mf.score(sol2))
    print ("Score mult:" , mf.scoreMult(sol2))


if __name__ == '__main__':
    test1()
    test2()
    test3()
    test4()
