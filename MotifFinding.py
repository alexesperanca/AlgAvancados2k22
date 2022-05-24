#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`MotifFinding` class, for finding conserved motifs, given a list of sequences.
This class includes diverse strategies for this, such as:
    - Exhaustive search type
    - Branch & Bound (search trees)
    - Consensus (heuristic algorithms)
    - Consensus (estochastic algorithms)
    - Gibbs sampling

"""

from Motifs import Motifs
from random import randint


class MotifFinding:
    
    def __init__(self, size: int = 8, seqs: list = None):
        '''Class that implements algorithms for finding conserved motifs, given a list of sequences,
        i.e. optimization algorithms for motif discovery/inference.

        Parameters
        ----------
        size : int, optional
            size of motif to find, by default 8
        seqs : _type_, optional
            list of sequences, by default None
        '''
        self.motifSize = size
        if (seqs != None):
            self.seqs = seqs
            self.alphabet = self.alphabet()
        else:
            self.seqs = []

    def __len__ (self) -> int:
        '''Method that returns length of the list of sequences.

        Returns
        -------
        int
            length of the list of sequences
        '''
        return len(self.seqs)
    
    def __getitem__(self, n: int) -> str:
        '''Method that returns a sequence with index n of the list of sequences

        Parameters
        ----------
        n : int
            number of the index

        Returns
        -------
        str
            sequence with index n
        '''
        return self.seqs[n]
    
    def seqSize (self, i: int) -> int:
        '''Method that returns the length of the sequence with index n of the list of sequences

        Parameters
        ----------
        i : int
            number of the index

        Returns
        -------
        int
            length of the sequence n
        '''
        return len(self.seqs[i])
    
    def readFile(self, fic: str):
        '''Method that reads a txt file with a list of sequences and appends them on the self.seqs

        Parameters
        ----------
        fic : str
            name of the file, example: 'file.txt'
        '''
        for s in open(fic, 'r'):
            self.seqs.append(s.strip().upper())
        self.alphabet = self.alphabet()
    
    def alphabet(self) -> str:
        '''Method that determines the type of alphabet of the sequences

        Returns
        -------
        str
            alphabet of all the sequences 
        '''
        if all (i in 'ACGT' for i in self.seqs[0]) is True: #DNA
            return 'ACGT'
        elif all (i in 'ACGU' for i in self.seqs[0]) is True: #RNA
            return 'ACGU'
        elif all (i in 'FLIMVSPTAY_HQNKDECWRG' for i in self.seqs[0]) is True: #AMINO
            return 'FLIMVSPTAY_HQNKDECWRG'
    
    def createMotifFromIndexes(self, indexes: list, pseudocount: int | float = 0) -> Motifs:
        '''Creates an instance of the Motifs Class - which calculates the probabilistic PWM (Position Weighted Matrix) profile of a given list of sequences.
        In this case, it is given also the indexes of the sequences to create this profile.

        Parameters
        ----------
        indexes : list
            list of indexes of the sequences
        pseudocount : int or float, optional
            pseudocount of the profile, by default 0

        Returns
        -------
        list
            Instance of the Class Motifs (a list of dictionaries)
        '''
        pseqs = []
        for i,ind in enumerate(indexes):
            sequence = self.seqs[i][ind:ind+self.motifSize]
            pseqs.append(sequence)
        return Motifs(pseqs, pseudo = pseudocount)

    # SCORES
        
    def score(self, s: list) -> float | int:
        '''Method that determines the additive score, given a list of indexes (calculating a profile)

        Parameters
        ----------
        s : list
            list of indexes

        Returns
        -------
        float or int
            score
        '''
        score = 0
        motif = self.createMotifFromIndexes(s)
        mat = motif.__create_prof__()
        for dic in mat:
            m = max(*dic.values())
            score += m
        return score
    
    
    def scoreMult(self, s: list) -> int | float:
        '''Method that determines the multiplicative score, given a list of indexes (calculating a profile)

        Parameters
        ----------
        s : list
            list of indexes

        Returns
        -------
        int or float
            score
        '''
        score = 1.0
        motif = self.createMotifFromIndexes(s)
        mat = motif.__create_prof__()
        for dic in mat:
            m = max(*dic.values())
            score *= m
        return score
 

    # EXHAUSTIVE SEARCH_type_
       
    def nextSol (self, s: list) -> list:
        '''Auxiliary method (to exhaustiveSearch) that gives the next vector of starting positions.

        Parameters
        ----------
        s : list
            current vector of starting positions s

        Returns
        -------
        list
            the next vector of starting positions s
        '''
        nextS = [0]*len(s)
        pos = len(s) - 1     
        while pos >= 0 and s[pos] == self.seqSize(pos) - self.motifSize:
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
        
    def exhaustiveSearch(self) -> list:
        '''Method that calculates the score of each possible vector of starting positions s.
        The best score determines the respective profile and consensus pattern (Method score).
        The objective is to maximize Score(s,seqs) by varying the starting positions.

        Returns
        -------
        list
            vector of starting positions s
        '''
        melhorScore = -1
        res = []
        s = [0]* len(self.seqs)
        while (s!= None):
            sc = self.score(s)
            if (sc > melhorScore):
                melhorScore = sc
                res = s
            s = self.nextSol(s) # auxiliary method 
        return res
     
    # BRANCH AND BOUND     
     
    def nextVertex (self, s: list) -> list:
        '''_summary_

        Parameters
        ----------
        s : list
            current vector of starting positions s (current Vertex)

        Returns
        -------
        list
            the next vector of starting positions s (next Vertex)
        '''
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
    
    
    def bypass (self, s: list) -> list:
        '''Method that skips ahead of the branches of a given vertex (bypass).

        Parameters
        ----------
        s : list
            current vector of starting positions s (current vertex)

        Returns
        -------
        list
            the next vector of starting positions s (next Vertex)
        '''
        res =  []
        pos = len(s) -1
        while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if pos < 0: res = None 
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos]+1)
        return res
        
    def branchAndBound (self) -> list:
        '''Method that implements Branch and Bound algorithm with four types of movements: next leaf; visit all leaves (similar to exhaustive search); visit next node, and bypass.

        Returns
        -------
        list
            the vector of starting positions s (representative of the best motif found)
        '''
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
        '''Method that computes the heuristic consensus algorithm: 
        - Considering only the first two sequences, choose the initial positions s1 and s2 that give a better score.
        - For each of the following sequences, iteratively, choose the best starting position in the sequence, in order to maximize the score.

        Returns
        -------
        list
            the vector of starting positions s (representative of the best motif found)
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

    def heuristicStochastic (self) -> list:
        '''Method that computes the heuristic stochastic consensus algorithm, using the most likely segments to adjust starting positions to achieve the best profile (motif).

        Returns
        -------
        list
            the vector of starting positions s (representative of the best motif found)
        '''
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

    def gibbs (self, iterations: int) -> list:
        '''Method that implements the Gibbs Sampling algorithm, by choosing new segments at random (increasing the possibilities of converging to a correct solution).

        Parameters
        ----------
        iterations : int
            maximum number of iterations

        Returns
        -------
        list
            the vector of starting positions s (representative of the best motif found)
        '''
        search = [randint(0, self.seqSize(i)-self.motifSize) for i in range(len(self.seqs))] # lista de index inicial
        best_search = list(search)
        best_score = self.score(best_search)
        for ite in range(iterations):
            ind_s1 = randint(0, len(self.seqs)-1)
            s1 = self.seqs[ind_s1]
            self.seqs.remove(s1)
            search.pop(ind_s1)
            motif = self.createMotifFromIndexes(search, pseudocount = 1)
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

    def roulette(self, f: list) -> int:
        '''Method auxiliary that determines the chosen position by its probability (probability of choosing a certain position is proportional to its score).

        Parameters
        ----------
        f : list
            list of probabilities 

        Returns
        -------
        int
            chosen position 
        '''
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
