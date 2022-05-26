#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.


import unittest
from MotifFinding import MotifFinding

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        '''Test that sets the variables ...
        '''
        #Create object with sequences or files:
        self.mf1 = MotifFinding(3, ["ATAGAGCTGA", "ACGTAGATGA", "AAGATAGGGG"]) #test2
                #test6
        self.mf2 = MotifFinding(4, ["AGGTACTT", "CCATACGT", "ACGTTAGT", "ACGTCCAT", "CCGTACGG"])
        self.mf3 = MotifFinding(seqs = ["GTAAACAATATTTATAGC", "AAAATTTACCTCGCAAGG", "CCGTACTGTCAAGCGTGG", "TGAGTAAACGACGTCCCA", "TACTTAACACCCTGTCAA"])
                #test1 / 3 / 4
        self.mf4 = MotifFinding()
        self.mf4.readFile("exemploMotifs.txt") 
        self.mf5 = MotifFinding()
        self.mf5.readFile("exemploMotifs2.txt")
        self.mf6 = MotifFinding()
        self.mf6.readFile("exemploMotifs3.txt")
        
        #Creating solutions:
        ## From Exhaustive Search (small sequences):
        self.sol1 = self.mf1.exhaustiveSearch()
        self.sol2 = self.mf2.exhaustiveSearch()
        self.sol3 = self.mf3.exhaustiveSearch()
        
        ## From Branch and Bound:
        self.sol4 = self.mf1.branchAndBound()
        self.sol5 = self.mf2.branchAndBound()
        self.sol6 = self.mf3.branchAndBound()
        # self.sol7 = self.mf4.branchAndBound()
        self.sol8 = self.mf5.branchAndBound()

        ## From Consensus (heuristic):
        self.sol9  = self.mf1.heuristicConsensus()
        self.sol10 = self.mf2.heuristicConsensus()
        self.sol11 = self.mf3.heuristicConsensus()
        self.sol12 = self.mf4.heuristicConsensus()
        self.sol13 = self.mf5.heuristicConsensus()
        self.sol14 = self.mf6.heuristicConsensus()

        ## From Consensus (stochastic):
        # self.sol15 = self.mf1.heuristicStochastic()
        self.sol16 = self.mf2.heuristicStochastic()
        # self.sol17 = self.mf3.heuristicStochastic()
        self.sol18 = self.mf4.heuristicStochastic()
        # self.sol19 = self.mf5.heuristicStochastic()
        # self.sol20 = self.mf6.heuristicStochastic()


        ## From Gibbs Sampling:
        self.sol21 = self.mf1.gibbs(1000)
        self.sol22 = self.mf2.gibbs(1000)
        self.sol23 = self.mf3.gibbs(1000)
        self.sol24 = self.mf4.gibbs(500)
        self.sol25 = self.mf5.gibbs(500)
        self.sol26 = self.mf6.gibbs(500)

        # Creating scores for each solution:
        self.sc1 = self.mf1.score(self.sol1)
        self.sc2 = self.mf2.score(self.sol2)
        self.sc3 = self.mf3.score(self.sol3)

        self.sc4 = self.mf1.score(self.sol4)
        self.sc5 = self.mf4.score(self.sol5)
        self.sc6 = self.mf4.score(self.sol6)
        # self.sc7 = self.mf4.score(self.sol7)
        self.sc8 = self.mf4.score(self.sol8)

        self.sc9  = self.mf1.score(self.sol9)
        self.sc10 = self.mf2.score(self.sol10)
        self.sc11 = self.mf3.score(self.sol11)
        self.sc12 = self.mf4.score(self.sol12)
        self.sc13 = self.mf5.score(self.sol13)
        self.sc14 = self.mf6.score(self.sol14)

        # self.sc15 = self.mf1.score(self.sol15)
        self.sc16 = self.mf2.score(self.sol16)
        # self.sc17 = self.mf3.score(self.sol17)
        self.sc18 = self.mf4.score(self.sol18)
        # self.sc19 = self.mf5.score(self.sol19)
        # self.sc20 = self.mf6.score(self.sol20)

        self.sc21 = self.mf1.score(self.sol21)
        self.sc22 = self.mf2.score(self.sol22)
        self.sc23 = self.mf3.score(self.sol23)
        self.sc24 = self.mf4.score(self.sol24)
        self.sc25 = self.mf5.score(self.sol25)
        self.sc26 = self.mf6.score(self.sol26)

        # Creating Consensus:
        self.cons1 = self.mf1.createMotifFromIndexes(self.sol1).consensus()
        self.cons2 = self.mf2.createMotifFromIndexes(self.sol2).consensus()
        self.cons3 = self.mf3.createMotifFromIndexes(self.sol3).consensus()


    def testSolution(self):
        ## exhaustive solution
        self.assertEqual(self.sol1, [1, 3, 4])
        self.assertEqual(self.sol2, [0, 4, 0, 0, 0])
        self.assertEqual(self.sol3, [0, 8, 4, 0, 10])
        
        ## branchAndBound
        self.assertEqual(self.sol4, [1, 3, 4])
        self.assertEqual(self.sol5, [0, 4, 0, 0, 0])
        self.assertEqual(self.sol6, [0, 8, 4, 0, 10])
        # self.assertEqual(self.sol7, [1, 4, 45, 5, 0])
        self.assertEqual(self.sol8, [25, 20, 2, 55, 59])
                
        ## consensus (heuristic) 
        self.assertEqual(self.sol9, [1, 3, 4])
        self.assertEqual(self.sol10, [0, 4, 0, 0, 0])
        self.assertEqual(self.sol11, [6, 0, 4, 2, 1])
        self.assertEqual(self.sol12, [0, 38, 14, 33, 1])
        self.assertEqual(self.sol13, [25, 20, 2, 55, 59])
        self.assertEqual(self.sol14, [0, 0, 235, 229, 72, 7, 90, 210, 239, 76])        
        
        ## consensus (stochastic) 
        # self.assertEqual(self.sol16, [1, 1, 2, 2, 2])
        # self.assertEqual(self.sol18, [27, 22, 0, 53, 57])  

        ## gibbs sampling
        # self.assertEqual(self.sol21, [1, 3, 4])
        # self.assertEqual(self.sol22, [0, 4, 0, 0, 4])
        # self.assertEqual(self.sol23, [0, 4, 0, 0, 4])

    def testScore(self):
        self.assertEqual(self.sc1, 9, 'should be')
        # self.assertEqual(self.sc5, 20, 'should be')


    def testConsensus(self):
        self.assertEqual(self.cons1, 'TAG', 'should be')
        # self.assertEqual(self.cons2, 'TAG', 'should be')


if __name__ == '__main__':
    unittest.main()

