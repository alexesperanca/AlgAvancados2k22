#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.


import unittest
import random
from MotifFinding import MotifFinding

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        '''Test that sets the variables ...
        '''
        #Create object with sequences or files:
        self.mf1 = MotifFinding(3, ["ATAGAGCTGA", "ACGTAGATGA", "AAGATAGGGG"]) 
        self.mf2 = MotifFinding(4, ["AGGTACTT", "CCATACGT", "ACGTTAGT", "ACGTCCAT", "CCGTACGG"])
        self.mf3 = MotifFinding(seqs = ["GTAAACAATATTTATAGC", "AAAATTTACCTCGCAAGG", "CCGTACTGTCAAGCGTGG", "TGAGTAAACGACGTCCCA", "TACTTAACACCCTGTCAA"])
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
        # self.sol7 = self.mf4.branchAndBound() #muito tempo
        self.sol8 = self.mf5.branchAndBound()

        ## From Consensus (heuristic):
        random.seed(7)
        self.sol9  = self.mf1.heuristicConsensus()
        self.sol10 = self.mf2.heuristicConsensus()
        self.sol11 = self.mf3.heuristicConsensus()
        self.sol12 = self.mf4.heuristicConsensus()
        self.sol13 = self.mf5.heuristicConsensus()
        self.sol14 = self.mf6.heuristicConsensus()

        ## From Consensus (stochastic):
        random.seed(7)
        # self.sol15 = self.mf1.heuristicStochastic() #muito tempo
        self.sol16 = self.mf2.heuristicStochastic()
        # self.sol17 = self.mf3.heuristicStochastic() #muito tempo
        self.sol18 = self.mf4.heuristicStochastic()
        # self.sol19 = self.mf5.heuristicStochastic() #muito tempo
        # self.sol20 = self.mf6.heuristicStochastic() #muito tempo


        ## From Gibbs Sampling:
        random.seed(7)
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

        # self.sc15 = self.mf1.score(self.sol15) #muito tempo
        self.sc16 = self.mf2.score(self.sol16)
        # self.sc17 = self.mf3.score(self.sol17) #muito tempo
        self.sc18 = self.mf4.score(self.sol18)
        # self.sc19 = self.mf5.score(self.sol19) #muito tempo
        # self.sc20 = self.mf6.score(self.sol20) #muito tempo

        self.sc21 = self.mf1.score(self.sol21)
        self.sc22 = self.mf2.score(self.sol22)
        self.sc23 = self.mf3.score(self.sol23)
        self.sc24 = self.mf4.score(self.sol24)
        self.sc25 = self.mf5.score(self.sol25)
        self.sc26 = self.mf6.score(self.sol26)

        self.ms1 = self.mf1.scoreMult(self.sol21)
        self.ms2 = self.mf2.scoreMult(self.sol22)
        self.ms3 = self.mf3.scoreMult(self.sol23)
        self.ms4 = self.mf4.scoreMult(self.sol24)
        self.ms5 = self.mf5.scoreMult(self.sol25)
        self.ms6 = self.mf6.scoreMult(self.sol26)

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
        self.assertEqual(self.sol8, [25, 20, 2, 55, 59])
                
        ## consensus (heuristic) 
        self.assertEqual(self.sol9, [1, 3, 4])
        self.assertEqual(self.sol10, [0, 4, 0, 0, 0])
        self.assertEqual(self.sol11, [6, 0, 4, 2, 1])
        self.assertEqual(self.sol12, [0, 38, 14, 33, 1])
        self.assertEqual(self.sol13, [25, 20, 2, 55, 59])
        self.assertEqual(self.sol14, [0, 0, 235, 229, 72, 7, 90, 210, 239, 76])        
        
        ## consensus (stochastic) 
        random.seed(7)
        self.assertEqual(self.sol16, [2, 0, 0, 4, 0])
        self.assertEqual(self.sol18, [46, 34, 10, 23, 25])  

        ## gibbs sampling
        random.seed(7)
        self.assertEqual(self.sol21, [2, 4, 1])
        self.assertEqual(self.sol22, [0, 4, 0, 0, 0])
        self.assertEqual(self.sol23, [10, 4, 8, 3, 2])

    def testScore(self):
        self.assertEqual(self.sc1, 9)
        self.assertEqual(self.sc2, 18)
        self.assertEqual(self.sc3, 28)
        self.assertEqual(self.sc8, 30)
        self.assertEqual(self.sc11, 25)
        self.assertEqual(self.sc12, 30)
        self.assertEqual(self.sc14, 66)
        self.assertEqual(self.sc16, 14)
        self.assertEqual(self.sc18, 26)  
        self.assertEqual(self.sc21, 9)
        self.assertEqual(self.sc23, 27)

        self.assertEqual(self.ms1, 27.0)
        self.assertEqual(self.ms2, 400.0)
        self.assertEqual(self.ms3, 12960.0) #isto faz sentido?
        self.assertEqual(self.ms4, 20736.0)
        self.assertEqual(self.ms5, 17280.0)
        self.assertEqual(self.ms6, 784000.0)

    def testConsensus(self):
        self.assertEqual(self.cons1, 'TAG')
        self.assertEqual(self.cons2, 'ACGT')


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.
