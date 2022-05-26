#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.


import unittest
from MotifFinding import MotifFinding

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        '''Test that it can ...
        '''
        #Create Class with sequences or files:
        self.mf1 = MotifFinding(3, ["ATAGAGCTGA", "ACGTAGATGA", "AAGATAGGGG"])
        self.mf2 = MotifFinding(3, ["AGGTACTT", "CCATACGT", "ACGTTAGT", "ACGTCCAT", "CCGTACGG"])
           
        self.mf3 = MotifFinding()
        self.mf3.readFile("exemploMotifs.txt")
        self.mf4 = MotifFinding()
        self.mf4.readFile("exemploMotifs2.txt")
        self.mf5 = MotifFinding()
        self.mf5.readFile("exemploMotifs3.txt")
        
        #Creating solutions:
        ## From Exhaustive Search:
        self.sol1 = self.mf1.exhaustiveSearch()

        ## Mannually:
        self.sol2 = [25,20,2,55,59]
        self.sol3 = [2,1,5,2]


        ## From Branch and Bound:
        self.sol4 = self.mf3.branchAndBound()
        # self.sol3 = self.mf3.readFile("exemploMotifs2.txt")
        # self.sol4 = self.mf4.readFile("exemploMotifs3.txt")

        ## From Consensus (heuristic):

        ## From Consensus (stochastic):

        ## From Gibbs Sampling:


        # Creating scores for each solution:
        self.sc1 = self.mf1.score(self.sol1)

        # Creating profiles from indexes:


        # Creating Consensus:





    def testSolution(self):
        self.assertEqual(self.sol1, [1, 3, 4], 'should be')
        self.assertEqual(self.sol2, [25, 20, 2, 55, 59], 'should be')

    def testScore(self):
        self.assertEqual(self.sc1, 9, 'should be')


if __name__ == '__main__':
    unittest.main()

