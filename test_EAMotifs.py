#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.


import unittest
from Indiv import Indiv, IndivInt, IndivReal
from Popul import Popul, PopulInt, PopulReal
from EvolAlgorithm import EvolAlgorithm
from EAMotifs import EAMotifsInt, EAMotifsReal
import random


class testEAMotifs(unittest.TestCase):
    def setUp(self):
        random.seed(8)
        self.EAmotif_int1 = EAMotifsInt(100, 1000, 50, 'exemploMotifs.txt')
        random.seed(8)
        self.EAmotif_int2 = EAMotifsInt(200, 2000, 60, 'exemploMotifs2.txt')
        random.seed(8)
        self.EAmotif_int3 = EAMotifsInt(50, 1000, 30, 'exemploMotifs3.txt')
        self.best_sol_int1, self.bestFit_int1 = self.EAmotif_int1.run()
        self.best_sol_int2, self.bestFit_int2 = self.EAmotif_int2.run()
        self.best_sol_int3, self.bestFit_int3 = self.EAmotif_int3.run()
        self.cons_int1 = self.EAmotif_int1.motifs.createMotifFromIndexes(self.best_sol_int1).consensus()
        self.cons_int2 = self.EAmotif_int2.motifs.createMotifFromIndexes(self.best_sol_int2).consensus()
        self.cons_int3 = self.EAmotif_int3.motifs.createMotifFromIndexes(self.best_sol_int3).consensus()


        # random.seed(8)
        self.EAmotif_real1 = EAMotifsReal(100, 1000, 50, 'exemploMotifs.txt')
        self.EAmotif_real2 = EAMotifsReal(200, 2000, 60, 'exemploMotifs2.txt')
        self.EAmotif_real3 = EAMotifsReal(50, 1000, 30, 'exemploMotifs3.txt')


    def testBestSolution(self):
        self.assertEqual(self.best_sol_int1, [0, 38, 24, 26, 1])
        self.assertEqual(self.best_sol_int2, [33, 24, 6, 59, 1])

    def testBestFitness(self):
        self.assertEqual(self.bestFit_int1, 28)
        self.assertEqual(self.bestFit_int2, 31)

    def testConsensus(self):
        self.assertEqual(self.cons_int1, 'CCTGATAC')
        self.assertEqual(self.cons_int2, 'ACGTACAC')


if __name__ == '__main__':
    unittest.main()