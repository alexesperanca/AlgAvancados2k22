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


        random.seed(8)
        self.EAmotif_real1 = EAMotifsReal(100, 1000, 50, 'exemploMotifs.txt')
        random.seed(8)
        self.EAmotif_real2 = EAMotifsReal(100, 1000, 50, 'exemploMotifs2.txt')

        self.best_sol_real1, self.bestFit_real1 = self.EAmotif_real1.run()
        self.best_sol_real2, self.bestFit_real2 = self.EAmotif_real2.run()
        self.o, self.prof1 = self.EAmotif_real1.profile(self.best_sol_real1)
        self.cons_real1 = self.EAmotif_real1.consensus(self.prof1)
        self.o2, self.prof2 = self.EAmotif_real2.profile(self.best_sol_real2)
        self.cons_real2 = self.EAmotif_real2.consensus(self.prof2)


    def testBestSolution(self):
        self.assertEqual(self.best_sol_int1, [48, 47, 1, 0, 44])
        self.assertEqual(self.best_sol_int2, [28, 43, 1, 0, 58])
        self.assertEqual(self.best_sol_int3, [0, 0, 60, 63, 281, 151, 0, 0, 178, 0])

        self.assertEqual(self.best_sol_real1, [48, 47, 1, 0, 44])
        self.assertEqual(self.best_sol_real2, [28, 43, 1, 0, 58])
        # self.assertEqual(self.best_sol_real3, [0, 0, 60, 63, 281, 151, 0, 0, 178, 0])

    def testBestFitness(self):
        self.assertEqual(self.bestFit_int1, 27)
        self.assertEqual(self.bestFit_int2, 33)
        self.assertEqual(self.bestFit_int3, 52)

        self.assertEqual(self.bestFit_real1, 27)
        self.assertEqual(self.bestFit_real2, 33)
        # self.assertEqual(self.bestFit_real3, 52)

    def testConsensus(self):
        self.assertEqual(self.cons_int1, 'AACGCTCG')
        self.assertEqual(self.cons_int2, 'AACGTACG')
        self.assertEqual(self.cons_int3, 'AAAATCTT')

        self.assertEqual(self.cons_real1, 'ACGTA')
        self.assertEqual(self.cons_real2, 'ACGTA')


if __name__ == '__main__':
    unittest.main()