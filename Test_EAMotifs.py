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
        random.seed(8)
        self.best_sol_int1, self.bestFit_int1 = self.EAmotif_int1.run()
        random.seed(8)
        self.best_sol_int2, self.bestFit_int2 = self.EAmotif_int2.run()
        random.seed(8)
        self.best_sol_int3, self.bestFit_int3 = self.EAmotif_int3.run()
        random.seed(8)
        self.cons_int1 = self.EAmotif_int1.motifs.createMotifFromIndexes(self.best_sol_int1).consensus()
        random.seed(8)
        self.cons_int2 = self.EAmotif_int2.motifs.createMotifFromIndexes(self.best_sol_int2).consensus()
        random.seed(8)
        self.cons_int3 = self.EAmotif_int3.motifs.createMotifFromIndexes(self.best_sol_int3).consensus()


        random.seed(8)
        self.EAmotif_real1 = EAMotifsReal(100, 1000, 50, 'exemploMotifs.txt')
        random.seed(8)
        self.EAmotif_real2 = EAMotifsReal(100, 1000, 50, 'exemploMotifs2.txt')

        random.seed(8)
        self.best_sol_real1, self.bestFit_real1 = self.EAmotif_real1.run()
        random.seed(8)
        self.best_sol_real2, self.bestFit_real2 = self.EAmotif_real2.run()
        random.seed(8)
        self.o, self.prof1 = self.EAmotif_real1.profile(self.best_sol_real1)
        random.seed(8)
        self.cons_real1 = self.EAmotif_real1.consensus(self.prof1)
        random.seed(8)
        self.o2, self.prof2 = self.EAmotif_real2.profile(self.best_sol_real2)
        random.seed(8)
        self.cons_real2 = self.EAmotif_real2.consensus(self.prof2)


    def testBestSolution(self):
        self.assertEqual(self.best_sol_int1, [0, 38, 24, 26, 1])
        self.assertEqual(self.best_sol_int2, [12, 19, 1, 7, 48])
        self.assertEqual(self.best_sol_int3, [1, 1, 1, 0, 1, 0, 69, 250, 31, 0])

        self.assertEqual(self.best_sol_real1, [58.52366296729691, 59.37254807355931, 6.333301053603342, 5.16372443049355, 12.572643370540177, 15.268068412756417, 1.802455768194473, 44.75382087215941, 30.494768165093532, 24.061487290146168, 39.12221782456808, 34.25608854316616, 34.14744884803779, 8.259515452458846, 18.683887371769845, 17.962615593712783, 30.295397708522948, 9.862257285755174, 39.14717319408232, 55.26851095002273])
        self.assertEqual(self.best_sol_real2, [26.40782019112007, 5.80989991514441, 14.148002147504322, 17.623161342898033, 17.900565581041967, 54.13807565887457, 1.9813272282447758, 40.53815707286813, 35.03135750051371, 16.119149697930496, 43.5192251835127, 17.582636848494882, 4.663511792805075, 33.79947458550591, 15.995763423489143, 26.671818192289606, 59.56041480812477, 8.459140357753194, 14.197536134704634, 35.62867435711404])

    def testBestFitness(self):
        self.assertEqual(self.bestFit_int1, 28)
        self.assertEqual(self.bestFit_int2, 27)
        self.assertEqual(self.bestFit_int3, 49)

        self.assertEqual(self.bestFit_real1, 34)
        self.assertEqual(self.bestFit_real2, 40)

    def testConsensus(self):
        self.assertEqual(self.cons_int1, 'CCTGATAC')
        self.assertEqual(self.cons_int2, 'TACGTACG')
        self.assertEqual(self.cons_int3, 'AATCATTA')

        self.assertEqual(self.cons_real1, 'CTGAG')
        self.assertEqual(self.cons_real2, 'ACGCA')


if __name__ == '__main__':
    unittest.main()