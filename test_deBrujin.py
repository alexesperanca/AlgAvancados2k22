# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from debruijn import DeBruijnGraph
from MyGraph import MyGraph

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        self.frags = ["ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
        self.frags1 = ['AATG', 'ATGC', 'ATGG', 'CAAT', 'GCAA', 'GGTC', 'GTCT', 'TCTG', 'TGCA', 'TGGT']

        self.dbgr = DeBruijnGraph(self.frags)
        self.dbgr1 = DeBruijnGraph(self.frags1)
        self.dbgr2 = MyGraph({1:[2], 2:[3,1], 3:[4], 4:[2,5], 5:[6], 6:[4]})

    def test_check_nearly_balanced_graph(self):
        self.assertEqual(self.dbgr.check_nearly_balanced_graph(), (None, None))
        self.assertEqual(self.dbgr1.check_nearly_balanced_graph(), ('ATG', 'CTG'))

    def test_eulerian_path(self):
        self.assertEqual(self.dbgr.eulerian_path(), None)
        self.assertEqual(self.dbgr1.eulerian_path(), ['ATG', 'TGC', 'GCA', 'CAA', 'AAT', 'ATG', 'TGG', 'GGT', 'GTC', 'TCT', 'CTG'])

    def test_seq_from_path(self):
        self.assertEqual(self.dbgr1.seq_from_path(['ATG', 'TGC', 'GCA', 'CAA', 'AAT', 'ATG', 'TGG', 'GGT', 'GTC', 'TCT', 'CTG']), "ATGCAATGGTCTG")

    def test_check_balanced_graph(self):
        self.assertTrue(self.dbgr2.check_balanced_graph())

    def test_eulerian_cycle(self):
        self.assertEqual(self.dbgr2.eulerian_cycle(), [1, 2, 3, 4, 5, 6, 4, 2, 1])

if __name__ == '__main__':
    unittest.main()