# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from overlap_graphs import OverlapGraph

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        self.frags1 = ['ACC-2', 'CCA-8', 'CAT-5', 'ATG-3']
        self.frags2 = ['ACC-2', 'CCA-8', 'CAT-5', 'ATG-3', 'TGG-13', 'GGC-10', 'GCA-9', 'CAT-6', 'ATT-4', 'TTT-15', 'TTC-14', 'TCA-12', 'CAT-7', 'ATA-1', 'TAA-11']
        self.frags3 = ['CAA-7', 'AAT-1', 'ATC-2', 'TCA-12', 'CAT-8', 'ATG-4', 'TGA-13', 'GAT-9', 'ATG-5', 'TGA-14', 'GAT-10', 'ATG-6', 'TGA-15', 'GAT-11', 'ATC-3']
        self.frags4 = ['ACC', 'CCA', 'CAT', 'ATG']
        self.frags5 = ['ACC', 'CCA', 'CAT', 'ATG', 'TGG', 'GGC', 'GCA', 'CAT', 'ATT', 'TTT', 'TTC', 'TCA', 'CAT', 'ATA', 'TAA']
        self.frags6 = ['ATT', 'TTA', 'GCA', 'CAA']
        self.frags7 = ['acC', 'CCA', 'cat', 'ATG']

        self.ovgr = OverlapGraph(["ATA",  "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA" , "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"], True)
        self.ovgr1 = OverlapGraph(['AAT', 'ATC', 'ATC', 'ATG', 'ATG', 'ATG', 'CAA', 'CAT', 'GAT', 'GAT', 'GAT', 'TCA', 'TGA', 'TGA', 'TGA'], True)
        self.ovgr2 = OverlapGraph(["ATA",  "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA" , "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"], False)

    def test_check_if_valid_path(self):
        self.assertTrue(self.ovgr.check_if_valid_path(self.frags1))
        self.assertTrue(self.ovgr.check_if_valid_path(self.frags2))
        self.assertTrue(self.ovgr2.check_if_valid_path(self.frags4))
        self.assertTrue(self.ovgr2.check_if_valid_path(self.frags5))
        self.assertFalse(self.ovgr2.check_if_valid_path(self.frags6))
        self.assertFalse(self.ovgr2.check_if_valid_path(self.frags7))

    def test_check_if_hamiltonian_path(self):
        self.assertFalse(self.ovgr.check_if_hamiltonian_path(self.frags1))
        self.assertTrue(self.ovgr.check_if_hamiltonian_path(self.frags2))
        self.assertTrue(self.ovgr1.check_if_hamiltonian_path(self.frags3))
        self.assertFalse(self.ovgr1.check_if_hamiltonian_path(self.frags4))
        self.assertFalse(self.ovgr1.check_if_hamiltonian_path(self.frags5))

    def test_seq_from_path(self):
        self.assertEqual(self.ovgr.seq_from_path(self.frags2), "ACCATGGCATTTCATAA")
        self.assertEqual(self.ovgr1.seq_from_path(self.frags3), "CAATCATGATGATGATC")
        self.assertEqual(self.ovgr2.seq_from_path(self.frags5), None)

    def test_search_hamiltonian_path(self):
        self.assertEqual(self.ovgr.search_hamiltonian_path(), self.frags2)
        self.assertEqual(self.ovgr1.search_hamiltonian_path(), self.frags3)
        self.assertEqual(self.ovgr2.search_hamiltonian_path(), None)

if __name__ == '__main__':
    unittest.main()
