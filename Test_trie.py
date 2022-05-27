# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from trie import SuffixTree, Trie
import pprint

class Testtrie(unittest.TestCase):
    def setUp(self):
        self.t1=Trie(["CTG", "CATA", "CAAGG"])
        self.t2=Trie(["AGAGAT", "AGC", "AGTCC", "CAGAT", "CCTA", "GAGAT", "GAT", "TC"])
    
    def test_insert(self):
        self.t1.insert('GGGA')
        self.assertEqual(self.t1.seq_list,['CTG', 'CATA', 'CAAGG', 'GGGA'])

    def test_print_tree(self):
        self.assertEqual(self.t1.print_tree(),pprint.pprint(self.t1, width = 1))
        self.assertEqual(self.t2.print_tree(),pprint.pprint(self.t2, width = 1))
    
    def test_match(self):
        self.assertEqual(self.t1.match('CTGCATACAAGG'),[('CTG', 0), ('CATA', 3), ('CAAGG', 7)])
        self.assertEqual(self.t2.match('GAGATCCTA'),[('GAGAT', 0), ('GAT', 2), ('TC', 4), ('CCTA', 5)])
        self.assertEqual(self.t1.match('GAGATCCTA'),'No match!')

class TestSuffixTree(unittest.TestCase):
    def setUp(self):
        self.t3=SuffixTree("TACTA")
        self.t4=SuffixTree('ACGT')

    def test_print_tree(self):
        self.assertEqual(self.t3.print_tree(),pprint.pprint(self.t3, width = 1))

    def test_find_pattern(self):
        self.assertEqual(self.t3.find_pattern_in_seq('TATA'),[('TA', 0), ('TA', 2)])
        self.assertEqual(self.t3.find_pattern_in_seq('ACG'),'No match!')
        self.assertEqual(self.t4.find_pattern_in_seq('AATTTCGACGTCGATTGAT'),[('ACGT', 7), ('CGT', 8), ('GT', 9)])
        
    def test_insert(self):
        self.t3.add_suffix('GGAT')
        self.assertEqual(self.t3.seq_list,['TACTA', 'ACTA', 'CTA', 'TA', 'A', 'GGAT'])


    def test_get_leafes_bellow(self):
        self.assertEqual(self.t3.get_leafes_below('C'),['CTA'])
        self.assertEqual(self.t3.get_leafes_below('A'),['ACTA', 'A#$#'])
        self.assertRaises(AssertionError,self.t3.get_leafes_below,'Node inputted not present in the tree')

    def test_repeats(self):
        self.assertEqual(self.t3.repeats('TA'),1)
        self.assertEqual(self.t3.repeats('CT'),9)


if __name__ == '__main__':
    unittest.main()