# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from BoyerMoore import BoyerMoore

class TestBoyerMoore (unittest.TestCase):
    def setUp(self):
        self.t1 = BoyerMoore('ACTG','ACCA')#.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC")
        self.t2 = BoyerMoore('ACTG','CTTA')
        self.t3 = BoyerMoore('ACTG','ACWA')
        self.t4 = BoyerMoore('actg','acca')
    
    def test_search_pattern(self):
        self.assertEqual(self.t1.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"),[5, 13, 23, 37])
        self.assertEqual(self.t2.search_pattern('CGTGCCTACTTACTTACTTACTTACGCGAA'),[8, 12, 16, 20])
        self.assertEqual(self.t3.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"),('No match!'))
        self.assertEqual(self.t4.search_pattern("atagaaccaatgaaccatgatgaaccatggatacccaaccacc"),[5, 13, 23, 37])
        self.assertEqual(self.t1.search_pattern("ACTGACTGACTGACTGACTGGTGTAGCAGGAGCGAGCAGGTATTATATGC"),('No match!'))

if __name__ == '__main__':
    unittest.main()
