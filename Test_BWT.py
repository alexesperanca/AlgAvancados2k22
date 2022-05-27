# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from BWT import BWT

class TestBWT(unittest.TestCase):
    def setUp(self):
        self.t1=BWT("TAGACAGAGA$")
        self.t2=BWT('ACG$GTAAAAC')

    def test_build_BWT(self):
        self.assertEqual(self.t1.bwt_dic,
        {'A0': '$0', 'G0': 'A0', 'G1': 'A1', 'G2': 'A2', 'T0': 'A3', 'C0': 'A4', 'A1': 'C0', 'A2': 'G0', 'A3': 'G1', 'A4': 'G2', '$0': 'T0'})
        self.assertEqual(self.t1.bwt_line,
        ['A0', 'G0', 'G1', 'G2', 'T0', 'C0', 'A1', 'A2', 'A3', 'A4', '$0'])
        self.assertEqual(self.t1.ord_line,
        ['$0', 'A0', 'A1', 'A2', 'A3', 'A4', 'C0', 'G0', 'G1', 'G2', 'T0'])
        self.assertEqual(self.t1.combinations,
        [('$TAGACAGAGA', 10), ('A$TAGACAGAG', 9), ('ACAGAGA$TAG', 3), ('AGA$TAGACAG', 7), ('AGACAGAGA$T', 1), ('AGAGA$TAGAC', 5), ('CAGAGA$TAGA', 4), ('GA$TAGACAGA', 8), ('GACAGAGA$TA', 2), ('GAGA$TAGACA', 6), ('TAGACAGAGA$', 0)])
        self.assertEqual(self.t2.bwt_dic,
        {'G0': '$0', 'T0': 'A0', 'A0': 'A1', 'A1': 'A2', 'A2': 'A3', 'C0': 'A4', 'A3': 'C0', 'A4': 'C1', 'C1': 'G0', '$0': 'G1', 'G1': 'T0'})
        self.assertEqual(self.t2.bwt_line,
        ['G0', 'T0', 'A0', 'A1', 'A2', 'C0', 'A3', 'A4', 'C1', '$0', 'G1'])
        self.assertEqual(self.t2.ord_line,
        ['$0', 'A0', 'A1', 'A2', 'A3', 'A4', 'C0', 'C1', 'G0', 'G1', 'T0'])
        self.assertEqual(self.t2.combinations,
        [('$GTAAAACACG', 3), ('AAAACACG$GT', 6), ('AAACACG$GTA', 7), ('AACACG$GTAA', 8), ('ACACG$GTAAA', 9), ('ACG$GTAAAAC', 0), ('CACG$GTAAAA', 10), ('CG$GTAAAACA', 1), ('G$GTAAAACAC', 2), 
('GTAAAACACG$', 4), ('TAAAACACG$G', 5)])

    def test_original_seq(self):
        self.assertEqual(self.t1.original_seq(),'TAGACAGAGA')
        self.assertEqual(self.t2.original_seq(),'GTAAAACACG')

    def test_find_pattern(self):
        self.assertEqual(self.t1.find_pattern("AGA"),[3, 4, 5])
        self.assertEqual(self.t2.find_pattern('GAG'),'No match!')

    def test_suffixarray(self):
        self.assertEqual(self.t1.suffixarray(),(10, 9, 3, 7, 1, 5, 4, 8, 2, 6, 0))
        self.assertEqual(self.t2.suffixarray(),(3, 6, 7, 8, 9, 0, 10, 1, 2, 4, 5))

if __name__ == '__main__':
    unittest.main()