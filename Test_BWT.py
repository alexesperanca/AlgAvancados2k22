import unittest
from BWT import BWT

class TestBWT(unittest.TestCase):
    def setUp(self):
        self.t1=BWT("TAGACAGAGA$")

    def test_build_BWT(self):
        self.assertEqual(self.t1.bwt_dic,
        {'A0': '$0', 'G0': 'A0', 'G1': 'A1', 'G2': 'A2', 'T0': 'A3', 'C0': 'A4', 'A1': 'C0', 'A2': 'G0', 'A3': 'G1', 'A4': 'G2', '$0': 'T0'})
        self.assertEqual(self.t1.bwt_line,
        ['A0', 'G0', 'G1', 'G2', 'T0', 'C0', 'A1', 'A2', 'A3', 'A4', '$0'])
        self.assertEqual(self.t1.ord_line,
        ['$0', 'A0', 'A1', 'A2', 'A3', 'A4', 'C0', 'G0', 'G1', 'G2', 'T0'])
        self.assertEqual(self.t1.combinations,
        [('$TAGACAGAGA', 10), ('A$TAGACAGAG', 9), ('ACAGAGA$TAG', 3), ('AGA$TAGACAG', 7), ('AGACAGAGA$T', 1), ('AGAGA$TAGAC', 5), ('CAGAGA$TAGA', 4), ('GA$TAGACAGA', 8), ('GACAGAGA$TA', 2), ('GAGA$TAGACA', 6), ('TAGACAGAGA$', 0)])

    def test_original_seq(self):
        self.assertEqual(self.t1.original_seq(),'TAGACAGAGA')

    def test_find_pattern(self):
        self.assertEqual(self.t1.find_pattern("AGA"),[3, 4, 5])

    def test_suffixarray(self):
        self.assertEqual(self.t1.suffixarray(),(10, 9, 3, 7, 1, 5, 4, 8, 2, 6, 0))

if __name__ == '__main__':
    unittest.main()