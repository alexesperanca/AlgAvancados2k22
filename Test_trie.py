import unittest
from trie import SuffixTree, Trie
import pprint

class Testtrie(unittest.TestCase):
    def setUp(self):
        self.t1=Trie(["CTG", "CATA", "CAAGG"])
        self.t2=Trie(["AGAGAT", "AGC", "AGTCC", "CAGAT", "CCTA", "GAGAT", "GAT", "TC"])
        self.t1i=self.t1.insert('GGGA')
    
    def test_insert(self):
        self.assertNotEqual(self.t1i,self.t1) 
# é uma forma de testar, pq se não está igual entao, é porque algo mudou, 
# mas não se consegue verificar se inseriu corretamente

    def test_print_tree(self):
        self.assertEqual(self.t1.print_tree(),pprint.pprint(self.t1, width = 1))
        self.assertEqual(self.t2.print_tree(),pprint.pprint(self.t2, width = 1))
        self.assertEqual(self.t1.print_tree(),pprint.pprint(self.t1i, width = 1))
#não sei se isto conta

    # def test_trie_matches(self):
        # self.assertEqual(self.t1.trie_matches('CTGCATACAAGG'),)
        # self.assertEqual(self.t2.trie_matches('GAGATCCTA'),)
    
    def test_match(self):
        pass


class TestSuffixTree(unittest.TestCase):
    def setUp(self):
        self.t3=SuffixTree("TACTA")
        self.t4=SuffixTree('AATTTCGATCGATTGAT')

    def test_print_tree(self):
        self.assertEqual(self.t3.print_tree(),pprint.pprint(self.t3, width = 1))

    def test_find_pattern(self):
        self.assertEqual(self.t3.find_pattern('TATA'),[('TA', 0), ('A', 1), ('TA', 2), ('A', 3)])
        self.assertEqual(self.t3.find_pattern('ACG'),'No match!')
        self.assertEqual(self.t3.find_pattern('ACGT'),'No match!')

#     def test_add_suffix(self):
#         self.assertEqual(self.t3.add_suffix('GGAT'),)

#     def test_get_leafes_bellow(self):
#         self.assertEqual(self.t3.get_leafes_below('C'),)

#     def test_find_pattern(self):
#         self.assertEqual(self.t3.find_pattern('ACTA'),
#         [('ACTA', 0), ('CTA', 1), ('TA', 2), ('A', 3)])

#     def test_repeats(self):
#         self.assertEqual(self.t3.repeats('TA'),2)


if __name__ == '__main__':
    unittest.main()