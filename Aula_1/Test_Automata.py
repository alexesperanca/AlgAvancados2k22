import unittest
from Automata import Automata

class TestAutomata (unittest.TestCase):
    def setUp(self):
        self.t1=Automata('AC','ACA')
        self.t2=Automata('ACTG','ACCA')
        self.t3=Automata('ACTG','CTTA')
        self.t4=Automata('ACTG','TCGA')

    
    def test_printAutomata(self):
        self.assertEqual(self.t1.numstates,4)
        self.assertEqual(self.t1.alphabet,('AC'))
        self.assertEqual(self.t1.transitionTable,{(0, 'A'): 1,
                                                    (0, 'C'): 0,
                                                    (1, 'A'): 1,
                                                    (1, 'C'): 2,
                                                    (2, 'A'): 3,
                                                    (2, 'C'): 0,
                                                    (3, 'A'): 1,
                                                    (3, 'C'): 2})
        self.assertEqual(self.t2.numstates,5)
        self.assertEqual(self.t2.alphabet,('ACTG'))
        self.assertEqual(self.t2.transitionTable,{(0, 'A'): 1,
                                                    (0, 'C'): 0,
                                                    (1, 'A'): 1,
                                                    (1, 'C'): 2,
                                                    (2, 'A'): 3,
                                                    (2, 'C'): 0,
                                                    (3, 'A'): 1,
                                                    (3, 'C'): 4,
                                                    (4, 'A'): 3,
                                                    (4, 'C'): 2})
        self.assertEqual(self.t3.numstates,5)
        self.assertEqual(self.t3.alphabet,('ACTG'))
        self.assertEqual(self.t3.transitionTable,{(0, 'C'): 1,
                                                    (0, 'T'): 0,
                                                    (1, 'C'): 1,
                                                    (1, 'T'): 2,
                                                    (2, 'C'): 3,
                                                    (2, 'T'): 0,
                                                    (3, 'C'): 1,
                                                    (3, 'T'): 4,
                                                    (4, 'C'): 3,
                                                    (4, 'T'): 2})
        self.assertEqual(self.t4.numstates,5)
        self.assertEqual(self.t4.alphabet,('ACTG'))
        self.assertEqual(self.t4.transitionTable,{(0, 'C'): 1,
                                                    (0, 'T'): 0,
                                                    (1, 'C'): 1,
                                                    (1, 'T'): 2,
                                                    (2, 'C'): 3,
                                                    (2, 'T'): 0,
                                                    (3, 'C'): 1,
                                                    (3, 'T'): 4,
                                                    (4, 'C'): 3,
                                                    (4, 'T'): 2})

    def test_applySeq(self):
        self.assertEqual(self.t1.applySeq('CACAACAA'),[0, 0, 1, 2, 3, 1, 2, 3, 1])
    
    def test_occurencesPattern(self):
        self.assertEqual(self.t1.occurencesPattern('CACAACAA'),(1,4))

if __name__ == '__main__':
    unittest.main()