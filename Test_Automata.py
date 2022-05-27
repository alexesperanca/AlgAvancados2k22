import unittest
from Aula_1.Automata import overlap
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
        self.assertEqual(self.t1.transitionTable,{(0, 'A'): 1, (0, 'C'): 0, (1, 'A'): 1, (1, 'C'): 2, (2, 'A'): 3, (2, 'C'): 0, (3, 'A'): 1, (3, 'C'): 2})
        
        self.assertEqual(self.t2.numstates,5)
        self.assertEqual(self.t2.alphabet,('ACTG'))
        self.assertEqual(self.t2.transitionTable,{(0, 'A'): 1, (0, 'C'): 0, (0, 'T'): 0, (0, 'G'): 0, (1, 'A'): 1, (1, 'C'): 2, (1, 'T'): 0, (1, 'G'): 0, (2, 'A'): 1, (2, 'C'): 3, (2, 'T'): 0, (2, 'G'): 0, (3, 'A'): 4, (3, 'C'): 0, (3, 'T'): 0, (3, 'G'): 0, (4, 'A'): 1, (4, 'C'): 2, (4, 'T'): 0, (4, 'G'): 0})
        
        self.assertEqual(self.t3.numstates,5)
        self.assertEqual(self.t3.alphabet,('ACTG'))
        self.assertEqual(self.t3.transitionTable,{(0, 'A'): 0, (0, 'C'): 1, (0, 'T'): 0, (0, 'G'): 0, (1, 'A'): 0, (1, 'C'): 1, (1, 'T'): 2, (1, 'G'): 0, (2, 'A'): 0, (2, 'C'): 1, (2, 'T'): 3, (2, 'G'): 0, (3, 'A'): 4, (3, 'C'): 1, (3, 'T'): 0, (3, 'G'): 0, (4, 'A'): 0, (4, 'C'): 1, (4, 'T'): 0, (4, 'G'): 0})
        
        self.assertEqual(self.t4.numstates,5)
        self.assertEqual(self.t4.alphabet,('ACTG'))
        self.assertEqual(self.t4.transitionTable,{(0, 'A'): 0, (0, 'C'): 0, (0, 'T'): 1, (0, 'G'): 0, (1, 'A'): 0, (1, 'C'): 2, (1, 'T'): 1, (1, 'G'): 0, (2, 'A'): 0, (2, 'C'): 0, (2, 'T'): 1, (2, 'G'): 3, (3, 'A'): 4, (3, 'C'): 0, (3, 'T'): 1, (3, 'G'): 0, (4, 'A'): 0, (4, 'C'): 0, (4, 'T'): 1, (4, 'G'): 0})

    def test_applySeq(self):
        self.assertEqual(self.t1.applySeq('CACAACAA'),[0, 0, 1, 2, 3, 1, 2, 3, 1])
        self.assertEqual(self.t4.applySeq('CACAACAA'),[0, 0, 0, 0, 0, 0, 0, 0, 0])
    
    def test_occurencesPattern(self):
        self.assertEqual(self.t1.occurencesPattern('CACAACAA'),[1,4])

    def test_nextState(self):
        self.assertEqual(self.t1.nextState(2,'A'),3)

    # def test_overlap(self):
    #     self.assertEqual(self.overlap('abb','ababab'),0)

if __name__ == '__main__':
    unittest.main()