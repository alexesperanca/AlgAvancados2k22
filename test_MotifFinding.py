import unittest
from MotifFinding import MotifFinding

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        '''Test that it can ...
        '''
        list_seq = ["ATAGAGCTGA", "ACGTAGATGA", "AAGATAGGGG"]
        mf = MotifFinding(3, list_seq)
        self.sol = mf.exhaustiveSearch()

    def test_assert(self):
        self.assertEqual(self.sol, [1, 3, 4], 'should be')



if __name__ == '__main__':
    unittest.main()

