#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.


import unittest
from Indiv import all
from Popul import all
from EvolAlgorithm import all
from EAMotifs import all


class testEAMotifs(unittest.TestCase):
    def setUp(self):
        '''Test that sets the variables ...
        '''
        self.indiv1 = Indiv(20, [0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1])
        self.indiv2 = Indiv(10, [0, 1, 1, 0, 1, 0, 0, 0, 0, 0])
        self.indiv3 = IndivReal(10, [0.9641751112376019, 0.6406499646742604, 0.8761019524930286, 0.76514749395904, 0.32785733853093124, 0.4643870498872552, 0.02180230783590298, 0.5668745944573185, 0.1070675037695521, 0.5250634027343329])
        self.indiv4 = Indiv(20, [0.6487115652713965, 0.24986920827276038, 0.8220391041437379, 0.5871829514358591, 0.8672132378484155, 0.7727937105790735, 0.4135297287345566, 0.42814310333507566, 0.343050484175356, 0.5409426647821827, 0.9485994307734075, 0.5993694607836115, 0.6776607423620291, 0.8992901984960642, 0.4212212120091351, 0.12104676252615232, 0.374693462705072, 0.6981942543221408, 0.8052405588946782, 0.8499835290774598])
        

        
    def testMutation(self):
        self.indiv1.mutation()
        self.assertFalse(self.indiv1.genes == [0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1])
        self.indiv2.mutation()
        self.assertFalse(self.indiv2.genes == [0, 1, 1, 0, 1, 0, 0, 0, 0, 0])
        self.indiv3.mutation()
        self.assertFalse(self.indiv3.genes == [0.9641751112376019, 0.6406499646742604, 0.8761019524930286, 0.76514749395904, 0.32785733853093124, 0.4643870498872552, 0.02180230783590298, 0.5668745944573185, 0.1070675037695521, 0.5250634027343329])
        self.indiv4.mutation()
        self.assertFalse(self.indiv4.genes == [0.6487115652713965, 0.24986920827276038, 0.8220391041437379, 0.5871829514358591, 0.8672132378484155, 0.7727937105790735, 0.4135297287345566, 0.42814310333507566, 0.343050484175356, 0.5409426647821827, 0.9485994307734075, 0.5993694607836115, 0.6776607423620291, 0.8992901984960642, 0.4212212120091351, 0.12104676252615232, 0.374693462705072, 0.6981942543221408, 0.8052405588946782, 0.8499835290774598])
        
    def testFitness(self):
        self.fit1 = self.indiv1.getFitness()
        print(self.fit1)
        self.fit2 = self.indiv1.setFitness(2)
        self.assertEqual(self.fit1, self.fit2)       

if __name__ == '__main__':
    unittest.main()