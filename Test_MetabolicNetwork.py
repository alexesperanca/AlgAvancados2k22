# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from MetabolicNetwork import MetabolicNetwork

class TestMetabolicNetwork(unittest.TestCase):
    def setUp(self):
        self.m = MetabolicNetwork("metabolite-reaction")
        self.m.add_vertex_type("R1","reaction")
        self.m.add_vertex_type("R2","reaction")
        self.m.add_vertex_type("R3","reaction")
        self.m.add_vertex_type("M1","metabolite")
        self.m.add_vertex_type("M2","metabolite")
        self.m.add_vertex_type("M3","metabolite")
        self.m.add_vertex_type("M4","metabolite")
        self.m.add_vertex_type("M5","metabolite")
        self.m.add_vertex_type("M6","metabolite")
        self.m.add_edge("M1","R1")
        self.m.add_edge("M2","R1")
        self.m.add_edge("R1","M3")
        self.m.add_edge("R1","M4")
        self.m.add_edge("M4","R2")
        self.m.add_edge("M6","R2")
        self.m.add_edge("R2","M3")
        self.m.add_edge("M4","R3")
        self.m.add_edge("M5","R3")
        self.m.add_edge("R3","M6")
        self.m.add_edge("R3","M4")
        self.m.add_edge("R3","M5")
        self.m.add_edge("M6","R3")

        self.mrn = MetabolicNetwork("metabolite-reaction")
        self.mrn.load_from_file("example-net.txt")

        self.mmn = MetabolicNetwork("metabolite-metabolite")
        self.mmn.load_from_file("example-net.txt")

    
    def test_get_nodes_type(self):
        self.assertEqual(self.m.get_nodes_type("reaction"), ['R1', 'R2', 'R3'])
        self.assertEqual(self.m.get_nodes_type("metabolite"), ['M1', 'M2', 'M3', 'M4', 'M5', 'M6'])
        self.assertEqual(self.mrn.get_nodes_type("reaction"), ["R1", "R2", "R3"])
        self.assertEqual(self.mrn.get_nodes_type("metabolite"), ["M1", "M2", "M3", "M4", "M5", "M6"])

    def test_degrees_centrality(self):
        self.assertEqual(self.mmn.degrees_centrality(), {'M1': 0.4, 'M2': 0.4, 'M3': 0.8, 'M4': 1.0, 'M5': 0.4, 'M6': 0.6})
    
    def test_closeness_centrality(self):
        self.assertEqual(self.mmn.closeness_centrality(), {'M1': 0.4, 'M2': 0.4, 'M3': 0.0, 'M4': 0.36, 'M5': 0.4, 'M6': 0.6})

    def test_betweenness_centrality(self):
        self.assertEqual(self.mmn.betweenness_centrality('M6', 'M2'), {'M1': 0, 'M2': 0, 'M3': 0, 'M4': 1, 'M5': 0, 'M6': 0})

    def test_active_reactions(self):
        self.assertEqual(self.mrn.active_reactions(["M3", "M2"]), ['R1', 'R2', 'R3'])

    def test_met_prod(self):
        self.assertEqual(self.mrn.met_prod(["M4"]), ['M3', 'M4', 'M5', 'M6'])

    def test_final_met(self):
        self.assertEqual(self.mrn.final_met(["M5", "M2"]), ['M4', 'M5', 'M6', 'M3'])

if __name__ == '__main__':
   unittest.main()