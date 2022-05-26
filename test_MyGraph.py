# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

import unittest
from MyGraph import MyGraph

class TestMotifFinding(unittest.TestCase):
    def setUp(self):
        self.gr = MyGraph({1:{2: None}, 2:{3: None}, 3:{2: None, 4: None}, 4:{2: None}})

        self.gr1 = MyGraph()
        self.gr1.add_vertex(1)
        self.gr1.add_vertex(2)
        self.gr1.add_vertex(3)
        self.gr1.add_vertex(4)
            
        self.gr1.add_edge(1,2)
        self.gr1.add_edge(2,3)
        self.gr1.add_edge(3,2)
        self.gr1.add_edge(3,4)
        self.gr1.add_edge(4,2)

    def test_getnodes(self):
        self.assertEqual(self.gr.get_nodes(), [1, 2, 3, 4])
        self.assertEqual(self.gr1.get_nodes(), [1, 2, 3, 4])

    def test_getedges(self):
        self.assertEqual(self.gr.get_edges(), [(1, 2), (2, 3), (3, 2), (3, 4), (4, 2)])
        self.assertEqual(self.gr1.get_edges(), [(1, 2), (2, 3), (3, 2), (3, 4), (4, 2)])
    
    def test_get_successors(self):
        self.assertEqual(self.gr.get_successors(2), [3])
        self.assertEqual(self.gr1.get_successors(2), [3])

    def test_get_predecessors(self):
        self.assertEqual(self.gr.get_predecessors(2), [1, 3, 4])
        self.assertEqual(self.gr1.get_predecessors(2), [1, 3, 4])

    def test_get_adjacents(self):
        self.assertEqual(self.gr.get_adjacents(2), [1, 3, 4])
        self.assertEqual(self.gr1.get_adjacents(2), [1, 3, 4])

    def test_in_degree(self):
        self.assertEqual(self.gr.in_degree(2), 3)
        self.assertEqual(self.gr1.in_degree(2), 3)

    def test_out_degree(self):
        self.assertEqual(self.gr.out_degree(2), 1)
        self.assertEqual(self.gr1.out_degree(2), 1)

    def test_degree(self):
        self.assertEqual(self.gr.degree(2), 3)
        self.assertEqual(self.gr1.degree(2), 3)

    def test_distance(self):
        self.assertEqual(self.gr.distance(1, 4), 3)
        self.assertEqual(self.gr.distance(3, 4), 1)
        self.assertEqual(self.gr1.distance(1, 4), 3)
        self.assertEqual(self.gr1.distance(3, 4), 1)

    def test_shortest_path(self):
        self.assertEqual(self.gr.shortest_path(1, 4), "1 -> 2 -> 3 -> 4")
        self.assertEqual(self.gr.shortest_path(4, 3), "4 -> 3")

    def test_reachable_with_dist(self):
        self.assertEqual(self.gr.reachable_with_dist(1), [(2, 1), (3, 2), (4, 3)])
        self.assertEqual(self.gr.reachable_with_dist(3), [(2, 1), (4, 1)])
        self.assertEqual(self.gr1.reachable_with_dist(1), [(2, 1), (3, 2), (4, 3)])
        self.assertEqual(self.gr1.reachable_with_dist(3), [(2, 1), (4, 1)])

    def test_node_has_cycle(self):
        self.assertEqual(self.gr.node_has_cycle(2), True)
        self.assertEqual(self.gr.node_has_cycle(1), False)
        self.assertEqual(self.gr1.node_has_cycle(2), True)
        self.assertEqual(self.gr1.node_has_cycle(1), False)

    def test_has_cycle(self):
        self.assertEqual(self.gr.has_cycle(), True)
        self.assertEqual(self.gr1.has_cycle(), True)

if __name__ == '__main__':
    unittest.main()