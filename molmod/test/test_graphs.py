# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# MolMod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


from __future__ import print_function

from builtins import range
import copy
import unittest

import pkg_resources
import numpy as np

from molmod import *



__all__ = ["GraphTestCase"]


class Case(object):
    def __init__(self, name, edges, sc=None, rings=None):
        self.name = name
        self.graph = Graph(edges)
        if sc is None:
            self.symmetry_cycles = None
        else:
            self.symmetry_cycles = set([tuple(cycles) for cycles in sc])
        if rings is None:
            self.rings = None
        else:
            self.rings = set(rings)


class GraphTestCase(unittest.TestCase):
    def iter_cases(self, disconnected=False):
        if disconnected:
            yield Case("2chains", [(0,1),(1,2),(3,4),(4,5)], rings=[])
            yield Case("singletons", [(2,3),(5,6),(6,7)], rings=[])
        yield Case("bond", [(0,1)], sc=[[], [(0,1),]], rings=[])
        yield Case("angle", [(0,1), (0,2)], sc=[[], [(1,2),]], rings=[])
        yield Case("star 3", [(0,1), (0,2), (0,3)], sc=[
                         # Relation to symmetries in 2D figure:
            [],          # identity transformation
            [(1,2,3),],  # rotation in plane
            [(1,3,2),],  # rotation in plane in the other direction
            [(1,2),],    # mirror in plane
            [(2,3),],    # mirror in plane
            [(1,3),],    # mirror in plane
        ], rings=[])
        yield Case("triangle", [(0,1), (1,2), (2,0)], sc=[
            [],
            [(0,1,2),],
            [(0,2,1),],
            [(0,1),],
            [(0,2),],
            [(1,2),],
        ], rings=[(0,1,2)])
        yield Case("chain 3", [(0,1), (1,2), (2,3)], sc=[
            [],
            [(0,3), (1,2)],
        ], rings=[])
        yield Case("star 4", [(0,1), (0,2), (0,3), (0,4)], sc=[
            [],
            [(1,2)],
            [(1,3)],
            [(1,4)],
            [(2,3)],
            [(2,4)],
            [(3,4)],
            [(1,2), (3,4)],
            [(1,3), (2,4)],
            [(1,4), (2,3)],
            [(1, 2, 3)],
            [(1, 3, 2)],
            [(1,2,4)],
            [(1,4,2)],
            [(1,3,4)],
            [(1,4,3)],
            [(2,3,4)],
            [(2,4,3)],
            [(1,2,3,4)],
            [(1,2,4,3)],
            [(1,3,2,4)],
            [(1,3,4,2)],
            [(1,4,2,3)],
            [(1,4,3,2)],
        ], rings=[])
        yield Case("2 star 3", [(0,1), (0,2), (0,3), (1,4), (1,5)], sc=[
            [],
            [(2,3)],
            [(4,5)],
            [(2,3), (4,5)],
            [(0,1), (2,5), (3,4)],
            [(0,1), (2,4), (3,5)],
            [(0,1), (2,5,3,4)],
            [(0,1), (2,4,3,5)],
        ], rings=[])
        yield Case("4 star 3", [
            (0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8),
            (3, 9)
        ], sc=[
            [],
            [(4,5),],
            [(6,7),],
            [(8,9),],
            [(4,5),(6,7)],
            [(4,5),(8,9)],
            [(6,7),(8,9)],
            [(4,5),(6,7),(8,9)],

            [(1,2),(4,6),(5,7)],
            [(1,2),(4,7),(5,6)],
            [(1,3),(4,8),(5,9)],
            [(1,3),(4,9),(5,8)],
            [(2,3),(6,8),(7,9)],
            [(2,3),(6,9),(7,8)],

            [(1,2),(4,6,5,7)],
            [(1,2),(4,7,5,6)],
            [(1,3),(4,8,5,9)],
            [(1,3),(4,9,5,8)],
            [(2,3),(6,8,7,9)],
            [(2,3),(6,9,7,8)],

            [(1,2),(4,6),(5,7),(8,9)],
            [(1,2),(4,7),(5,6),(8,9)],
            [(1,3),(4,8),(5,9),(6,7)],
            [(1,3),(4,9),(5,8),(6,7)],
            [(2,3),(4,5),(6,8),(7,9)],
            [(2,3),(4,5),(6,9),(7,8)],

            [(1,2),(4,6,5,7),(8,9)],
            [(1,2),(4,7,5,6),(8,9)],
            [(1,3),(4,8,5,9),(6,7)],
            [(1,3),(4,9,5,8),(6,7)],
            [(2,3),(4,5),(6,8,7,9)],
            [(2,3),(4,5),(6,9,7,8)],

            [(1,2,3),(4,6,8,5,7,9)],
            [(1,2,3),(4,6,9,5,7,8)],
            [(1,2,3),(4,7,8,5,6,9)],
            [(1,2,3),(4,7,9,5,6,8)],

            [(1,3,2),(4,8,6,5,9,7)],
            [(1,3,2),(4,8,7,5,9,6)],
            [(1,3,2),(4,9,6,5,8,7)],
            [(1,3,2),(4,9,7,5,8,6)],

            [(1,2,3),(4,6,8),(5,7,9)],
            [(1,2,3),(4,6,9),(5,7,8)],
            [(1,2,3),(4,7,8),(5,6,9)],
            [(1,2,3),(4,7,9),(5,6,8)],
            [(1,3,2),(4,8,6),(5,9,7)],
            [(1,3,2),(4,8,7),(5,9,6)],
            [(1,3,2),(4,9,6),(5,8,7)],
            [(1,3,2),(4,9,7),(5,8,6)],
        ], rings=[])
        yield Case("pentagon", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], sc=[
            [],
            [(0,1,2,3,4)],
            [(0,2,4,1,3)],
            [(0,3,1,4,2)],
            [(0,4,3,2,1)],
            [(1,4),(2,3)],
            [(0,3),(1,2)],
            [(0,1),(2,4)],
            [(0,4),(1,3)],
            [(0,2),(3,4)],
        ], rings=[(0,1,2,3,4)])
        yield Case("orthobenzene", [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (3, 7)
        ], rings=[(0,1,2,3,4,5)])
        yield Case("square", [(0, 1), (1, 2), (2, 3), (3, 0)], sc=[
            [],
            [(0,1,2,3)],
            [(0,3,2,1)],
            [(0,2)],
            [(1,3)],
            [(0,1),(2,3)],
            [(0,2),(1,3)],
            [(0,3),(1,2)],
        ], rings=[(0,1,2,3)])
        yield Case("benzene", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
            rings=[(0,1,2,3,4,5)]
        )
        yield Case("0-2-4 hexane", [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (2, 7),
            (4,8)
        ], rings=[(0,1,2,3,4,5)])
        yield Case("ethene", [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)], sc=[
            [],
            [(2,3)],
            [(4,5)],
            [(2,3),(4,5)],
            [(0,1),(2,4),(3,5)],
            [(0,1),(2,5),(3,4)],
            [(0,1),(2,4,3,5)],
            [(0,1),(2,5,3,4)],
        ])
        yield Case("difficult", [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)], sc=[
            [],
            [(3,4)],
            [(5,6)],
            [(3,4),(5,6)],
            [(1,2),(3,5),(4,6)],
            [(1,2),(3,6),(4,5)],
            [(1,2),(3,5,4,6)],
            [(1,2),(3,6,4,5)],
        ])
        yield Case("naphthalene", [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (2, 6), (6, 7),
            (7, 8), (8, 9), (9, 3)
        ], rings=[(0,1,2,3,4,5),(2,3,9,8,7,6)])
        yield Case("cage", [
            (0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 7),
            (6, 7)
        ], rings=[])
        yield Case("tetraeder",
            [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)],
            rings=[(0,1,2),(0,1,3),(0,2,3),(1,2,3)]
        )
        yield Case("cube", [
            (0, 1), (0, 2), (0, 3), (1, 4), (1, 6), (2, 4), (2, 5), (3, 5),
            (3, 6), (4, 7), (5, 7), (6, 7)
        ],rings=[(0,1,4,2),(0,2,5,3),(0,1,6,3),(1,4,7,6),(3,5,7,6),(2,4,7,5)])
        yield Case("pentane", [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (0, 6), (1, 7),
            (1, 8), (2, 9), (2, 10), (3, 11), (3, 12), (4, 13), (4, 14),
            (5, 15), (5, 16)
        ], rings=[(0,1,2,3,4)])
        yield Case("weird", [
            (28,14), (28,4), (3,31), (32,10), (27,37), (27,38), (37,22),
            (33,15), (4,31), (24,39), (1,29), (32,22), (33,23), (26,36),
            (17,31), (2,38), (18,36), (33,42), (2,30), (33,12), (8,35),
            (29,5), (11,30), (32,14), (24,34), (25,35), (34,43), (29,15),
            (13,31), (32,40), (26,39), (16,30), (16,34), (0,36), (5,30),
            (3,39), (20,38), (36,23), (17,35), (34,6), (28,7), (25,38),
            (41,35), (0,28), (21,39), (9,29), (19,37), (1,37)
        ])

    def test_distances(self):
        # normal graph
        graph = Graph([(0,1), (1,2), (2,3), (3,4), (4,0), (2,5), (3,6), (3,7)])
        expecting = np.array([
            [0, 1, 2, 2, 1, 3, 3, 3],
            [1, 0, 1, 2, 2, 2, 3, 3],
            [2, 1, 0, 1, 2, 1, 2, 2],
            [2, 2, 1, 0, 1, 2, 1, 1],
            [1, 2, 2, 1, 0, 3, 2, 2],
            [3, 2, 1, 2, 3, 0, 3, 3],
            [3, 3, 2, 1, 2, 3, 0, 2],
            [3, 3, 2, 1, 2, 3, 2, 0],
        ],dtype=np.int32)
        self.assertEqual(expecting.shape,graph.distances.shape)
        self.assert_((expecting==graph.distances).all())
        # disconnected graph
        graph = Graph([(0,1), (1,2), (3,4)])
        expecting = np.array([
            [ 0, 1, 2, 0, 0],
            [ 1, 0, 1, 0, 0],
            [ 2, 1, 0, 0, 0],
            [ 0, 0, 0, 0, 1],
            [ 0, 0, 0, 1, 0],
        ],dtype=np.int32)
        self.assertEqual(expecting.shape,graph.distances.shape)
        self.assert_((expecting==graph.distances).all())
        # graph with singletons
        graph = Graph([(2,3), (3,4)])
        expecting = np.array([
            [ 0, 0, 0, 0, 0],
            [ 0, 0, 0, 0, 0],
            [ 0, 0, 0, 1, 2],
            [ 0, 0, 1, 0, 1],
            [ 0, 0, 2, 1, 0],
        ],dtype=np.int32)
        self.assertEqual(expecting.shape,graph.distances.shape)
        self.assert_((expecting==graph.distances).all())

    def test_neighbors(self):
        for case in self.iter_cases():
            g = case.graph
            counter = 0
            for central, neighbors in g.neighbors.items():
                for neighbor in neighbors:
                    counter += 1
                    self.assert_(frozenset([central,neighbor]) in g.edges)
            self.assertEqual(counter, len(g.edges)*2)

    def test_central_vertices(self):
        for case in self.iter_cases():
            g = case.graph
            max_distances = g.distances.max(axis=1)
            max_distances_min = max_distances[max_distances>0].min()
            self.assert_(len(g.central_vertices>0))
            for c in g.central_vertices:
                self.assert_(g.distances[c].max() == max_distances_min)
            self.assert_(g.central_vertex in g.central_vertices)

    def test_independent_vertices(self):
        edges = [(0,1), (0,2), (0,3), (1,4), (1,5), (6,7), (6,8), (6,9), (7,10), (7,11)]
        g = Graph(edges)
        self.assertEqual(g.independent_vertices, [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]])

    def test_fingerprints(self):
        for case in self.iter_cases():
            g0 = case.graph
            permutation = np.random.permutation(g0.num_vertices)
            new_edges = tuple((permutation[i], permutation[j]) for i,j in g0.edges)
            g1 = Graph(new_edges, g0.num_vertices)
            self.assert_((g0.fingerprint==g1.fingerprint).all())
            for i in range(g0.num_vertices):
                self.assert_((g0.vertex_fingerprints[i]==g1.vertex_fingerprints[permutation[i]]).all())

    def test_symmetries(self):
        cases = self.iter_cases()
        for case in cases:
            try:
                foo = case.graph.symmetry_cycles
                del foo
                if len(case.graph.independent_vertices) != 1:
                    self.fail("Sould have raised an error for disconnected graphs.")
            except ValueError:
                if len(case.graph.independent_vertices) == 1:
                    raise
                continue
            # generic tests:
            if True:
                for cycles in case.graph.symmetry_cycles:
                    if len(cycles) > 2:
                        # all cycles lengths must have a common divisor > 1
                        max_cycle_len = max(len(cycle) for cycle in cycles)
                        #print max_cycle_len
                        for divisor in range(2, max_cycle_len+1):
                            ok = True
                            for cycle in cycles:
                                if len(cycle)%divisor != 0:
                                    ok = False
                                    break
                            if ok:
                                break
                        if not ok:
                            self.fail_("The cycle lengths must have a common divisor. Not found in %s: %s" % (case.name, cycles))
            # checking expected data
            if case.symmetry_cycles is None: continue
            common = set([])
            unexpected = set([])
            unsatisfied = copy.deepcopy(case.symmetry_cycles)

            for cycles in case.graph.symmetry_cycles:
                if cycles in unsatisfied:
                    unsatisfied.remove(cycles)
                    common.add(cycles)
                else:
                    unexpected.add(cycles)

            def message():
                return ("-- Graphs %s ---\n" % case.name) + \
                       ("common (%i): %s\n" % (len(common), common)) + \
                       ("unexpected (%i): %s\n" % (len(unexpected), unexpected)) + \
                       ("unsatisfied (%i): %s\n" % (len(unsatisfied), unsatisfied))

            self.assert_(len(unexpected) == 0, message())
            self.assert_(len(unsatisfied) == 0, message())

    def test_equivalent_vertices(self):
        for case in self.iter_cases(disconnected=False):
            g = case.graph
            cycles = g.symmetry_cycles
            equivalent_vertices = {}
            for cycle in cycles:
                for sub in cycle:
                    for vertex in sub:
                        s = equivalent_vertices.setdefault(vertex, set([]))
                        s.update(sub)
            for vertex in range(g.num_vertices):
                equivalent_vertices.setdefault(vertex, set([vertex]))
            self.assertEqual(equivalent_vertices, g.equivalent_vertices)

    # auxiliary graph routines

    def test_iter_breadth_first(self):
        cases = self.iter_cases()
        for case in cases:
            g = case.graph
            start = np.random.randint(g.num_vertices)
            result = list(g.iter_breadth_first(start))
            if len(case.graph.independent_vertices) != 1:
                continue
            self.assertEqual(len(result), g.num_vertices)
            l_last = 0
            visited = np.zeros(g.num_vertices, int)
            for n,l in result:
                if l < l_last:
                    self.fail_("Distances in iter_breadth_first must be monotomically increasing.")
                self.assertEqual(l, g.distances[n,start])
                l_last = l
                visited[n] = 1
            self.assert_((visited==1).all())

    def test_iter_breadth_first_edges(self):
        cases = self.iter_cases()
        for case in cases:
            g = case.graph
            start = np.random.randint(g.num_vertices)
            result = list(g.iter_breadth_first_edges(start))
            if len(case.graph.independent_vertices) != 1:
                continue
            self.assertEqual(len(result), g.num_edges)
            distance_last = 0
            visited = np.zeros(g.num_edges, int)
            edge_index = dict((edge,i) for i,edge in enumerate(g.edges))
            for edge,distance,flag in result:
                if distance < distance_last:
                    self.fail_("Distances in iter_breadth_first_edge must be monotomically increasing.")
                self.assertEqual(distance, g.distances[edge[0],start])
                if flag:
                    self.assertEqual(distance, g.distances[edge[1],start])
                else:
                    self.assertEqual(distance+1, g.distances[edge[1],start])
                distance_last = distance
                visited[edge_index[frozenset(edge)]] = 1
            self.assert_((visited==1).all())

    def test_iter_shortest_paths(self):
        # a few exotic cases
        cases = [
            (
                Graph([(0,1),(0,2),(1,3),(2,3),(3,4),(3,5),(4,6),(5,6)]),
                0,6,5,set([(0,1,3,4,6),(0,1,3,5,6),(0,2,3,4,6),(0,2,3,5,6)])
            ),
            (
                Graph([(0,1),(1,2),(2,3),(0,4),(4,5),(5,3),(3,6)]),
                0,6,5,set([(0,1,2,3,6),(0,4,5,3,6)])
            ),
            (
                Graph([(0,1),(0,2),(1,3),(2,4),(3,4),(3,5),(4,5)]),
                0,5,4,set([(0,1,3,5),(0,2,4,5)])
            ),
            (
                Graph([(0,1),(0,2),(1,3),(2,3),(1,4),(2,5),(3,6),(4,6),(5,6)]),
                0,6,4,set([(0,1,4,6),(0,1,3,6),(0,2,3,6),(0,2,5,6)])
            ),
        ]
        for graph, begin, end, length, expected_paths in cases:
            #print graph.edges
            for path in graph.iter_shortest_paths(begin,end):
                self.assertEqual(len(path), length)
                self.assert_(path in expected_paths)
                expected_paths.discard(path)
                #print path
            self.assertEqual(len(expected_paths), 0)

    def test_get_subgraph(self):
        for case in self.iter_cases():
            #print case.name
            graph = case.graph
            # non-normalized case
            if graph.num_edges<3:
                continue
            while True:
                subvertices = np.random.permutation(graph.num_vertices)
                subvertices = subvertices[:graph.num_vertices-1]
                try:
                    subgraph = graph.get_subgraph(subvertices)
                    break
                except GraphError:
                    pass
            self.assertEqual(graph.num_vertices, subgraph.num_vertices)
            for i,j in subgraph.edges:
                self.assert_(i in subvertices)
                self.assert_(j in subvertices)
            for i_new, i_old in enumerate(subgraph._old_edge_indexes):
                self.assertEqual(subgraph.edges[i_new], graph.edges[i_old])

            # normalized case
            subgraph = graph.get_subgraph(subvertices, normalize=True)
            self.assert_(subgraph.num_vertices <= len(subvertices))
            for p0,p1 in enumerate(subgraph._old_edge_indexes):
                i0,j0 = subgraph.edges[p0]
                #i1,j1 = subgraph.edges[p1]
                edge_old = frozenset([subgraph._old_vertex_indexes[i0],subgraph._old_vertex_indexes[j0]])
                self.assertEqual(edge_old, graph.edges[p1])

    def test_halfs(self):
        # Tests both get_halfs and get_halfs_double
        edges1 = [(0,1), (1,2), (2,3), (3,4), (4,0), (2,5), (3,6), (3,7)]
        graph1 = Graph(edges1)
        edges2 = [(0,1), (1,2), (2,3), (1,4), (1,5), (5,6), (3,7), (7,8), (8,9), (9,3)]
        graph2 = Graph(edges2)

        try:
            graph1.get_halfs(2,3)
            self.fail("graph1.get_halfs(2,3) should have raised a GraphError.")
        except GraphError:
            pass

        try:
            graph2.get_halfs(7,8)
            self.fail("graph2.get_halfs(7,8) should have raised a GraphError.")
        except GraphError:
            pass

        part1, part2 = graph2.get_halfs(2,3)
        self.assertEqual(part1, set([0,1,2,4,5,6]))
        self.assertEqual(part2, set([3,7,8,9]))
        part1, part2 = graph2.get_halfs(5,6)
        self.assertEqual(part1, set([0,1,2,3,4,5,7,8,9]))
        self.assertEqual(part2, set([6]))
        part1, part2 = graph2.get_halfs(1,2)
        self.assertEqual(part1, set([0,1,4,5,6]))
        self.assertEqual(part2, set([2,3,7,8,9]))
        part1, part2 = graph1.get_halfs(3,7)
        self.assertEqual(part1, set([0,1,2,3,4,5,6]))
        self.assertEqual(part2, set([7]))

        part1, part2, hinges = graph1.get_halfs_double(0,1,2,3)
        self.assertEqual(part1, set([0,4,3,6,7]))
        self.assertEqual(part2, set([1,2,5]))
        self.assertEqual(hinges, (0,1,3,2))
        part1, part2, hinges = graph2.get_halfs_double(3,7,8,9)
        self.assertEqual(part1, set([0,1,2,3,4,5,6,9]))
        self.assertEqual(part2, set([7,8]))
        self.assertEqual(hinges, (3,7,9,8))

        try:
            part1, part2, hinges = graph1.get_halfs_double(0,2,1,3)
            self.fail("graph1.get_halfs_double(0,2,1,3) should raise a GraphError")
        except GraphError:
            pass

        try:
            part1, part2, hinges = graph1.get_halfs_double(0,1,2,5)
            self.fail("graph1.get_halfs_double(0,1,2,5) should raise a GraphError")
        except GraphError:
            pass

        try:
            part1, part2, hinges = graph2.get_halfs_double(0,1,2,3)
            self.fail("graph2.get_halfs_double(0,1,2,3) should raise a GraphError")
        except GraphError:
            pass

        try:
            part1, part2, hinges = graph2.get_halfs_double(1,0,2,3)
            self.fail("graph2.get_halfs_double(1,0,2,3) should raise a GraphError")
        except GraphError:
            pass

        edges3 = [
            (9,4),(8,1),(3,12),(2,6),(1,4),(3,4),(2,3),
            (1,7),(2,5),(11,3),(10,4),(0,2),(0,1)
        ]
        graph3 = Graph(edges3)

        part1, part2, hinges = graph3.get_halfs_double(1,4,2,3)
        self.assertEqual(part1, set([0,1,2,5,6,7,8]))
        self.assertEqual(part2, set([3,4,9,10,11,12]))
        self.assertEqual(hinges, (1,4,2,3))

    # match generator related tests

    def check_graph_search(self, pattern, verbose=False, debug=False, callback=None):
        if verbose: print()
        graph_search = GraphSearch(
            pattern,
            debug=debug
        )
        cases = self.iter_cases()
        for case in cases:
            if verbose: print()
            if verbose: print()
            if verbose: print("GRAPH %s" % case.name)
            matches = []
            for match in graph_search(case.graph):
                matches.append(match)
                if verbose:
                    print("_ _ _ _ _", match, "_ _ _ _ _")
            if callback is not None:
                callback(case, matches)

    def test_ring_pattern(self):
        def callback(case, matches):
            if case.rings is None:
                return
            self.assertEqual(len(case.rings), len(matches))
            for match in matches:
                self.assert_(match.ring_vertices in case.rings)
        self.check_graph_search(RingPattern(10), callback=callback)

    def test_custom_pattern(self):
        # just run through the code
        for case in self.iter_cases():
            try:
                graph_search = GraphSearch(CustomPattern(case.graph))
                if len(case.graph.independent_vertices) != 1:
                    self.fail("Sould have raised an error for disconnected graphs.")
            except SubgraphPatternError:
                if len(case.graph.independent_vertices) == 1:
                    raise

            next(graph_search(case.graph))

    def test_start_vertex(self):
        subject_graph = Graph([(0,1),(1,2),(1,3),(2,4)])
        pattern_graph = Graph([(0,1),(1,2)])
        pattern0 = CustomPattern(pattern_graph)
        matches0 = list(GraphSearch(pattern0)(subject_graph))
        pattern1 = CustomPattern(pattern_graph, start_vertex=0)
        matches1 = list(GraphSearch(pattern1)(subject_graph))
        self.assertEqual(len(matches0), len(matches1))
        # compare all matches
        for match0 in matches0:
            found = False
            for match1 in matches1:
                if match0.forward == match1.forward:
                    found = True
                    break
            self.assert_(found)

    def test_pattern_lower_symmetry(self):
        subject_graph = Graph([
            (0,1),(1,2),(2,3),(3,4),(4,5),(5,0),
            (0,6),(1,7),(2,8),(3,9),(4,10),(5,11),
        ])
        pattern_graph = Graph([(0,1),(0,2),(0,3)])
        pattern = CustomPattern(pattern_graph, vertex_tags={1:1})
        collection = {}
        for match in GraphSearch(pattern)(subject_graph):
            a = match.forward[0]
            b = match.forward[1]
            s = collection.setdefault(a, set([]))
            s.add(b)
        for bs in collection.values():
            self.assertEqual(len(bs), 3)

    def test_rings_zeolite(self):
        cases = [
            ("opt_5ring10T.xyz", (10, 10, 12)),
            ("opt_5ring12T.xyz", (8, 10, 10, 10, 10)),
            ("opt_6_6ring.xyz", (12, 12)),
            ("opt_6ring.xyz", (12,)),
        ]
        gs = GraphSearch(RingPattern(12), debug=False)
        for fn_xyz, expected_sizes in cases:
            mol = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/%s" % fn_xyz))
            mol.set_default_graph()
            sizes = []
            for match in gs(mol.graph):
                sizes.append(len(match))
            sizes.sort()
            sizes = tuple(sizes)
            self.assertEqual(sizes, expected_sizes)
