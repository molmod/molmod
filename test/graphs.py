# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from molmod.graphs import OneToOne, Graph, MatchGenerator, EgoMatchDefinition,\
    RingMatchDefinition, SubgraphMatchDefinition

import unittest, copy


__all__ = ["GraphsTestCase"]


class GraphsTestCase(unittest.TestCase):
    def setUp(self):
        self.graphs = [
            (
                "bond",
                Graph([(0,1)]),
                set([
                    frozenset([]),
                    frozenset([(0,1),])
                ])
            ), (
                "angle",
                Graph([(0,1), (0,2)]),
                set([
                    frozenset([]),
                    frozenset([(1,2),])
                ])
            ), (
                "star 3",
                Graph([(0,1), (0,2), (0,3)]),
                set([
                    frozenset([]),
                    frozenset([(1,2,3),]),
                    frozenset([(1,3,2),]),
                    frozenset([(1,2),]),
                    frozenset([(2,3),]),
                    frozenset([(1,3),]),
                ])
            ), (
                "triangle",
                Graph([(0,1), (1,2), (2,0)]),
                set([
                    frozenset([]),
                    frozenset([(0,1,2),]),
                    frozenset([(0,2,1),]),
                    frozenset([(0,1),]),
                    frozenset([(0,2),]),
                    frozenset([(1,2),]),
                ])
            ), (
                "chain 3",
                Graph([(0,1), (1,2), (2,3)]),
                set([
                    frozenset([]),
                    frozenset([(0,3), (1,2)]),
                ])
            ), (
                "star 4",
                Graph([(0,1), (0,2), (0,3), (0,4)]),
                set([
                    frozenset([]),
                    frozenset([(1,2)]),
                    frozenset([(1,3)]),
                    frozenset([(1,4)]),
                    frozenset([(2,3)]),
                    frozenset([(2,4)]),
                    frozenset([(3,4)]),
                    frozenset([(1,2), (3,4)]),
                    frozenset([(1,3), (2,4)]),
                    frozenset([(1,4), (2,3)]),
                    frozenset([(1, 2, 3)]),
                    frozenset([(1, 3, 2)]),
                    frozenset([(1,2,4)]),
                    frozenset([(1,4,2)]),
                    frozenset([(1,3,4)]),
                    frozenset([(1,4,3)]),
                    frozenset([(2,3,4)]),
                    frozenset([(2,4,3)]),
                    frozenset([(1,2,3,4)]),
                    frozenset([(1,2,4,3)]),
                    frozenset([(1,3,2,4)]),
                    frozenset([(1,3,4,2)]),
                    frozenset([(1,4,2,3)]),
                    frozenset([(1,4,3,2)])
                ])
            ), (
                "2 star 3",
                Graph([(0,1), (0,2), (0,3), (1,4), (1,5)]),
                set([
                    frozenset([]),
                    frozenset([(2,3)]),
                    frozenset([(4,5)]),
                    frozenset([(2,3), (4,5)]),
                    frozenset([(0,1), (2,5), (3,4)]),
                    frozenset([(0,1), (2,4), (3,5)]),
                    frozenset([(0,1), (2,5,3,4)]),
                    frozenset([(0,1), (2,4,3,5)])
                ])
            ), (
                "4 star 3",
                Graph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9)]),
                set([
                   frozenset([]),
                   frozenset([(4,5),]),
                   frozenset([(6,7),]),
                   frozenset([(8,9),]),
                   frozenset([(4,5),(6,7)]),
                   frozenset([(4,5),(8,9)]),
                   frozenset([(6,7),(8,9)]),
                   frozenset([(4,5),(6,7),(8,9)]),

                   frozenset([(1,2),(4,6),(5,7)]),
                   frozenset([(1,2),(4,7),(5,6)]),
                   frozenset([(1,3),(4,8),(5,9)]),
                   frozenset([(1,3),(4,9),(5,8)]),
                   frozenset([(2,3),(6,8),(7,9)]),
                   frozenset([(2,3),(6,9),(7,8)]),

                   frozenset([(1,2),(4,6,5,7)]),
                   frozenset([(1,2),(4,7,5,6)]),
                   frozenset([(1,3),(4,8,5,9)]),
                   frozenset([(1,3),(4,9,5,8)]),
                   frozenset([(2,3),(6,8,7,9)]),
                   frozenset([(2,3),(6,9,7,8)]),

                   frozenset([(1,2),(4,6),(5,7),(8,9)]),
                   frozenset([(1,2),(4,7),(5,6),(8,9)]),
                   frozenset([(1,3),(4,8),(5,9),(6,7)]),
                   frozenset([(1,3),(4,9),(5,8),(6,7)]),
                   frozenset([(2,3),(4,5),(6,8),(7,9)]),
                   frozenset([(2,3),(4,5),(6,9),(7,8)]),

                   frozenset([(1,2),(4,6,5,7),(8,9)]),
                   frozenset([(1,2),(4,7,5,6),(8,9)]),
                   frozenset([(1,3),(4,8,5,9),(6,7)]),
                   frozenset([(1,3),(4,9,5,8),(6,7)]),
                   frozenset([(2,3),(4,5),(6,8,7,9)]),
                   frozenset([(2,3),(4,5),(6,9,7,8)]),

                   frozenset([(1,2,3),(4,6,8,5,7,9)]),
                   frozenset([(1,2,3),(4,6,9,5,7,8)]),
                   frozenset([(1,2,3),(4,7,8,5,6,9)]),
                   frozenset([(1,2,3),(4,7,9,5,6,8)]),

                   frozenset([(1,3,2),(4,8,6,5,9,7)]),
                   frozenset([(1,3,2),(4,8,7,5,9,6)]),
                   frozenset([(1,3,2),(4,9,6,5,8,7)]),
                   frozenset([(1,3,2),(4,9,7,5,8,6)]),

                   frozenset([(1,2,3),(4,6,8),(5,7,9)]),
                   frozenset([(1,2,3),(4,6,9),(5,7,8)]),
                   frozenset([(1,2,3),(4,7,8),(5,6,9)]),
                   frozenset([(1,2,3),(4,7,9),(5,6,8)]),
                   frozenset([(1,3,2),(4,8,6),(5,9,7)]),
                   frozenset([(1,3,2),(4,8,7),(5,9,6)]),
                   frozenset([(1,3,2),(4,9,6),(5,8,7)]),
                   frozenset([(1,3,2),(4,9,7),(5,8,6)])
                ])
            )
        ]
        self.todo = [
            (
                "pentagon",
                Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
                []
            ), (
                "ortho benzene",
                Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (3, 7)]),
                []
            ), (
                "square",
                Graph([(0, 1), (1, 2), (2, 3), (3, 0)]),
                []
            ), (
                "benzene",
                Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
                []
            ), (
                "0-2-4 hexane",
                Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (2, 7), (4,8)]),
                []
            ), (
                "ethene",
                Graph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]),
                []
            ), (
                "difficult",
                Graph([(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]),
                []
            ), (
                "naphthalene",
                Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (2, 6), (6, 7), (7, 8), (8, 9), (9, 3)]),
                []
            ), (
                "cage",
                Graph([(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 7), (6, 7)]),
                []
            ), (
                "tetraeder",
                Graph([(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]),
                []
            ), (
                "cube",
                Graph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 6), (2, 4), (2, 5), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]),
                []
            ), (
                "pentane",
                Graph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (0, 6), (1, 7), (1, 8), (2, 9), (2, 10), (3, 11), (3, 12), (4, 13), (4, 14), (5, 15), (5, 16)]),
                []
            )
        ]

    def test_symmetries(self):
        for name, graph, expected_symmetries in self.graphs:
            common = set([])
            unexpected = set([])
            unsatisfied = copy.deepcopy(expected_symmetries)

            graph.init_symmetries()
            for cycles in graph.symmetry_cycles:
                if cycles in unsatisfied:
                    unsatisfied.remove(cycles)
                    common.add(cycles)
                else:
                    unexpected.add(cycles)

            def message():
                return ("-- Graphs %s ---\n" % name) + \
                       ("common (%i): %s\n" % (len(common), common)) + \
                       ("unexpected (%i): %s\n" % (len(unexpected), unexpected)) + \
                       ("unsatisfied (%i): %s\n" % (len(unsatisfied), unsatisfied))

            self.assert_(len(unexpected) == 0, message())
            self.assert_(len(unsatisfied) == 0, message())

    def do_match_generator_test(self, match_definition, verbose=False, debug=False, callback=None):
        if verbose: print
        match_generator = MatchGenerator(
            match_definition,
            debug=debug
        )
        for name, graph, foo in self.graphs + self.todo:
            if verbose: print
            if verbose: print
            if verbose: print "GRAPH %s" % name
            matches = []
            for match in match_generator(graph):
                matches.append(match)
                if verbose: print "_ _ _ _ _", match, "_ _ _ _ _"
            if callback is not None:
                callback(name, graph, matches)

    def test_ego_match_definition(self):
        def callback(name, graph, matches):
            self.assert_(len(matches) > 0, "Expected at least one match (graph=%s), got %i." % (name, len(matches)))

        self.do_match_generator_test(EgoMatchDefinition(), callback=callback)

    def test_ring_match_definition(self):
        self.do_match_generator_test(RingMatchDefinition(10))

    def test_subgraph_match_definition(self):
        for name, graph, foo in self.graphs + self.todo:
            match_generator = MatchGenerator(SubgraphMatchDefinition(graph))
            match_generator(graph).next()

    def test_symmetries(self):
        g = Graph(set([frozenset([28, 14]), frozenset([4, 28]), frozenset([20, 38]), frozenset([3, 31]), frozenset([32, 10]), frozenset([27, 38]), frozenset([37, 22]), frozenset([17, 31]), frozenset([4, 31]), frozenset([24, 39]), frozenset([1, 29]), frozenset([32, 22]), frozenset([33, 23]), frozenset([26, 36]), frozenset([33, 15]), frozenset([2, 38]), frozenset([18, 36]), frozenset([33, 42]), frozenset([2, 30]), frozenset([33, 12]), frozenset([8, 35]), frozenset([29, 5]), frozenset([11, 30]), frozenset([32, 14]), frozenset([24, 34]), frozenset([1, 37]), frozenset([25, 35]), frozenset([34, 43]), frozenset([29, 15]), frozenset([13, 31]), frozenset([32, 40]), frozenset([26, 39]), frozenset([16, 30]), frozenset([16, 34]), frozenset([41, 35]), frozenset([0, 36]), frozenset([5, 30]), frozenset([3, 39]), frozenset([27, 37]), frozenset([36, 23]), frozenset([17, 35]), frozenset([34, 6]), frozenset([28, 7]), frozenset([21, 39]), frozenset([0, 28]), frozenset([9, 29]), frozenset([19, 37]), frozenset([25, 38])]))
        for match in MatchGenerator(EgoMatchDefinition(), debug=False)(g):
            pass

    def test_nodes_per_independent_graph(self):
        g = Graph(
            [(0,1), (0,2), (0,3), (1,4), (1,5), (6,7), (6,8), (6,9), (7,10), (7,11)],
            [0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9, 11],
        )
        result = g.get_nodes_per_independent_graph()
        self.assert_(result == [[0, 2, 4, 1, 3, 5], [6, 8, 10, 7, 9, 11]])
