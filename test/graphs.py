# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


from molmod.graphs import OneToOne, Graph, MatchGenerator, EgoMatchDefinition,\
    RingMatchDefinition, SubgraphMatchDefinition

import unittest, copy


__all__ = ["GraphsTestCase"]


class Case(object):
    def __init__(self, name, pairs, symmetries=None):
        self.name = name
        self.graph = Graph(set([frozenset([a,b]) for a,b in pairs]))
        if symmetries is None:
            self.symmetries = None
        else:
            self.symmetries = set([frozenset([frozenset(group) for group in symmetry]) for symmetry in symmetries])


class GraphsTestCase(unittest.TestCase):
    def setUp(self):
        self.cases = [
            Case("bond", [(0,1)], [[], [(0,1),]]),
            Case("angle", [(0,1), (0,2)], [[], [(1,2),]]),
            Case("star 3", [(0,1), (0,2), (0,3)], [
                    [],
                    [(1,2,3),],
                    [(1,3,2),],
                    [(1,2),],
                    [(2,3),],
                    [(1,3),],
            ]), Case("triangle", [(0,1), (1,2), (2,0)], [
                    [],
                    [(0,1,2),],
                    [(0,2,1),],
                    [(0,1),],
                    [(0,2),],
                    [(1,2),],
            ]), Case("chain 3", [(0,1), (1,2), (2,3)], [
                    [],
                    [(0,3), (1,2)],
            ]), Case("star 4", [(0,1), (0,2), (0,3), (0,4)], [
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
            ]), Case("2 star 3", [(0,1), (0,2), (0,3), (1,4), (1,5)], [
                    [],
                    [(2,3)],
                    [(4,5)],
                    [(2,3), (4,5)],
                    [(0,1), (2,5), (3,4)],
                    [(0,1), (2,4), (3,5)],
                    [(0,1), (2,5,3,4)],
                    [(0,1), (2,4,3,5)],
            ]), Case("4 star 3", [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9)], [
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
            ]),
            Case("pentagon", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
            Case("ortho benzene", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (3, 7)]),
            Case("square", [(0, 1), (1, 2), (2, 3), (3, 0)]),
            Case("benzene", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
            Case("0-2-4 hexane", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (2, 7), (4,8)]),
            Case("ethene", [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]),
            Case("difficult", [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]),
            Case("naphthalene", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (2, 6), (6, 7), (7, 8), (8, 9), (9, 3)]),
            Case("cage", [(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 7), (6, 7)]),
            Case("tetraeder", [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]),
            Case("cube", [(0, 1), (0, 2), (0, 3), (1, 4), (1, 6), (2, 4), (2, 5), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]),
            Case("pentane", [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (0, 6), (1, 7), (1, 8), (2, 9), (2, 10), (3, 11), (3, 12), (4, 13), (4, 14), (5, 15), (5, 16)]),
        ]

    def test_symmetries(self):
        for case in self.cases:
            if case.symmetries is None: continue
            common = set([])
            unexpected = set([])
            unsatisfied = copy.deepcopy(case.symmetries)

            case.graph.init_symmetries()
            for cycles in case.graph.symmetry_cycles:
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
        for case in self.cases:
            if verbose: print
            if verbose: print
            if verbose: print "GRAPH %s" % case.name
            matches = []
            for match in match_generator(case.graph):
                matches.append(match)
                if verbose: print "_ _ _ _ _", match, "_ _ _ _ _"
            if callback is not None:
                callback(case.name, case.graph, matches)

    def test_ego_match_definition(self):
        def callback(name, graph, matches):
            self.assert_(len(matches) > 0, "Expected at least one match (graph=%s), got %i." % (name, len(matches)))

        self.do_match_generator_test(EgoMatchDefinition(), callback=callback)

    def test_ring_match_definition(self):
        self.do_match_generator_test(RingMatchDefinition(10))

    def test_subgraph_match_definition(self):
        for case in self.cases:
            match_generator = MatchGenerator(SubgraphMatchDefinition(case.graph))
            match_generator(case.graph).next()

    def test_symmetries(self):
        pairs = [(28,14), (28,4), (3,31), (32,10), (27,37), (27,38), (37,22), (33,15), (4,31), (24,39), (1,29), (32,22), (33,23), (26,36), (17,31), (2,38), (18,36), (33,42), (2,30), (33,12), (8,35), (29,5), (11,30), (32,14), (24,34), (25,35), (34,43), (29,15), (13,31), (32,40), (26,39), (16,30), (16,34), (0,36), (5,30), (3,39), (20,38), (36,23), (17,35), (34,6), (28,7), (25,38), (41,35), (0,28), (21,39), (9,29), (19,37), (1,37)]
        g = Graph(set([frozenset([a,b]) for a,b in pairs]))
        for match in MatchGenerator(EgoMatchDefinition(), debug=False)(g):
            pass

    def test_nodes_per_independent_graph(self):
        pairs = [(0,1), (0,2), (0,3), (1,4), (1,5), (6,7), (6,8), (6,9), (7,10), (7,11)]
        g = Graph(
            set([frozenset([a,b]) for a,b in pairs]),
            [0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9, 11],
        )
        result = g.get_nodes_per_independent_graph()
        self.assert_(result == [[0, 2, 4, 1, 3, 5], [6, 8, 10, 7, 9, 11]])


