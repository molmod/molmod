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


from molmod.graphs2 import OneToOne, Graph, MatchGenerator, ExactMatchFilter

import unittest, copy


__all__ = ["ExampleGraphs2"]


class ExampleGraphs2(unittest.TestCase):
    def setUp(self):
        self.graphs = [
            (
                "bond",
                Graph([(0,1)]),
                [(), ((0,1),)]
            ), (
                "angle",
                Graph([(0,1), (0,2)]),
                [(), ((1,2),)]
            ), (
                "star 3",
                Graph([(0,1), (0,2), (0,3)]),
                [(), ((1,2,3),), ((1,3,2),), ((1,2),), ((2,3),), ((1,3),)]
            ), (
                "triangle",
                Graph([(0,1), (1,2), (2,0)]),
                [(), ((0,1,2),), ((0,2,1),), ((0,1),), ((1,2),), ((0,2),)]
            ), (
                "chain 3",
                Graph([(0,1), (1,2), (2,3)]),
                [(), ((0,3), (1,2))]
            ), (
                "star 4",
                Graph([(0,1), (0,2), (0,3), (0,4)]),
                [(), ((1,2),), ((1,3),), ((1,4),), ((2,3),), ((2,4),), ((3,4),),
                 ((1,2),(3,4)), ((1,3),(2,4)), ((1,4),(2,3)),
                 ((1, 2, 3),), ((1, 3, 2),), ((1,2,4),), ((1,4,2),), ((1,3,4),),
                 ((1,4,3),), ((2,3,4),), ((2,4,3),), ((1,2,3,4),), ((1,2,4,3),),
                 ((1,3,2,4),), ((1,3,4,2),), ((1,4,2,3),), ((1,4,3,2),)]
            ), (
                "2 star 3",
                Graph([(0,1), (0,2), (0,3), (1,4), (1,5)]),
                [(),
                ((2,3),), ((4,5),), ((2,3),(4,5)),

                ((0,1), (2,5), (3,4)), ((0,1), (2,4), (3,5)),
                ((0,1), (2,5,3,4)), ((0,1), (2,4,3,5))
                ]
            ), (
                "4 star 3",
                Graph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9)]),
                [(),
                ((4,5),), ((6,7),), ((8,9),),
                ((4,5),(6,7)), ((4,5),(8,9)), ((6,7),(8,9)),
                ((4,5),(6,7),(8,9)),

                ((1,2),(4,6),(5,7)), ((1,2),(4,7),(5,6)),
                ((1,3),(4,8),(5,9)), ((1,3),(4,9),(5,8)),
                ((2,3),(6,8),(7,9)), ((2,3),(6,9),(7,8)),

                ((1,2),(4,6,5,7)), ((1,2),(4,7,5,6)),
                ((1,3),(4,8,5,9)), ((1,3),(4,9,5,8)),
                ((2,3),(6,8,7,9)), ((2,3),(6,9,7,8)),

                ((1,2),(4,6),(5,7),(8,9)), ((1,2),(4,7),(5,6),(8,9)),
                ((1,3),(4,8),(5,9),(6,7)), ((1,3),(4,9),(5,8),(6,7)),
                ((2,3),(4,5),(6,8),(7,9)), ((2,3),(4,5),(6,9),(7,8)),

                ((1,2),(4,6,5,7),(8,9)), ((1,2),(4,7,5,6),(8,9)),
                ((1,3),(4,8,5,9),(6,7)), ((1,3),(4,9,5,8),(6,7)),
                ((2,3),(4,5),(6,8,7,9)), ((2,3),(4,5),(6,9,7,8)),

                ((1,2,3),(4,6,8,5,7,9)), ((1,2,3),(4,6,9,5,7,8)),
                ((1,2,3),(4,7,8,5,6,9)), ((1,2,3),(4,7,9,5,6,8)),

                ((1,3,2),(4,8,6,5,9,7)), ((1,3,2),(4,8,7,5,9,6)),
                ((1,3,2),(4,9,6,5,8,7)), ((1,3,2),(4,9,7,5,8,6)),

                ((1,2,3),(4,6,8),(5,7,9)), ((1,2,3),(4,6,9),(5,7,8)),
                ((1,2,3),(4,7,8),(5,6,9)), ((1,2,3),(4,7,9),(5,6,8)),
                ((1,3,2),(4,8,6),(5,9,7)), ((1,3,2),(4,8,7),(5,9,6)),
                ((1,3,2),(4,9,6),(5,8,7)), ((1,3,2),(4,9,7),(5,8,6))
                ]
            )
        ]
        self.todo = [
            (
                "ortho hexane",
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
            )
        ]

    def tst_symmetries(self):
        for name, graph, expected_symmetries in self.graphs:
            common = []
            unexpected = []
            unsatisfied = copy.deepcopy(expected_symmetries)

            for symmetry in graph.symmetries:
                try:
                    i = unsatisfied.index(symmetry)
                    del unsatisfied[i]
                    common.append(symmetry)
                except:
                    unexpected.append(symmetry)

            def message():
                return ("-- Graphs %s ---\n" % name) + \
                       ("common (%i): %s\n" % (len(common), common)) + \
                       ("unexpected (%i): %s\n" % (len(unexpected), unexpected)) + \
                       ("unsatisfied (%i): %s\n" % (len(unsatisfied), unsatisfied))

            self.assert_(len(unexpected) == 0, message())
            self.assert_(len(unsatisfied) == 0, message())


    def test_match_generator(self):
        #print
        for name, graph, expected_symmetries in self.graphs:
            #print
            #print
            #print "GRAPH %s" % name
            match_generator = MatchGenerator(
                ExactMatchFilter(),
                graph, graph,
                debug=False
            )
            #print match_generator().next()
            for match in match_generator():
                pass
                #print "_ _ _ _ _", match, "_ _ _ _ _"
