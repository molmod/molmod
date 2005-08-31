# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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


from pychem.graphs import SymmetricGraph, MatchFilterMolecular, OneToOne, Criterium
from pychem.molecules import molecule_from_xyz

import unittest, copy


__all__ = ["suite"]

suite = unittest.TestSuite()        


class TestExampleGraphs(unittest.TestCase):
    def setUp(self):
        self.graphs = [
            (   
                "bond", 
                SymmetricGraph([(0,1)]), 
                [(), ((0,1),)]
            ), (
                "angle",
                SymmetricGraph([(0,1), (0,2)]),
                [(), ((1,2),)]
            ), (
                "star 3", 
                SymmetricGraph([(0,1), (0,2), (0,3)]),
                [(), ((1,2,3),), ((1,3,2),), ((1,2),), ((2,3),), ((1,3),)]
            ), (
                "triangle", 
                SymmetricGraph([(0,1), (1,2), (2,0)]),
                [(), ((0,1,2),), ((0,2,1),), ((0,1),), ((1,2),), ((0,2),)]
            ), (
                "chain 3", 
                SymmetricGraph([(0,1), (1,2), (2,3)]),
                [(), ((0,3), (1,2))]
            ), (
                "star 4", 
                SymmetricGraph([(0,1), (0,2), (0,3), (0,4)]),
                [(), ((1,2),), ((1,3),), ((1,4),), ((2,3),), ((2,4),), ((3,4),),
                 ((1,2),(3,4)), ((1,3),(2,4)), ((1,4),(2,3)),
                 ((1, 2, 3),), ((1, 3, 2),), ((1,2,4),), ((1,4,2),), ((1,3,4),),
                 ((1,4,3),), ((2,3,4),), ((2,4,3),), ((1,2,3,4),), ((1,2,4,3),),
                 ((1,3,2,4),), ((1,3,4,2),), ((1,4,2,3),), ((1,4,3,2),)]
            ), (
                "2 star 3", 
                SymmetricGraph([(0,1), (0,2), (0,3), (1,4), (1,5)]),
                [(),
                ((2,3),), ((4,5),), ((2,3),(4,5)), 
                
                ((0,1), (2,5), (3,4)), ((0,1), (2,4), (3,5)), 
                ((0,1), (2,5,3,4)), ((0,1), (2,4,3,5))
                ]
            ), (
                "4 star 3", 
                SymmetricGraph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9)]),
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
                SymmetricGraph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (3, 7)]),
                []
            ), (
                "square",
                SymmetricGraph([(0, 1), (1, 2), (2, 3), (3, 0)]),
                []
            ), (
                "benzene",
                SymmetricGraph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
                []
            ), (
                "0-2-4 hexane",
                SymmetricGraph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (2, 7), (4,8)]),
                []
            ), (
                "ethene",
                SymmetricGraph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]),
                []
            ), (
                "difficult",
                SymmetricGraph([(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]),
                []
            ), (
                "naphthalene",
                SymmetricGraph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (2, 6), (6, 7), (7, 8), (8, 9), (9, 3)]),
                []
            ), (
                "cage",
                SymmetricGraph([(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 7), (6, 7)]),
                []
            ), (
                "tetraeder",
                SymmetricGraph([(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]),
                []
            ), (
                "cube",
                SymmetricGraph([(0, 1), (0, 2), (0, 3), (1, 4), (1, 6), (2, 4), (2, 5), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]),
                []
            )
        ]
        
    def test_symmetries(self):
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


class TestOneToOne(unittest.TestCase):
    def test_multiplication(self):
        A = OneToOne([(1, "a"), (4, "z"), (2, "b")])
        B = OneToOne([("a", 7), ("z", 2), ("b", 0)])
        C = OneToOne([(1, 7), (4, 2), (2, 0)])
        self.assertEqual((B*A).forward, C.forward)


class AtomRequire(Criterium):
    def __init__(self, number, molecule):
        self.number = number
        self.molecule = molecule
        Criterium.__init__(self, number)
        
    def __call__(self, index):
        return self.molecule.numbers[index] == self.number


class AtomDeny(Criterium):
    def __init__(self, number, molecule):
        self.number = number
        self.molecule = molecule
        Criterium.__init__(self, number)
        
    def __call__(self, index):
        return self.molecule.numbers[index] != self.number


class TestMolecularGraphs(unittest.TestCase):        
    def test_tpa(self):
        molecule = molecule_from_xyz("input/tpa.xyz")
        graph, bonds = molecule.get_graph()
        
        subgraph = SymmetricGraph([(0, 1), (0, 2), (0, 3)], 0)
        graph_filter = MatchFilterMolecular(
            subgraph, 
            {0: 0, 1: 1, 1: 1, 3: 1},
            bonds,
            atom_criteria = {0: AtomRequire(6, molecule),
                             1: AtomRequire(1, molecule),
                             2: AtomRequire(1, molecule),
                             3: AtomRequire(1, molecule)}
        )
        
        for match in subgraph.yield_matching_subgraphs(graph):
            for parsed in graph_filter.parse(match):
                pass
                
        
suite.addTests([
    unittest.makeSuite(TestExampleGraphs),
    unittest.makeSuite(TestOneToOne)
])
