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

from molmod.molecular_graphs2 import *
from molmod.graphs2 import MatchGenerator, CriteriaSet
from molmod.molecules import molecule_xyz_from_filename
from molmod.data import BOND_SINGLE

import unittest, copy

__all__ = ["MolecularGraphs2"]


class MolecularGraphs2(unittest.TestCase):
    def load_graph(self, filename):
        self.molecule = molecule_xyz_from_filename(filename)
        self.molecular_graph = MolecularGraph(self.molecule)

    def verify(self, expected_results, test_results, yield_alternatives):
        for key in test_results.iterkeys():
            unsatisfied = expected_results[key]
            test = test_results[key]
            correct = []
            unexpected = []
            for test_item in test:
                item_correct = False
                for alternative in yield_alternatives(test_item):
                    if (alternative in unsatisfied):
                        correct.append(test_item)
                        unsatisfied.remove(alternative)
                        item_correct = True
                        break
                if not item_correct:
                    unexpected.append(test_item)
            message  = "Incorrect matches (%s):\n" % key
            message += "unexpected  (%i): %s\n" % (len(unexpected), unexpected)
            message += "unsatisfied (%i): %s\n" % (len(unsatisfied), unsatisfied)
            self.assertEqual(len(unexpected), 0, message)
            self.assertEqual(len(unsatisfied), 0, message)

    def test_bonds_tpa(self):
        self.load_graph("input/tpa.xyz")
        match_definition = BondMatchDefinition([
            CriteriaSet("HC", atom_criteria(1, 6)),
            CriteriaSet("CC", atom_criteria(6, 6)),
            CriteriaSet("CN", atom_criteria(6, 7)),
            CriteriaSet("C-sp3", atom_criteria(6, HasNumNeighbors(4))),
            CriteriaSet("C-[CN]", atom_criteria(6, MolecularOr(HasAtomNumber(6), HasAtomNumber(7)))),
            CriteriaSet("long", relation_criteria={frozenset([0,1]): BondLongerThan(2.1)}),
        ])
        expected_results = {
            "HC": set([(14, 1), (13, 1), (16, 2), (15, 2), (17, 3), (19, 3), (18, 3), (20, 4), (21, 4), (23, 5), (22, 5), (25, 6), (26, 6), (24, 6), (27, 7), (28, 7), (29, 8), (30, 8), (31, 9), (33, 9), (32, 9), (35, 10), (34, 10), (37, 11), (36, 11), (39, 12), (40, 12), (38, 12)]),
            "CC": set([(2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)]),
            "CN": set([(10, 0), (1, 0), (4, 0), (7, 0)]),
            "C-sp3": set([(10, 0), (1, 0), (4, 0), (7, 0), (2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)]),
            "C-[CN]": set([(10, 0), (1, 0), (4, 0), (7, 0), (2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)]),
            "long": set([(10, 0), (1, 0), (4, 0), (7, 0), (2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)]),
        }
        test_results = dict((key, []) for key in expected_results)
        match_generator = MatchGenerator(match_definition, debug=False)
        for match in match_generator(self.molecular_graph):
            test_results[match.tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))

        def yield_alternatives(test_item):
            yield test_item
            a, b = test_item
            if self.molecule.numbers[a] == self.molecule.numbers[b]:
                yield b, a

        self.verify(expected_results, test_results, yield_alternatives)

    def test_bending_angles_tpa(self):
        self.load_graph("input/tpa.xyz")
        match_definition = BendingAngleMatchDefinition([
            CriteriaSet("HCH", atom_criteria(1, 6, 1)),
            CriteriaSet("HCC", atom_criteria(1, 6, 6)),
            CriteriaSet("CCC", atom_criteria(6, 6, 6)),
            CriteriaSet("CCN", atom_criteria(6, 6, 7)),
            CriteriaSet("CNC", atom_criteria(6, 7, 6)),
            CriteriaSet("HCN", atom_criteria(1, 6, 7)),
        ])
        expected_results = {
            'HCC': set([(14, 1, 2), (13, 1, 2), (16, 2, 1), (15, 2, 1), (16, 2, 3), (15, 2, 3), (17, 3, 2), (19, 3, 2), (18, 3, 2), (20, 4, 5), (21, 4, 5), (23, 5, 6), (23, 5, 4), (22, 5, 6), (22, 5, 4), (25, 6, 5), (26, 6, 5), (24, 6, 5), (27, 7, 8), (28, 7, 8), (29, 8, 9), (30, 8, 9), (29, 8, 7), (30, 8, 7), (31, 9, 8), (33, 9, 8), (32, 9, 8), (35, 10, 11), (34, 10, 11), (37, 11, 12), (36, 11, 12), (37, 11, 10), (36, 11, 10), (39, 12, 11), (40, 12, 11), (38, 12, 11)]),
            'CCN': set([(2, 1, 0), (5, 4, 0), (8, 7, 0), (11, 10, 0)]),
            'CNC': set([(1, 0, 10), (4, 0, 10), (7, 0, 10), (4, 0, 1), (7, 0, 1), (7, 0, 4)]),
            'HCH': set([(13, 1, 14), (15, 2, 16), (19, 3, 17), (18, 3, 17), (18, 3, 19), (21, 4, 20), (22, 5, 23), (26, 6, 25), (24, 6, 25), (24, 6, 26), (28, 7, 27), (30, 8, 29), (33, 9, 31), (32, 9, 31), (32, 9, 33), (34, 10, 35), (36, 11, 37), (40, 12, 39), (38, 12, 39), (38, 12, 40)]),
            'HCN': set([(14, 1, 0), (13, 1, 0), (20, 4, 0), (21, 4, 0), (27, 7, 0), (28, 7, 0), (35, 10, 0), (34, 10, 0)]),
            'CCC': set([(3, 2, 1), (4, 5, 6), (7, 8, 9), (10, 11, 12)])
        }
        test_results = dict((key, []) for key in expected_results)
        match_generator = MatchGenerator(match_definition, debug=False)
        for match in match_generator(self.molecular_graph):
            test_results[match.tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))

        def yield_alternatives(test_item):
            yield test_item
            a, b, c = test_item
            if self.molecule.numbers[a] == self.molecule.numbers[c]:
                yield c, b, a

        self.verify(expected_results, test_results, yield_alternatives)

    def test_dihedral_angles_tpa(self):
        self.load_graph("input/tpa.xyz")
        match_definition = DihedralAngleMatchDefinition([
            CriteriaSet("HCCH", atom_criteria(1, 6, 6, 1)),
            CriteriaSet("HCCN", atom_criteria(1, 6, 6, 7)),
            CriteriaSet("HCNC", atom_criteria(1, 6, 7, 6)),
        ])
        expected_results = {
            'HCCH': set([(16, 2, 1, 14), (15, 2, 1, 14), (16, 2, 1, 13), (15, 2, 1, 13), (17, 3, 2, 16), (19, 3, 2, 16), (18, 3, 2, 16), (17, 3, 2, 15), (19, 3, 2, 15), (18, 3, 2, 15), (23, 5, 4, 20), (22, 5, 4, 20), (23, 5, 4, 21), (22, 5, 4, 21), (25, 6, 5, 23), (26, 6, 5, 23), (24, 6, 5, 23), (25, 6, 5, 22), (26, 6, 5, 22), (24, 6, 5, 22), (29, 8, 7, 27), (30, 8, 7, 27), (29, 8, 7, 28), (30, 8, 7, 28), (31, 9, 8, 29), (33, 9, 8, 29), (32, 9, 8, 29), (31, 9, 8, 30), (33, 9, 8, 30), (32, 9, 8, 30), (37, 11, 10, 35), (36, 11, 10, 35), (37, 11, 10, 34), (36, 11, 10, 34), (39, 12, 11, 37), (40, 12, 11, 37), (38, 12, 11, 37), (39, 12, 11, 36), (40, 12, 11, 36), (38, 12, 11, 36)]),
            'HCCN': set([(16, 2, 1, 0), (15, 2, 1, 0), (23, 5, 4, 0), (22, 5, 4, 0), (29, 8, 7, 0), (30, 8, 7, 0), (37, 11, 10, 0), (36, 11, 10, 0)]),
            'HCNC': set([(14, 1, 0, 10), (13, 1, 0, 10), (20, 4, 0, 10), (21, 4, 0, 10), (27, 7, 0, 10), (28, 7, 0, 10), (35, 10, 0, 1), (34, 10, 0, 1), (20, 4, 0, 1), (21, 4, 0, 1), (27, 7, 0, 1), (28, 7, 0, 1), (35, 10, 0, 4), (34, 10, 0, 4), (14, 1, 0, 4), (13, 1, 0, 4), (27, 7, 0, 4), (28, 7, 0, 4), (35, 10, 0, 7), (34, 10, 0, 7), (14, 1, 0, 7), (13, 1, 0, 7), (20, 4, 0, 7), (21, 4, 0, 7)])
        }
        test_results = dict((key, []) for key in expected_results)
        match_generator = MatchGenerator(match_definition, debug=False)
        for match in match_generator(self.molecular_graph):
            test_results[match.tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))

        def yield_alternatives(test_item):
            yield test_item
            a, b, c, d = test_item
            if (self.molecule.numbers[a] == self.molecule.numbers[d]) and (self.molecule.numbers[b] == self.molecule.numbers[c]):
                yield d, c, b, a

        self.verify(expected_results, test_results, yield_alternatives)

    def test_tetra_tpa(self):
        self.load_graph("input/tpa.xyz")
        match_definition = TetraMatchDefinition([
            CriteriaSet("C-(HCCH)", atom_criteria(6, 1, 6, 6, 1))
        ], node_tags={0: 1, 1: 1})
        expected_results = {
            'C-(HCCH)': set([
                ( 8, 29,  9,  7, 30), ( 8, 30,  9,  7, 29), ( 2,  1, 15, 16,  3), ( 2,  3, 15, 16,  1),
                (11, 10, 36, 37, 12), (11, 12, 36, 37, 10), ( 5,  4, 22, 23,  6), ( 5,  6, 22, 23,  4),
            ]),
        }
        test_results = dict((key, []) for key in expected_results)
        match_generator = MatchGenerator(match_definition, debug=False)
        for match in match_generator(self.molecular_graph):
            test_results[match.tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))


        def yield_alternatives(test_item):
            a, b, c, d, e = test_item
            yield test_item
            yield a, c, b, e, d
            if (self.molecule.numbers[b] == self.molecule.numbers[e]):
                yield a, e, c, d, b
                yield a, c, e, b, d
            if (self.molecule.numbers[c] == self.molecule.numbers[d]):
                yield a, b, d, c, e
                yield a, d, b, e, c
            if (self.molecule.numbers[b] == self.molecule.numbers[e]) and (self.molecule.numbers[c] == self.molecule.numbers[d]):
                yield a, e, d, c, b
                yield a, d, e, b, c

        self.verify(expected_results, test_results, yield_alternatives)

    def test_dihedral_angles_precursor(self):
        self.load_graph("input/precursor.xyz")
        self.molecular_graph.init_neighbors()
        match_definition = DihedralAngleMatchDefinition([
            CriteriaSet("all"),
        ])
        # construct all dihedral angles:
        all_dihedrals = set([])
        for b, c in self.molecular_graph.pairs:
            for a in self.molecular_graph.neighbors[b]:
                if a != c:
                    for d in self.molecular_graph.neighbors[c]:
                        if d != b:
                            all_dihedrals.add((a, b, c, d))
        expected_results = {
            'all': all_dihedrals,
        }
        test_results = dict((key, []) for key in expected_results)
        match_generator = MatchGenerator(match_definition, debug=False)
        for match in match_generator(self.molecular_graph):
            test_results[match.tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))

        def yield_alternatives(test_item):
            yield test_item
            a, b, c, d = test_item
            yield d, c, b, a

        self.verify(expected_results, test_results, yield_alternatives)
