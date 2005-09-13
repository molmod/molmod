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

from pychem.molecular_graphs import *
from pychem.molecules import molecule_from_xyz_filename
from pychem.moldata import BOND_SINGLE

import unittest, copy


__all__ = ["suite"]

suite = unittest.TestSuite()   

class TestMolecularGraphsTPA(unittest.TestCase):        
    def setUp(self):
        self.molecule = molecule_from_xyz_filename("input/tpa.xyz")
        self.molecular_graph = MolecularGraph(self.molecule)
        
    def verify(self, expected_results, test_results, generate_alternatives):
        for key in test_results.iterkeys():
            unsatisfied = expected_results[key]
            test = test_results[key]
            correct = []
            unexpected = []
            for test_item in test:
                alternatives = generate_alternatives(test_item)
                item_correct = False
                for alternative in alternatives:
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
            
    def test_bonds(self):
        criteria_sets = BondSets([
            CriteriaSet("HC", ((1, 6), None)),
            CriteriaSet("CC", ((6, 6), None)),
            CriteriaSet("CN", ((6, 7), None)),
            CriteriaSet("C-sp3", ((6, None), None), ({1: NumNeighboursRequire(4)}, None), False),
            CriteriaSet("C-[CN]", ((6, None), None), ({1: MolecularOr(AtomNumberRequire(6), AtomNumberRequire(7))}, None), False)
        ])
        expected_results = {
            "HC": set([(14, 1), (13, 1), (16, 2), (15, 2), (17, 3), (19, 3), (18, 3), (20, 4), (21, 4), (23, 5), (22, 5), (25, 6), (26, 6), (24, 6), (27, 7), (28, 7), (29, 8), (30, 8), (31, 9), (33, 9), (32, 9), (35, 10), (34, 10), (37, 11), (36, 11), (39, 12), (40, 12), (38, 12)]),
            "CC": set([(2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)]), 
            "CN": set([(10, 0), (1, 0), (4, 0), (7, 0)]), 
            "C-sp3": set([(10, 0), (1, 0), (4, 0), (7, 0), (2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)]), 
            "C-[CN]": set([(10, 0), (1, 0), (4, 0), (7, 0), (2, 1), (3, 2), (5, 4), (6, 5), (8, 7), (9, 8), (11, 10), (12, 11)])
        }
        test_results = dict((key, []) for key in expected_results)
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            test_results[tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))
            
        def generate_alternatives(test_item):
            a, b = test_item
            if self.molecule.numbers[a] == self.molecule.numbers[b]:
                return [test_item, (b, a)]
            else:
                return [test_item]
                
        self.verify(expected_results, test_results, generate_alternatives)

    def test_bond_angles(self):
        criteria_sets = BondAngleSets([
            CriteriaSet("HCH", ((1, 6, 1), None)),
            CriteriaSet("HCC", ((1, 6, 6), None)),
            CriteriaSet("CCC", ((6, 6, 6), None)),
            CriteriaSet("CCN", ((6, 6, 7), None)),
            CriteriaSet("CNC", ((6, 7, 6), None)),
            CriteriaSet("HCN", ((1, 6, 7), None))
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
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            test_results[tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))
        
        def generate_alternatives(test_item):
            a, b, c = test_item
            if self.molecule.numbers[a] == self.molecule.numbers[c]:
                return [test_item, (c, b, a)]
            else:
                return [test_item]
                
        self.verify(expected_results, test_results, generate_alternatives)        

    def test_dihedral_angles(self):
        criteria_sets = DihedralAngleSets([
            CriteriaSet("HCCH", ((1, 6, 6, 1), None)),
            CriteriaSet("HCCN", ((1, 6, 6, 7), None)),
            CriteriaSet("HCNC", ((1, 6, 7, 6), None))
        ])
        expected_results = {
            'HCCH': set([(16, 2, 1, 14), (15, 2, 1, 14), (16, 2, 1, 13), (15, 2, 1, 13), (17, 3, 2, 16), (19, 3, 2, 16), (18, 3, 2, 16), (17, 3, 2, 15), (19, 3, 2, 15), (18, 3, 2, 15), (23, 5, 4, 20), (22, 5, 4, 20), (23, 5, 4, 21), (22, 5, 4, 21), (25, 6, 5, 23), (26, 6, 5, 23), (24, 6, 5, 23), (25, 6, 5, 22), (26, 6, 5, 22), (24, 6, 5, 22), (29, 8, 7, 27), (30, 8, 7, 27), (29, 8, 7, 28), (30, 8, 7, 28), (31, 9, 8, 29), (33, 9, 8, 29), (32, 9, 8, 29), (31, 9, 8, 30), (33, 9, 8, 30), (32, 9, 8, 30), (37, 11, 10, 35), (36, 11, 10, 35), (37, 11, 10, 34), (36, 11, 10, 34), (39, 12, 11, 37), (40, 12, 11, 37), (38, 12, 11, 37), (39, 12, 11, 36), (40, 12, 11, 36), (38, 12, 11, 36)]),
            'HCCN': set([(16, 2, 1, 0), (15, 2, 1, 0), (23, 5, 4, 0), (22, 5, 4, 0), (29, 8, 7, 0), (30, 8, 7, 0), (37, 11, 10, 0), (36, 11, 10, 0)]),
            'HCNC': set([(14, 1, 0, 10), (13, 1, 0, 10), (20, 4, 0, 10), (21, 4, 0, 10), (27, 7, 0, 10), (28, 7, 0, 10), (35, 10, 0, 1), (34, 10, 0, 1), (20, 4, 0, 1), (21, 4, 0, 1), (27, 7, 0, 1), (28, 7, 0, 1), (35, 10, 0, 4), (34, 10, 0, 4), (14, 1, 0, 4), (13, 1, 0, 4), (27, 7, 0, 4), (28, 7, 0, 4), (35, 10, 0, 7), (34, 10, 0, 7), (14, 1, 0, 7), (13, 1, 0, 7), (20, 4, 0, 7), (21, 4, 0, 7)])
        }
        test_results = dict((key, []) for key in expected_results)
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            test_results[tag].append(tuple(match.get_destination(index) for index in xrange(len(match))))
        
        def generate_alternatives(test_item):
            a, b, c, d = test_item
            if (self.molecule.numbers[a] == self.molecule.numbers[d]) and (self.molecule.numbers[b] == self.molecule.numbers[c]):
                return [test_item, (d, c, b, a)]
            else:
                return [test_item]
                
        self.verify(expected_results, test_results, generate_alternatives)     
            
            
suite.addTests([
    unittest.makeSuite(TestMolecularGraphsTPA)
])
