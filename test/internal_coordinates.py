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

from pychem.internal_coordinates import Collection
from pychem.molecular_graphs import BondSets, BondAngleSets, DihedralAngleSets, CriteriaSet
from pychem.molecules import molecule_from_xyz_filename
from pychem.moldata import BOND_SINGLE
from pychem.units import from_angstrom

import unittest, math


__all__ = ["suite"]

suite = unittest.TestSuite()   

class TestInternalCoordinatesTPA(unittest.TestCase):        
    def setUp(self):
        self.molecule = molecule_from_xyz_filename("input/tpa_optimized.xyz")
        self.collection = Collection(self.molecule)
        
        
    def verify(self, expected_results, internal_coordinates, yield_alternatives):
        for internal_coordinate in internal_coordinates:
            value, derivates = internal_coordinate(self.molecule.coordinates)
            id = None
            for id in yield_alternatives(internal_coordinate.id):
                expected_value = expected_results.get(id)
                if expected_value != None:
                    break
            #if (expected_value == None):
            #    print id, None
            #elif (abs(value - expected_value) > 1e-10):
            #    print id, abs(value-expected_value)
            #else:
            #    print id, "OK"
            if (expected_value != None):
                self.assertAlmostEqual(value, expected_value, 3)

    def load_expected_results(self, filename, conversion):
        result = {}
        f = file(filename, 'r')
        for line in f:
            words = line.split()
            value = float(words[0])
            id = tuple([int(index) for index in words[1:]])
            result[id] = conversion(value)
        return result
            
    def test_bonds(self):
        bond_sets = BondSets([
            CriteriaSet("All bonds", (None, None)),
        ])
        self.collection.add_bond_lengths(bond_sets)
        expected_results = self.load_expected_results(
            "input/tpa_stretch.csv", 
            from_angstrom
        )
            
        def yield_alternatives(test_item):
            yield test_item
            a, b = test_item
            yield (b, a)
                
        self.verify(
            expected_results, 
            self.collection["All bonds"], 
            yield_alternatives
        )

    def test_bond_angles(self):
        bond_angle_sets = BondAngleSets([
            CriteriaSet("All bond angles", (None, None))
        ])
        self.collection.add_bond_cosines(bond_angle_sets)
        expected_results = self.load_expected_results(
            "input/tpa_bend.csv", 
            lambda x: math.cos(math.pi*x/180.0)
        )
        
        def yield_alternatives(test_item):
            yield test_item
            a, b, c = test_item
            yield (c, b, a)
                
        self.verify(
            expected_results, 
            self.collection["All bond angles"], 
            yield_alternatives
        )

    def test_dihedral_angles(self):
        dihedral_angle_sets = DihedralAngleSets([
            CriteriaSet("All dihedral angles", (None, None)),
        ])
        self.collection.add_dihedral_cosines(dihedral_angle_sets)
        expected_results = self.load_expected_results(
            "input/tpa_tors.csv", 
            lambda x: math.cos(math.pi*x/180.0)
        )
        def yield_alternatives(test_item):
            yield test_item
            a, b, c, d = test_item
            yield (d, c, b, a)
                
        self.verify(
            expected_results, 
            self.collection["All dihedral angles"], 
            yield_alternatives
        )     
            
            
suite.addTests([
    unittest.makeSuite(TestInternalCoordinatesTPA)
])
