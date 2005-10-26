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

from pychem.internal_coordinates import InternalCoordinatesCache
from pychem.molecular_graphs import BondSets, BendSets, DihedralSets, OutOfPlaneSets, CriteriaSet
from pychem.molecules import molecule_from_xyz_filename
from pychem.moldata import BOND_SINGLE
from pychem.units import from_angstrom

import unittest, math, copy, Numeric

__all__ = ["InternalCoordinatesTPA", "Chainrule"]


class InternalCoordinatesTPA(unittest.TestCase):        
    def setUp(self):
        self.molecule = molecule_from_xyz_filename("input/tpa_optimized.xyz")
        self.ic_cache = InternalCoordinatesCache(self.molecule)
        
        
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
        self.ic_cache.add_bond_lengths(bond_sets)
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
            self.ic_cache["All bonds"], 
            yield_alternatives
        )

    def test_bond_angles(self):
        bend_sets = BendSets([
            CriteriaSet("All bends", (None, None))
        ])
        self.ic_cache.add_bend_cosines(bend_sets)
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
            self.ic_cache["All bends"], 
            yield_alternatives
        )

    def test_dihedral_angles(self):
        dihedral_sets = DihedralSets([
            CriteriaSet("All dihedrals", (None, None)),
        ])
        self.ic_cache.add_dihedral_cosines(dihedral_sets)
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
            self.ic_cache["All dihedrals"], 
            yield_alternatives
        )     
            


class Chainrule(unittest.TestCase):        
    def setUp(self):
        self.ethene = molecule_from_xyz_filename("input/ethene.xyz")
        # Define the two (buggy) internal coordinates.
        self.ic_cache = InternalCoordinatesCache(self.ethene)
    
    def pair_test(self, internal_coordinate, ethene1, ethene2, expected_cos1, expected_cos2):
        test_cos1, gradient1 = internal_coordinate(ethene1.coordinates)
        test_cos2, gradient2 = internal_coordinate(ethene2.coordinates)
        if abs(test_cos1 - expected_cos1) > 1e-5:
            self.errors.append("Ethene1 problem: test cosine (%s) and expected cosine (%s) differ: %s" % (test_cos1, expected_cos1, test_cos1 - expected_cos1))
        if abs(test_cos2 - expected_cos2) > 1e-5:
            self.errors.append("Ethene2 problem: test cosine (%s) and expected cosine (%s) differ: %s" % (test_cos2, expected_cos2, test_cos2 - expected_cos2))
        delta = ethene2.coordinates - ethene1.coordinates
        gradient = 0.5*(gradient1+gradient2)
        delta_cos_estimate = Numeric.dot(Numeric.ravel(gradient), Numeric.ravel(delta))
        if abs(delta_cos_estimate - (expected_cos2 - expected_cos1)) > 1e-4:
            self.errors.append("Chain rule problem: delta_cos_estimate (%s) and delta_cos_expected (%s) differ: %s" % (delta_cos_estimate, (expected_cos2 - expected_cos1), delta_cos_estimate - (expected_cos2 - expected_cos1)))
        if abs(delta_cos_estimate - (test_cos2 - test_cos1)) > 1e-4:
            self.errors.append("Chain rule problem: delta_cos_estimate (%s) and delta_cos_test (%s) differ: %s" % (delta_cos_estimate, (test_cos2 - test_cos1), delta_cos_estimate - (test_cos2 - test_cos1)))

    def test_dihedral(self):
        self.ic_cache.add_dihedral_cosines(DihedralSets([CriteriaSet("HCCH", ((1, 6, 6, 1), None))]))
        dihedral_cos = self.ic_cache["HCCH"][0]

        self.errors = []
        
        def mutate_ethene(angle):
            result = copy.deepcopy(self.ethene)
            result.coordinates[1,1] = self.ethene.coordinates[1,1]*math.cos(angle)
            result.coordinates[1,2] = self.ethene.coordinates[1,1]*math.sin(angle)
            result.coordinates[2,1] = self.ethene.coordinates[2,1]*math.cos(angle)
            result.coordinates[2,2] = self.ethene.coordinates[2,1]*math.sin(angle)
            return result
        
        
        number = 100
        for index in xrange(number):
            angle1 = float(index)/number*2*math.pi
            angle2 = float(index+1)/number*2*math.pi
            ethene1 = mutate_ethene(angle1)
            ethene2 = mutate_ethene(angle2)
            self.pair_test(dihedral_cos, ethene1, ethene2, math.cos(angle1), math.cos(angle2))
        
        self.assertEqual(len(self.errors), 0, "\n".join(self.errors))
    
    def test_out_of_plane(self):
        self.ic_cache.add_out_of_plane_cosines(OutOfPlaneSets([CriteriaSet("CC(HCl)", ((6, 6, 1, 17), None))]))
        out_of_plane_cos = self.ic_cache["CC(HCl)"][0]

        self.errors = []
        
        def mutate_ethene(angle):
            result = copy.deepcopy(self.ethene)
            result.coordinates[1,0] = self.ethene.coordinates[1,0]*math.cos(angle)
            result.coordinates[1,2] = self.ethene.coordinates[1,0]*math.sin(angle)
            result.coordinates[2,0] = self.ethene.coordinates[2,0]*math.cos(angle)
            result.coordinates[2,2] = self.ethene.coordinates[2,0]*math.sin(angle)
            return result
        
        number = 50
        for index in xrange(number):
            angle1 = (float(index)/number-0.5)*math.pi
            angle2 = (float(index+1)/number-0.5)*math.pi
            ethene1 = mutate_ethene(angle1)
            ethene2 = mutate_ethene(angle2)
            self.pair_test(out_of_plane_cos, ethene1, ethene2, math.cos(angle1), math.cos(angle2))
        
        self.assertEqual(len(self.errors), 0, "\n".join(self.errors))            

