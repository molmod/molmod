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

from pychem.internal_coordinates import InternalCoordinatesCache, Select, Delta, Dot, Cos, Mul, Sub, Distance, DistanceSqr, Sqrt, Div, Sqr, Scale
from pychem.molecular_graphs import BondSets, BendSets, DihedralSets, OutOfPlaneSets, CriteriaSet
from pychem.molecules import molecule_from_xyz_filename
from pychem.moldata import BOND_SINGLE
from pychem.units import from_angstrom

import unittest, math, copy, Numeric, RandomArray

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
        self.ic_cache = InternalCoordinatesCache(self.ethene)

    def create_external_basis(self, coordinates):
        x_rotation = Numeric.zeros((3,3), Numeric.Float)
        x_rotation[1,2] = -1
        x_rotation[2,1] = 1
        y_rotation = Numeric.zeros((3,3), Numeric.Float)
        y_rotation[2,0] = -1
        y_rotation[0,2] = 1
        z_rotation = Numeric.zeros((3,3), Numeric.Float)
        z_rotation[0,1] = -1
        z_rotation[1,0] = 1
        return Numeric.array([
            [1, 0, 0]*len(coordinates),
            [0, 1, 0]*len(coordinates),
            [0, 0, 1]*len(coordinates),
            Numeric.ravel(Numeric.dot(coordinates, x_rotation)),
            Numeric.ravel(Numeric.dot(coordinates, y_rotation)),
            Numeric.ravel(Numeric.dot(coordinates, z_rotation))
        ], Numeric.Float)
    
    def check_sanity(self, coordinates, tangent, internal_coordinate):
        external_basis = Numeric.transpose(self.create_external_basis(coordinates))
        components = Numeric.dot(Numeric.ravel(tangent), external_basis)
        self.assert_(sum(components**2) < 1e-10, "Chain rule problem: failed on sanity_check for '%s' (%s)" % (internal_coordinate.tag, components))
        #print internal_coordinate.tag, sum(Numeric.dot(Numeric.ravel(tangent), external_basis)**2)
            
    def pair_test(self, internal_coordinate, ethene1, ethene2, expected_cos1, expected_cos2):
        test_cos1, tangent1 = internal_coordinate(ethene1.coordinates)
        test_cos2, tangent2 = internal_coordinate(ethene2.coordinates)
        
        # validate values
        if abs(test_cos1 - expected_cos1) > 1e-10:
            self.errors.append("Ethene1 problem: test cosine (%s) and expected cosine (%s) differ: %s" % (test_cos1, expected_cos1, test_cos1 - expected_cos1))
        if abs(test_cos2 - expected_cos2) > 1e-10:
            self.errors.append("Ethene2 problem: test cosine (%s) and expected cosine (%s) differ: %s" % (test_cos2, expected_cos2, test_cos2 - expected_cos2))

        # validate tangents
        delta = ethene2.coordinates - ethene1.coordinates
        tangent = 0.5*(tangent1+tangent2)
        delta_cos_estimate = Numeric.dot(Numeric.ravel(tangent), Numeric.ravel(delta))
        if abs(delta_cos_estimate - (expected_cos2 - expected_cos1)) > 1e-4:
            self.errors.append("Chain rule problem: delta_cos_estimate (%s) and delta_cos_expected (%s) differ: %s" % (delta_cos_estimate, (expected_cos2 - expected_cos1), delta_cos_estimate - (expected_cos2 - expected_cos1)))
        if abs(delta_cos_estimate - (test_cos2 - test_cos1)) > 1e-4:
            self.errors.append("Chain rule problem: delta_cos_estimate (%s) and delta_cos_test (%s) differ: %s" % (delta_cos_estimate, (test_cos2 - test_cos1), delta_cos_estimate - (test_cos2 - test_cos1)))

        # sanity check on tangents
        self.check_sanity(ethene1.coordinates, tangent1, internal_coordinate)
        self.check_sanity(ethene2.coordinates, tangent2, internal_coordinate)

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

    def load_internal_coordinates(self):
        self.ic_cache.add_out_of_plane_cosines(OutOfPlaneSets([CriteriaSet("CC(HCl)", ((6, 6, 1, 17), None))]))
        self.ic_cache.add_dihedral_cosines(DihedralSets([CriteriaSet("HCCH", ((1, 6, 6, 1), None))]))
        self.ic_cache.add_bond_lengths(BondSets([CriteriaSet("All bonds", (None, None))]))
        
        result = []
        # a) test wether the source is ok
        result.append(self.ic_cache["All bonds"][0])
        d1, d2 = self.ic_cache["All bonds"][0:2]
        s0 = self.ic_cache.add(Select, 0)
        s1 = self.ic_cache.add(Select, 1)
        s2 = self.ic_cache.add(Select, 2)
        e1 = self.ic_cache.add(Delta, s1, s0)
        e2 = self.ic_cache.add(Delta, s1, s2)        
        # b) part tests
        # b.1) unary, protocal A
        ic = self.ic_cache.add(Sqrt, d1)
        self.ic_cache.add_internal_coordinate("sqrt d1", ic)
        result.append(ic)
        ic = self.ic_cache.add(Sqr, d1)
        self.ic_cache.add_internal_coordinate("sqr d1", ic)
        result.append(ic)
        # b.2) binary, protocal A
        ic = self.ic_cache.add(Mul, d1, d2)
        self.ic_cache.add_internal_coordinate("mul d1 d2", ic)
        result.append(ic)
        ic = self.ic_cache.add(Sub, d1, d2)
        self.ic_cache.add_internal_coordinate("sub d1 d2", ic)
        result.append(ic)
        ic = self.ic_cache.add(Div, d1, d2)
        self.ic_cache.add_internal_coordinate("div d1 d2", ic)
        result.append(ic)
        # b.3) unary, protocal B
        ic = self.ic_cache.add(DistanceSqr, e1)
        self.ic_cache.add_internal_coordinate("distance_sqr e1", ic)
        result.append(ic)
        # b.4) binary, protocal B
        temp_ic = self.ic_cache.add(Sub, e1, e2)
        ic = self.ic_cache.add(DistanceSqr, temp_ic)
        self.ic_cache.add_internal_coordinate("distance_sqr sub e1 e2", ic)
        result.append(ic)
        ic = self.ic_cache.add(Dot, e1, e2)
        self.ic_cache.add_internal_coordinate("dot e1 e2", ic)
        result.append(ic)
        ic = self.ic_cache.add(Cos, d1, d2)
        self.ic_cache.add_internal_coordinate("cos d1 d2", ic)
        result.append(ic)
        # b.5) exotics
        temp_ic = self.ic_cache.add(Scale, d1, e1)
        ic = self.ic_cache.add(DistanceSqr, temp_ic)
        self.ic_cache.add_internal_coordinate("distance_sqr scale d1 e1", ic)
        result.append(ic)
        # c) application tests
        result.append(self.ic_cache["CC(HCl)"][0])
        result.append(self.ic_cache["HCCH"][0])
        return result

    def sanity_test(self, internal_coordinate, mod_ethene):
        foo, tangent = internal_coordinate(mod_ethene.coordinates)
        self.check_sanity(mod_ethene.coordinates, tangent, internal_coordinate)

    def test_sanity_random(self):
        internal_coordinates = self.load_internal_coordinates()

        def mutate_ethene():
            result = copy.deepcopy(self.ethene)
            result.coordinates += RandomArray.uniform(-3, 3, result.coordinates.shape)
            return result
        
        for internal_coordinate in internal_coordinates:
            for index in xrange(100):
                mod_ethene = mutate_ethene()
                self.sanity_test(internal_coordinate, mod_ethene)

    def consistent_test(self, internal_coordinate, ethene1, ethene2):
        value1, tangent1 = internal_coordinate(ethene1.coordinates)
        value2, tangent2 = internal_coordinate(ethene2.coordinates)
        
        # validate tangents
        delta = ethene2.coordinates - ethene1.coordinates
        tangent = 0.5*(tangent1+tangent2)
        delta_value_estimate = Numeric.dot(Numeric.ravel(tangent), Numeric.ravel(delta))
        if abs(delta_value_estimate - (value2 - value1)) > 1e-10:
            self.errors.append(
                "Chain rule problem: delta_value_estimate (%s) and value2 - value1 (%s) differ: %s" % (
                    delta_value_estimate, 
                    (value2 - value1), 
                    delta_value_estimate - (value2 - value1)
                )
            )

    def test_random(self):
        internal_coordinates = self.load_internal_coordinates()

        def mutate_ethene(strength, ethene):
            result = copy.deepcopy(ethene)
            result.coordinates += RandomArray.uniform(-strength, strength, result.coordinates.shape)
            return result
        
        self.errors = []
        
        for internal_coordinate in internal_coordinates:
            for index in xrange(100):
                mod1_ethene = mutate_ethene(3, self.ethene)
                mod2_ethene = mutate_ethene(1e-4, mod1_ethene)
                self.consistent_test(internal_coordinate, mod1_ethene, mod2_ethene)

        self.assertEqual(len(self.errors), 0, "\n".join(self.errors))            
