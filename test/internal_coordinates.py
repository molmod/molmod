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
# --


from molmod.internal_coordinates import InternalCoordinatesCache, Select, Delta, Dot, Mul, Sub, Distance, DistanceSqr, Sqrt, Div, Sqr, Scale
from molmod.molecular_graphs import generate_molecular_graph, atom_criteria
from molmod.graphs import CriteriaSet
from molmod.data.bonds import BOND_SINGLE
from molmod.units import angstrom

from ccio.xyz import XYZFile

import unittest, math, copy, numpy, sys


__all__ = ["InternalCoordinatesTestCase", "ChainruleTestCase"]


class InternalCoordinatesTestCase(unittest.TestCase):
    def setUp(self):
        self.molecule = XYZFile("input/tpa_optimized.xyz").get_molecule()
        graph = generate_molecular_graph(self.molecule)
        self.ic_cache = InternalCoordinatesCache(graph)

    def verify(self, expected_results, internal_coordinates, yield_alternatives):
        for internal_coordinate in internal_coordinates:
            value, derivates = internal_coordinate.value_tangent(self.molecule.coordinates)
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
        self.ic_cache.add_bond_lengths([CriteriaSet(tag="All bonds")])
        expected_results = self.load_expected_results(
            "input/tpa_stretch.csv",
            (lambda x: x*angstrom),
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
        self.ic_cache.add_bending_cosines([CriteriaSet(tag="All bends")])
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
        self.ic_cache.add_dihedral_cosines([CriteriaSet(tag="All dihedrals")])
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



class ChainruleTestCase(unittest.TestCase):
    def setUp(self):
        self.ethene = XYZFile("input/ethene.xyz").get_molecule()
        graph = generate_molecular_graph(self.ethene)
        self.ic_cache = InternalCoordinatesCache(graph)

    def create_external_basis(self, coordinates):
        x_rotation = numpy.zeros((3,3), float)
        x_rotation[1,2] = -1
        x_rotation[2,1] = 1
        y_rotation = numpy.zeros((3,3), float)
        y_rotation[2,0] = -1
        y_rotation[0,2] = 1
        z_rotation = numpy.zeros((3,3), float)
        z_rotation[0,1] = -1
        z_rotation[1,0] = 1
        return numpy.array([
            [1, 0, 0]*len(coordinates),
            [0, 1, 0]*len(coordinates),
            [0, 0, 1]*len(coordinates),
            numpy.ravel(numpy.dot(coordinates, x_rotation)),
            numpy.ravel(numpy.dot(coordinates, y_rotation)),
            numpy.ravel(numpy.dot(coordinates, z_rotation))
        ], float)

    def check_sanity(self, coordinates, tangent, internal_coordinate):
        external_basis = numpy.transpose(self.create_external_basis(coordinates))
        components = numpy.dot(numpy.ravel(tangent), external_basis)
        self.assert_(sum(components**2) < 1e-10, "Chain rule problem: failed on sanity_check for '%s' (%s)" % (internal_coordinate.tag, components))
        #print internal_coordinate.tag, sum(numpy.dot(numpy.ravel(tangent), external_basis)**2)

    def pair_test(self, internal_coordinate, ethene1, ethene2, expected_value1, expected_value2):
        test_value1, tangent1 = internal_coordinate.value_tangent(ethene1.coordinates)
        test_value2, tangent2 = internal_coordinate.value_tangent(ethene2.coordinates)

        # validate values
        if abs(test_value1 - expected_value1) > 1e-10:
            self.errors.append("Ethene1 problem: test value (%s) and expected value (%s) differ: %s" % (test_value1, expected_value1, test_value1 - expected_value1))
        if abs(test_value2 - expected_value2) > 1e-10:
            self.errors.append("Ethene2 problem: test value (%s) and expected value (%s) differ: %s" % (test_value2, expected_value2, test_value2 - expected_value2))

        # validate tangents
        delta = ethene2.coordinates - ethene1.coordinates
        tangent = 0.5*(tangent1+tangent2)
        delta_value_estimate = numpy.dot(numpy.ravel(tangent), numpy.ravel(delta))
        if abs(delta_value_estimate - (expected_value2 - expected_value1)) > 1e-4:
            self.errors.append("Chain rule problem: delta_value_estimate (%s) and delta_value_expected (%s) differ: %s" % (delta_value_estimate, (expected_value2 - expected_value1), delta_value_estimate - (expected_value2 - expected_value1)))
        if abs(delta_value_estimate - (test_value2 - test_value1)) > 1e-4:
            self.errors.append("Chain rule problem: delta_value_estimate (%s) and delta_value_test (%s) differ: %s" % (delta_value_estimate, (test_value2 - test_value1), delta_value_estimate - (test_value2 - test_value1)))

        # sanity check on tangents
        self.check_sanity(ethene1.coordinates, tangent1, internal_coordinate)
        self.check_sanity(ethene2.coordinates, tangent2, internal_coordinate)

    def test_dihedral(self):
        self.ic_cache.add_dihedral_cosines([CriteriaSet(atom_criteria(1, 6, 6, 1), tag="HCCH")])
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

    def test_out_of_plane_cosine(self):
        self.ic_cache.add_out_of_plane_cosines([CriteriaSet(atom_criteria(6, 6, 1, 17), tag="CC(HCl)")])
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

    def test_out_of_plane_distance(self):
        self.ic_cache.add_out_of_plane_distances([CriteriaSet(atom_criteria(6, 6, 1, 17), tag="C(CHCl)")])
        out_of_plane_distance = self.ic_cache["C(CHCl)"][0]

        self.errors = []

        def mutate_ethene(distance):
            result = copy.deepcopy(self.ethene)
            result.coordinates[0,2] = distance
            return result

        number = 50
        for index in xrange(number):
            distance1 = float(index+0.1)/number
            distance2 = float(index+0.3)/number
            ethene1 = mutate_ethene(distance1)
            ethene2 = mutate_ethene(distance2)
            self.pair_test(out_of_plane_distance, ethene1, ethene2, distance1, distance2)

        self.assertEqual(len(self.errors), 0, "\n".join(self.errors))

    def load_internal_coordinates(self):
        self.ic_cache.add_out_of_plane_cosines([CriteriaSet(atom_criteria(6, 6, 1, 17), tag="CC(HCl)")])
        self.ic_cache.add_dihedral_cosines([CriteriaSet(atom_criteria(1, 6, 6, 1), tag="HCCH")])
        self.ic_cache.add_bond_lengths([CriteriaSet(tag="All bonds")])

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
        ic = self.ic_cache.add(DistanceSqr, e1)
        self.ic_cache.add_internal_coordinate("e1*e1", ic)
        result.append(ic)
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
        foo, tangent = internal_coordinate.value_tangent(mod_ethene.coordinates)
        self.check_sanity(mod_ethene.coordinates, tangent, internal_coordinate)

    def test_sanity_random(self):
        internal_coordinates = self.load_internal_coordinates()

        def mutate_ethene():
            result = copy.deepcopy(self.ethene)
            result.coordinates += numpy.random.uniform(-3, 3, result.coordinates.shape)
            return result

        for internal_coordinate in internal_coordinates:
            for index in xrange(100):
                mod_ethene = mutate_ethene()
                self.sanity_test(internal_coordinate, mod_ethene)

    def tangent_test(self, icn, index, largest, internal_coordinate, ethene1, ethene2):
        dof = len(ethene1.coordinates)*3
        value1, tangent1 = internal_coordinate.value_tangent(ethene1.coordinates)
        value2, tangent2 = internal_coordinate.value_tangent(ethene2.coordinates)
        tangent1.shape = (dof,)
        tangent2.shape = (dof,)

        delta = ethene2.coordinates - ethene1.coordinates
        delta.shape = (dof,)
        tangent = 0.5*(tangent1 + tangent2)

        # validate tangents
        delta_value_estimate = numpy.dot(tangent, delta)
        error = abs(delta_value_estimate - (value2 - value1))
        if error > 1e-8:
            self.errors.append(
                "%2i  %2i/%2i %s Chain rule problem: delta_value_estimate1 (%10.5f) and value2 - value1 (%10.5f) differ: % 3.2e" % (
                    icn, index, largest,
                    internal_coordinate.label(),
                    delta_value_estimate1,
                    (value2 - value1),
                    error
                )
            )

    def curvature_test(self, icn, index, largest, internal_coordinate, ethene):
        dof = len(ethene.coordinates)*3
        value, tangent = internal_coordinate.value_tangent(ethene.coordinates)
        tangent.shape = (dof,)
        curvature = internal_coordinate.curvature(ethene.coordinates)
        curvature.shape = (dof, dof)

        def str_matrix(matrix):
            result = ""
            for row in matrix:
                for element in row:
                    if element == 0.0:
                        result += " ."
                    else:
                        log10_element = int(math.log10(abs(element)))
                        if log10_element >= 10:
                            result += "XX"
                        elif log10_element <= -10:
                            result += ",,"
                        else:
                            result += "%+2i" % log10_element
                result += "\n"
            return result

        epsilon = 1e-5
        numeric_curvature = numpy.zeros(curvature.shape, float)
        #print internal_coordinate.label()
        for i in xrange(dof):
            excenter_coordinates = ethene.coordinates.copy()
            excenter_coordinates[i/3, i%3] += epsilon
            excenter_value, excenter_tangent = internal_coordinate.value_tangent(excenter_coordinates)
            excenter_tangent.shape = (dof,)
            numeric_curvature[i,:] += 0.5*(excenter_tangent - tangent)/epsilon
            numeric_curvature[:,i] += 0.5*(excenter_tangent - tangent)/epsilon

        #print str_matrix(abs(curvature - numeric_curvature))

        #print "analytic symmetric:", sum(numpy.ravel(curvature - numpy.transpose(curvature))**2)
        #print "numeric  symmetric:", sum(numpy.ravel(numeric_curvature - numpy.transpose(numeric_curvature))**2)
        #print "% 3.2e (% 3.2e, % 3.2e)" % (error, scale, nscale)

        #symmetry_error = sum(numpy.ravel(numeric_curvature - numpy.transpose(numeric_curvature))**2)
        #if symmetry_error > 1e-5:
        #    self.errors.append(
        #        "%2i/%2i %s Chain rule symmetry_error: % 3.2e" % (
        #            index, largest,
        #            internal_coordinate.label(),
        #            symmetry_error
        #        )
        #    )

        curvature_error = sum(numpy.ravel(curvature - numeric_curvature)**2)
        curvature_scale = sum(numpy.ravel(curvature)**2)
        ncurvature_scale = sum(numpy.ravel(numeric_curvature)**2)
        if curvature_error > 1e-6 and curvature_error * 1e6 > curvature_scale:
            self.errors.append(
                "%2i  %2i/%2i %s Chain rule curvature_error: % 3.2e (% 3.2e)\n%s" % (
                    icn, index, largest,
                    internal_coordinate.label(),
                    curvature_error, curvature_scale,
                    str_matrix(curvature - numeric_curvature)
                )
            )

    def test_random(self):
        internal_coordinates = self.load_internal_coordinates()

        def mutate_ethene(strength, ethene):
            result = copy.deepcopy(ethene)
            result.coordinates += numpy.random.uniform(-strength, strength, result.coordinates.shape)
            return result

        self.errors = []

        for icn, internal_coordinate in enumerate(internal_coordinates):
            for index in xrange(10):
                mod1_ethene = mutate_ethene(3, self.ethene)
                mod2_ethene = mutate_ethene(1e-5, mod1_ethene)
                self.tangent_test(icn, index, 99, internal_coordinate, mod1_ethene, mod2_ethene)
                self.curvature_test(icn, index, 99, internal_coordinate, mod1_ethene)

        self.assertEqual(len(self.errors), 0, "\n".join(self.errors))

