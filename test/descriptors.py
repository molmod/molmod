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


from molmod.data.periodic import periodic
from molmod.descriptors import MolecularDescriptorTV1, sph_harmonics
from molmod.transformations import Complete

from molmod.io.xyz import XYZFile

import numpy

import unittest, math, random


__all__ = ["DescriptorTestCase"]


class DescriptorTestCase(unittest.TestCase):
    def load_molecule(self, filename):
        m = XYZFile("input/"+filename).get_molecule()
        return m.coordinates, numpy.array([
            periodic[number].mass
            for number in m.numbers
        ], float)

    def do_test(self, coordinates, masses, expected_is):
        reference = MolecularDescriptorTV1(coordinates, masses)
        self.assert_(reference.inversion_symmetric == expected_is)

        for counter in xrange(3):
            transformation = Complete()
            transformation.set_rotation_properties(
                random.uniform(0, 2*math.pi),
                numpy.random.uniform(-1, 1, (3,)),
                False,
            )
            transformation.t = numpy.random.uniform(-5, 5, (3,))
            new_coordinates = numpy.array([
                transformation.vector_apply(coordinate)
                for coordinate
                in coordinates
            ], float)
            new_descriptor = MolecularDescriptorTV1(new_coordinates, masses)
            self.assert_(reference.compare_structure(new_descriptor))
            self.assert_(not reference.compare_global_rotation(new_descriptor))
            self.assert_(not reference.compare_global_translation(new_descriptor))

    def test_precursor(self):
        coordinates, masses = self.load_molecule("precursor.xyz")
        self.do_test(coordinates, masses, expected_is=False)

    def test_water(self):
        coordinates, masses = self.load_molecule("water.xyz")
        self.do_test(coordinates, masses, expected_is=False)

    def test_tpa(self):
        coordinates, masses = self.load_molecule("tpa.xyz")
        self.do_test(coordinates, masses, expected_is=False)

    def test_mfi_cluster(self):
        coordinates, masses = self.load_molecule("mfi_cluster.xyz")
        self.do_test(coordinates, masses, expected_is=True)

    def test_sph_harmonics(self):
        def Yl0m00(phi, theta):
            return 0.5*math.sqrt(1/math.pi)


        def Yl1m1n(phi, theta):
            return 0.5*math.sqrt(1.5/math.pi)*numpy.exp(-1j*phi)*math.sin(theta)

        def Yl1m00(phi, theta):
            return 0.5*math.sqrt(3.0/math.pi)*math.cos(theta)

        def Yl1m1p(phi, theta):
            return -0.5*math.sqrt(1.5/math.pi)*numpy.exp(+1j*phi)*math.sin(theta)


        def Yl2m2n(phi, theta):
            return 0.25*math.sqrt(7.5/math.pi)*numpy.exp(-2j*phi)*math.sin(theta)**2

        def Yl2m1n(phi, theta):
            return 0.5*math.sqrt(7.5/math.pi)*numpy.exp(-1j*phi)*math.sin(theta)*math.cos(theta)

        def Yl2m00(phi, theta):
            return 0.25*math.sqrt(5.0/math.pi)*(3*math.cos(theta)**2-1)

        def Yl2m1p(phi, theta):
            return -0.5*math.sqrt(7.5/math.pi)*numpy.exp(+1j*phi)*math.sin(theta)*math.cos(theta)

        def Yl2m2p(phi, theta):
            return 0.25*math.sqrt(7.5/math.pi)*numpy.exp(+2j*phi)*math.sin(theta)**2

        phi = 0.1
        theta = 0.4

        values = sph_harmonics(2, phi, theta)
        test_values = [
            numpy.array([Yl0m00(phi, theta)]),
            numpy.array([Yl1m1n(phi, theta), Yl1m00(phi, theta), Yl1m1p(phi, theta)]),
            numpy.array([Yl2m2n(phi, theta), Yl2m1n(phi, theta), Yl2m00(phi, theta), Yl2m1p(phi, theta), Yl2m2p(phi, theta)]),
        ]

        self.assertAlmostEqual(abs(values[0][0] - test_values[0][0]), 0, 5, "Yl0m00 is wrong.")
        self.assertAlmostEqual(abs(values[1][0] - test_values[1][0]), 0, 5, "Yl1m1n is wrong.")
        self.assertAlmostEqual(abs(values[1][1] - test_values[1][1]), 0, 5, "Yl1m00 is wrong.")
        self.assertAlmostEqual(abs(values[1][2] - test_values[1][2]), 0, 5, "Yl1m1p is wrong.")
        self.assertAlmostEqual(abs(values[2][0] - test_values[2][0]), 0, 5, "Yl2m2n is wrong.")
        self.assertAlmostEqual(abs(values[2][1] - test_values[2][1]), 0, 5, "Yl2m1n is wrong.")
        self.assertAlmostEqual(abs(values[2][2] - test_values[2][2]), 0, 5, "Yl2m00 is wrong.")
        self.assertAlmostEqual(abs(values[2][3] - test_values[2][3]), 0, 5, "Yl2m1p is wrong.")
        self.assertAlmostEqual(abs(values[2][4] - test_values[2][4]), 0, 5, "Yl2m2p is wrong.")



