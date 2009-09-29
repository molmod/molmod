# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from common import BaseTestCase

from molmod.molecules import *
from molmod.unit_cells import UnitCell
from molmod.io.xyz import XYZFile

import unittest, numpy


__all__ = ["MoleculesTestCase"]


class MoleculesTestCase(BaseTestCase):
    def test_distance_matrix(self):
        molecule = XYZFile("input/tpa.xyz").get_molecule()
        dm = 0
        for i in 0,1,2:
            dm += numpy.subtract.outer(molecule.coordinates[:,i],molecule.coordinates[:,i])**2
        dm = numpy.sqrt(dm)
        self.assert_((abs(molecule.distance_matrix - dm) < 1e-5).all(), "Wrong distance matrix")

    def test_distance_matrix_periodic(self):
        for i in xrange(10):
            N = 6
            unit_cell = UnitCell(
                numpy.random.uniform(0,1,(3,3)),
                numpy.random.randint(0,2,3).astype(bool),
            )
            fractional = numpy.random.uniform(0,1,(N,3))
            coordinates = unit_cell.to_cartesian(fractional)
            from molmod.ext import molecules_distance_matrix
            dm = molecules_distance_matrix(coordinates, unit_cell.matrix,
                                           unit_cell.reciprocal_zero)
            for i in xrange(N):
                for j in xrange(i,N):
                    delta = coordinates[j]-coordinates[i]
                    delta = unit_cell.shortest_vector(delta)
                    distance = numpy.linalg.norm(delta)
                    self.assertAlmostEqual(dm[i,j], distance)

    def test_read_only(self):
        numbers = [8, 1]
        coordinates = [
            [0, 1, 2],
            [-2, 3, 1],
        ]
        mol = Molecule(numbers, coordinates)
        try:
            mol.coordinates[0,1] = 5.0
            self.fail("Should have raised an error")
        except RuntimeError:
            pass
        self.assert_(isinstance(mol.coordinates[0,0], float))

        coordinates = numpy.array(coordinates)
        mol = Molecule(numbers, coordinates)
        coordinates[0,1] = 5.0
        self.assertAlmostEqual(mol.coordinates[0,1], 1.0)

    def test_mass(self):
        molecule = XYZFile("input/water.xyz").get_molecule()
        molecule.set_default_masses()
        self.assertAlmostEqual(molecule.mass, 32839.849030585152)

    def test_com(self):
        molecule = XYZFile("input/water.xyz").get_molecule()
        molecule.set_default_masses()
        expected_com = 0
        for i in xrange(3):
            expected_com += molecule.masses[i]*molecule.coordinates[i]
        expected_com /= molecule.mass
        self.assertArraysAlmostEqual(molecule.com, expected_com)

    def test_inertia_tensor(self):
        molecule = XYZFile("input/water.xyz").get_molecule()
        molecule.set_default_masses()
        expected_result = sum(
            m*(numpy.identity(3)*(r**2).sum()-numpy.outer(r,r))
            for m, r in zip(molecule.masses, (molecule.coordinates-molecule.com))
        )
        self.assertArraysAlmostEqual(molecule.inertia_tensor, expected_result)

    def test_chemical_formula(self):
        molecule = XYZFile("input/water.xyz").get_molecule()
        self.assertEqual(molecule.chemical_formula, "OH2")

    def test_copy(self):
        molecule = XYZFile("input/water.xyz").get_molecule()
        import copy
        tmp = copy.copy(molecule)
        self.assertEqual(id(tmp), id(molecule))
        tmp = copy.deepcopy(molecule)
        self.assertEqual(id(tmp), id(molecule))

