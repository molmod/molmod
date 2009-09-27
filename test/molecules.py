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


from molmod.molecules import *
from molmod.unit_cells import UnitCell
from molmod.io.xyz import XYZFile

import unittest, numpy


__all__ = ["MoleculesTestCase"]


class MoleculesTestCase(unittest.TestCase):
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
