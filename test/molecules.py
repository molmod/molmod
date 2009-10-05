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

import unittest, numpy


__all__ = ["MoleculesTestCase"]


class MoleculesTestCase(BaseTestCase):
    def test_distance_matrix(self):
        molecule = Molecule.from_file("input/tpa.xyz")
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
        molecule = Molecule.from_file("input/water.xyz")
        molecule.set_default_masses()
        self.assertAlmostEqual(molecule.mass, 32839.849030585152)

    def test_com(self):
        molecule = Molecule.from_file("input/water.xyz")
        molecule.set_default_masses()
        expected_com = 0
        for i in xrange(3):
            expected_com += molecule.masses[i]*molecule.coordinates[i]
        expected_com /= molecule.mass
        self.assertArraysAlmostEqual(molecule.com, expected_com)

    def test_inertia_tensor(self):
        molecule = Molecule.from_file("input/water.xyz")
        molecule.set_default_masses()
        expected_result = sum(
            m*(numpy.identity(3)*(r**2).sum()-numpy.outer(r,r))
            for m, r in zip(molecule.masses, (molecule.coordinates-molecule.com))
        )
        self.assertArraysAlmostEqual(molecule.inertia_tensor, expected_result)

    def test_chemical_formula(self):
        molecule = Molecule.from_file("input/water.xyz")
        self.assertEqual(molecule.chemical_formula, "OH2")

    def test_copy(self):
        molecule = Molecule.from_file("input/water.xyz")
        import copy
        tmp = copy.copy(molecule)
        self.assertEqual(id(tmp), id(molecule))
        tmp = copy.deepcopy(molecule)
        self.assertEqual(id(tmp), id(molecule))

    def test_default_graph(self):
        molecule = Molecule.from_file("input/water.xyz")
        molecule.set_default_graph()
        self.assertEqual(molecule.graph.num_edges, 2)

    def test_from_file_cml(self):
        molecule = Molecule.from_file("input/caplayer.cml")
        self.assertEqual(molecule.size, 81)
        self.assertEqual(molecule.graph.num_edges, 90)

    def test_from_file_fchk(self):
        molecule = Molecule.from_file("input/1TOH.b3lyp.fchk")
        self.assertEqual(molecule.size, 9)

    def test_from_file_pdb(self):
        molecule = Molecule.from_file("input/il2.pdb")
        self.assertEqual(molecule.size, 2084)

    def test_from_file_sdf(self):
        molecule = Molecule.from_file("input/CID_22898828.sdf")
        self.assertEqual(molecule.size, 14)

    def test_from_file_xyz(self):
        molecule = Molecule.from_file("input/water.xyz")
        self.assertEqual(molecule.size, 3)

    def test_to_cml(self):
        mol0 = Molecule.from_file("input/caplayer.cml")
        mol0.write_to_file("output/caplayer.cml")
        mol1 = Molecule.from_file("output/caplayer.cml")
        self.assertEqual(mol0.title, mol1.title)
        self.assertArraysEqual(mol0.numbers, mol1.numbers)
        self.assertArraysAlmostEqual(mol0.coordinates, mol1.coordinates)
        self.assertEqual(mol0.graph.num_edges, mol1.graph.num_edges)

    def test_to_xyz(self):
        mol0 = Molecule.from_file("input/water.xyz")
        mol0.write_to_file("output/water.xyz")
        mol1 = Molecule.from_file("output/water.xyz")
        self.assertEqual(mol0.title, mol1.title)
        self.assertArraysEqual(mol0.numbers, mol1.numbers)
        self.assertArraysAlmostEqual(mol0.coordinates, mol1.coordinates)

    def test_copy_with(self):
        mol0 = Molecule.from_file("input/water.xyz")
        mol1 = mol0.copy_with(numbers=[3,4,5])
        self.assertArraysEqual(mol0.numbers, numpy.array([8,1,1]))
        self.assertArraysEqual(mol1.numbers, numpy.array([3,4,5]))
        self.assertArraysAlmostEqual(mol0.coordinates, mol1.coordinates)


