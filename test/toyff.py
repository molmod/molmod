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


from molmod.toyff import guess_geometry, ToyFF
from molmod.io.xyz import XYZFile
from molmod.io.sdf import SDFReader
from molmod.molecular_graphs import MolecularGraph

import unittest, numpy, os


__all__ = ["ToyFFTestCase"]


class ToyFFTestCase(unittest.TestCase):
    def load_molecule(self, xyz_fn):
        molecule = XYZFile(os.path.join("input", xyz_fn)).get_molecule()
        molecule.graph = MolecularGraph.from_geometry(molecule)
        return molecule

    def iter_molecules(self, allow_multi=False):
        xyz_fns = [
          "water.xyz", "cyclopentane.xyz", "ethene.xyz", "funny.xyz",
          "tea.xyz", "tpa.xyz", "thf_single.xyz", "precursor.xyz",
          "butane.xyz", "octane.xyz",
        ]
        for xyz_fn in xyz_fns:
            molecule = self.load_molecule(xyz_fn)
            if allow_multi or len(molecule.graph.independent_vertices) == 1:
                yield molecule
        sdf_fns = [
            "example.sdf", "CID_22898828.sdf", "SID_55127927.sdf",
            "SID_56274343.sdf", "SID_40363570.sdf", "SID_40363571.sdf",
            "SID_31646548.sdf", "SID_31646545.sdf", "SID_41893278.sdf",
            "SID_41893280.sdf", "SID_54258192.sdf", "SID_55488598.sdf",
        ]
        for sdf_fn in sdf_fns:
            for i, molecule in enumerate(SDFReader(os.path.join("input", sdf_fn))):
                if allow_multi or len(molecule.graph.independent_vertices) == 1:
                    yield molecule

    def test_guess_geometry(self):
        for input_mol in self.iter_molecules(allow_multi=False):
            output_mol = guess_geometry(input_mol.graph)
            output_mol.title = input_mol.title
            output_mol.write_to_file("output/guess_%s.xyz" % input_mol.title)

    def get_random_ff(self):
        N = 6
        edges = set([])
        while len(edges) < 2*N:
            v1 = numpy.random.randint(N)
            while True:
                v2 = numpy.random.randint(N)
                if v2 != v1:
                    break
            edges.add(frozenset([v1,v2]))
        edges = tuple(edges)
        numbers = numpy.random.randint(6, 10, N)
        graph = MolecularGraph(edges, numbers)
        ff = ToyFF(graph)

        from molmod.ext import molecules_distance_matrix
        while True:
            coordinates = numpy.random.uniform(0,3,(N,3))
            dm = molecules_distance_matrix(coordinates)
            if dm.max() > 1.0:
                break

        mask = numpy.zeros(dm.shape, bool)
        for i in xrange(N):
            for j in xrange(i):
                mask[i,j] = True

        return ff, coordinates, dm, mask

    def check_toyff_gradient(self, ff, coordinates):
        energy0, gradient0 = ff(coordinates, True)
        eps = numpy.random.uniform(-1e-6, 1e-6, coordinates.shape)
        energy1, gradient1 = ff(coordinates+eps, True)
        self.assertAlmostEqual(
            energy1 - energy0,
            0.5*numpy.dot(gradient0 + gradient1, eps.ravel())
        )

    def test_dm_quad_energy(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.dm_quad = 1.0
            energy = ff(coordinates, False)
            my_terms = (dm - ff.dm0)**2*ff.dmk
            my_terms[ff.dm0==0] = 0.0
            my_terms[mask] = 0.0
            self.assertAlmostEqual(energy, my_terms.sum())

    def test_dm_quad_gradient(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.dm_quad = 1.0
            self.check_toyff_gradient(ff, coordinates)

    def test_dm_reci_energy(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.dm_reci = 1.0
            energy = ff(coordinates, False)
            r0 = numpy.add.outer(ff.vdw_radii, ff.vdw_radii)
            d = dm/r0
            my_terms = (d-1)*(d-1)/d
            my_terms[ff.dm<=1] = 0.0
            my_terms[d>=1] = 0.0
            my_terms[mask] = 0.0
            self.assertAlmostEqual(energy, my_terms.sum())

    def test_dm_reci_gradient(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.dm_reci = 1.0
            self.check_toyff_gradient(ff, coordinates)

    def test_bond_quad_energy(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.bond_quad = 1.0
            energy = ff(coordinates, False)
            lengths = numpy.sqrt(((coordinates[ff.bond_edges[:,0]] - coordinates[ff.bond_edges[:,1]])**2).sum(axis=1))
            my_terms = (lengths - ff.bond_lengths)**2
            self.assertAlmostEqual(energy, my_terms.sum())

    def test_bond_quad_gradient(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.bond_quad = 1.0
            self.check_toyff_gradient(ff, coordinates)

    def test_bond_hyper_energy(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.bond_hyper = 1.0
            energy = ff(coordinates, False)
            lengths = numpy.sqrt(((coordinates[ff.bond_edges[:,0]] - coordinates[ff.bond_edges[:,1]])**2).sum(axis=1))
            my_terms = numpy.cosh((lengths - ff.bond_lengths)*ff.bond_hyper_scale)-1
            self.assertAlmostEqual(energy, my_terms.sum())

    def test_bond_hyper_gradient(self):
        for i in xrange(10):
            ff, coordinates, dm, mask = self.get_random_ff()
            ff.bond_hyper = 1.0
            self.check_toyff_gradient(ff, coordinates)


