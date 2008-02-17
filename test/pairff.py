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

import molmod.pairff

import unittest, numpy, math

__all__ = ["PairFFTestCase", "CoulombFFTestCase"]


class Debug1FF(molmod.pairff.PairFF):
    def yield_pair_energies(self, index1, index2):
        yield self.distances[index1,index2]**2, 1

    def yield_pair_gradients(self, index1, index2):
        yield 2*self.distances[index1,index2], numpy.zeros(3)

    def yield_pair_hessians(self, index1, index2):
        yield 2, numpy.zeros((3,3))


class Debug2FF(molmod.pairff.PairFF):
    def yield_pair_energies(self, index1, index2):
        yield 1, sum((self.coordinates[index1] - self.coordinates[index2])**2)

    def yield_pair_gradients(self, index1, index2):
        yield 0, 2*(self.coordinates[index1] - self.coordinates[index2])

    def yield_pair_hessians(self, index1, index2):
        yield 0, 2*numpy.identity(3, float)


class Debug3FF(molmod.pairff.PairFF):
    def yield_pair_energies(self, index1, index2):
        yield self.distances[index1,index2]**3, sum((self.coordinates[index1] - self.coordinates[index2])**2)

    def yield_pair_gradients(self, index1, index2):
        yield 3*self.distances[index1,index2]**2, 2*(self.coordinates[index1] - self.coordinates[index2])

    def yield_pair_hessians(self, index1, index2):
        yield 6*self.distances[index1,index2], 2*numpy.identity(3, float)


class Debug4FF(molmod.pairff.PairFF):
    def yield_pair_energies(self, index1, index2):
        yield self.distances[index1,index2]**2, (self.coordinates[index1] - self.coordinates[index2])[0]**2

    def yield_pair_gradients(self, index1, index2):
        yield 2*self.distances[index1,index2], 2*(self.coordinates[index1] - self.coordinates[index2])*numpy.array([1, 0, 0])

    def yield_pair_hessians(self, index1, index2):
        yield 2, numpy.array([[2, 0, 0], [0, 0, 0], [0, 0, 0]], float)


class PairFFTestCase(unittest.TestCase):
    def make_coulombff(self, do_charges, do_dipoles, do_excludes):
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]], float
        )

        if do_charges:
            charges = numpy.array([0.3, 0.5, -0.8], float)
        else:
            charges = None

        if do_dipoles:
            dipoles = numpy.array([
                [ 0.2, -0.1, -0.7],
                [-0.8,  0.5,  0.3],
                [ 0.0, -0.6, -0.1]], float
            )
        else:
            dipoles = None

        if do_excludes:
            exclude_pairs = [set([0,2]), set([0,1])]
        else:
            exclude_pairs = []

        return molmod.pairff.CoulombFF(coordinates, charges=charges, dipoles=dipoles, exclude_pairs=exclude_pairs)

    def test_coulombff_c(self):
        self.ff_test(self.make_coulombff(do_charges=True,  do_dipoles=False, do_excludes=False))

    def test_coulombff_d(self):
        self.ff_test(self.make_coulombff(do_charges=False, do_dipoles=True,  do_excludes=False))

    def test_coulombff_cd(self):
        self.ff_test(self.make_coulombff(do_charges=True,  do_dipoles=True,  do_excludes=False))

    def test_coulombff_c_ex(self):
        self.ff_test(self.make_coulombff(do_charges=True,  do_dipoles=False, do_excludes=True))

    def test_coulombff_d_ex(self):
        self.ff_test(self.make_coulombff(do_charges=False, do_dipoles=True,  do_excludes=True))

    def test_coulombff_cd_ex(self):
        self.ff_test(self.make_coulombff(do_charges=True,  do_dipoles=True,  do_excludes=True))

    def make_dispersionff(self, do_excludes):
        atom_strengths = numpy.array([0.3, 0.5, 0.8], float)
        strengths = numpy.outer(atom_strengths, atom_strengths)
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]], float
        )

        if do_excludes:
            exclude_pairs = [set([0,2]), set([0,1])]
        else:
            exclude_pairs = []

        return molmod.pairff.DispersionFF(coordinates, strengths, exclude_pairs=exclude_pairs)

    def test_dispersionff(self):
        self.ff_test(self.make_dispersionff(do_excludes=False))

    def test_dispersionff_ex(self):
        self.ff_test(self.make_dispersionff(do_excludes=True))

    def make_debug1ff(self, do_excludes):
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]], float
        )

        if do_excludes:
            exclude_pairs = [set([0,2]), set([0,1])]
        else:
            exclude_pairs = []

        return Debug1FF(coordinates, exclude_pairs=exclude_pairs)

    def test_debug1ff(self):
        self.ff_test(self.make_debug1ff(do_excludes=False))

    def make_debug2ff(self, do_excludes):
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]], float
        )

        if do_excludes:
            exclude_pairs = [set([0,2]), set([0,1])]
        else:
            exclude_pairs = []

        return Debug2FF(coordinates, exclude_pairs=exclude_pairs)

    def test_debug2ff(self):
        self.ff_test(self.make_debug2ff(do_excludes=False))

    def make_debug3ff(self, do_excludes):
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]], float
        )

        if do_excludes:
            exclude_pairs = [set([0,2]), set([0,1])]
        else:
            exclude_pairs = []

        return Debug3FF(coordinates, exclude_pairs=exclude_pairs)

    def test_debug3ff(self):
        self.ff_test(self.make_debug3ff(do_excludes=False))

    def make_debug4ff(self, do_excludes):
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]], float
        )

        if do_excludes:
            exclude_pairs = [set([0,2]), set([0,1])]
        else:
            exclude_pairs = []

        return Debug4FF(coordinates, exclude_pairs=exclude_pairs)

    def test_debug4ff(self):
        self.ff_test(self.make_debug4ff(do_excludes=False))

    def ff_test(self, ff):
        coordinates = ff.coordinates
        numc = len(coordinates)

        energy = ff.energy()
        gradient = ff.gradient()
        hessian = ff.hessian()

        #print ff.hessian_flat()

        delta = 1e-5
        # test the yield_pair_gradient and yield_pair_hessian generators:
        for atom1 in xrange(len(coordinates)):
            for atom2 in xrange(atom1):
                ff.update_coordinates(coordinates)
                distance = ff.distances[atom1, atom2]
                an_se = 0.0
                an_ve = 0.0
                an_sg = 0.0
                an_vg = 0.0
                an_sh = 0.0
                an_vh = 0.0
                for (se, ve), (sg, vg), (sh, vh) in zip(
                    ff.yield_pair_energies(atom1, atom2),
                    ff.yield_pair_gradients(atom1, atom2),
                    ff.yield_pair_hessians(atom1, atom2)
                ):
                    an_se += se
                    an_ve += ve
                    an_sg += sg
                    an_vg += vg
                    an_sh += sh
                    an_vh += vh

                num_vg = numpy.zeros(an_vg.shape, float)
                num_vh = numpy.zeros(an_vh.shape, float)
                for i in xrange(3):
                    delta_coordinates = coordinates.copy()
                    delta_coordinates[atom1,i] += delta
                    ff.update_coordinates(delta_coordinates)
                    num_vg[i] = (sum(pair[1] for pair in ff.yield_pair_energies(atom1, atom2)) - an_ve)/delta
                    num_vh[i] = (sum(pair[1] for pair in ff.yield_pair_gradients(atom1, atom2)) - an_vg)/delta

                    num_sg = (sum(pair[0] for pair in ff.yield_pair_energies(atom1, atom2)) - an_se)/(ff.distances[atom1,atom2] - distance)
                    error = sum((num_sg - an_sg).ravel()**2)
                    reference = sum(num_sg.ravel()**2)
                    self.assertAlmostEqual(error, 0.0, 3, "num_sg: % 12.8f / % 12.8f" % (error, reference))
                    num_sh = (sum(pair[0] for pair in ff.yield_pair_gradients(atom1, atom2)) - an_sg)/(ff.distances[atom1,atom2] - distance)
                    error = sum((num_sh - an_sh).ravel()**2)
                    reference = sum(num_sh.ravel()**2)
                    self.assertAlmostEqual(error, 0.0, 3, "num_sh: % 12.8f / % 12.8f" % (error, reference))

                error = sum((num_vg - an_vg).ravel()**2)
                reference = sum(num_vg.ravel()**2)
                self.assertAlmostEqual(error, 0.0, 3, "num_vg: % 12.8f / % 12.8f" % (error, reference))

                error = sum((num_vh - an_vh).ravel()**2)
                reference = sum(num_vh.ravel()**2)
                self.assertAlmostEqual(error, 0.0, 3, "num_vh: % 12.8f / % 12.8f" % (error, reference))


        # 1) hessian should be symmetric
        hessian_flat = ff.hessian_flat()
        error = sum((hessian_flat - hessian_flat.transpose()).ravel()**2)
        reference = sum(hessian_flat.ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "1) The hessian is not symmetric: % 12.8f / % 12.8f" % (error, reference))

        # 1a) test the diagonal hessian blocks
        for atom in xrange(numc):
            error = sum((hessian[atom,atom] - hessian[atom,atom].transpose()).ravel()**2)
            reference = sum(hessian[atom,atom].ravel()**2)
            self.assertAlmostEqual(error, 0.0, 3, "1a) Diagonal hessian block %i is not symmetric: % 12.8f / % 12.8f" % (atom, error, reference))

        # 1b) test the off-diagonal hessian blocks
        for atom1 in xrange(numc):
            for atom2 in xrange(atom1):
                error = sum((hessian[atom1,atom2] - hessian[atom2,atom1].transpose()).ravel()**2)
                reference = sum(hessian[atom1,atom2].ravel()**2)
                self.assertAlmostEqual(error, 0.0, 3, "1a) Off-diagonal hessian block (%i,%i) is not symmetric: % 12.8f / % 12.8f" % (atom1, atom2, error, reference))

        # 2) test the cartesian gradient/hessian

        # 2a) test the analytical gradient
        numerical_gradient = numpy.zeros(gradient.shape, float)
        for atom in xrange(len(coordinates)):
            for index in xrange(3):
                delta_coordinates = coordinates.copy()
                delta_coordinates[atom,index] += delta
                ff.update_coordinates(delta_coordinates)
                numerical_gradient[atom,index] = (ff.energy() - energy) / delta
        error = sum((numerical_gradient - gradient).ravel()**2)
        reference = sum((numerical_gradient).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2a) The analytical gradient is incorrect: % 12.8f / % 12.8f" % (error, reference))

        # 2 pre_b) create a mask for the diagonal
        diagonal_mask = numpy.zeros(hessian.shape, float)
        for index in xrange(numc):
            diagonal_mask[index, index] = 1
        off_diagonal_mask = 1 - diagonal_mask

        # 2b) test the analytical hessian
        numerical_hessian = numpy.zeros(hessian.shape, float)
        for atom1 in xrange(len(coordinates)):
            for atom2 in xrange(len(coordinates)):
                for index1 in xrange(3):
                    for index2 in xrange(3):
                        delta_coordinates = coordinates.copy()
                        delta_coordinates[atom1,index1] += delta
                        delta_coordinates[atom2,index2] += delta
                        ff.update_coordinates(delta_coordinates)
                        numerical_hessian[atom1,atom2,index1,index2] = (ff.energy() - energy - delta*numerical_gradient[atom1,index1] - delta*numerical_gradient[atom2,index2])/(delta*delta)
        #print hessian
        #print ((numerical_hessian - hessian)*diagonal_mask)
        error = sum(((numerical_hessian - hessian)*diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))

        error = sum(((numerical_hessian - hessian)*off_diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*off_diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The off-diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))

        # 2c) test the analytical hessian in another way
        numerical_hessian = numpy.zeros(hessian.shape, float)
        for atom in xrange(len(coordinates)):
            for index in xrange(3):
                delta_coordinates = coordinates.copy()
                delta_coordinates[atom,index] += delta
                ff.update_coordinates(delta_coordinates)
                numerical_hessian[atom,:,index,:] = (ff.gradient() - gradient) / delta

        #print  numerical_hessian

        error = sum(((numerical_hessian - hessian)*diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))

        error = sum(((numerical_hessian - hessian)*off_diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*off_diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The off-diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))


class CoulombFFTestCase(unittest.TestCase):
    def test_cc1(self):
        coordinates = numpy.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = numpy.array([-1, 1], float)
        ff = molmod.pairff.CoulombFF(coordinates, charges)
        self.assertAlmostEqual(ff.energy(), -1.0, 5, "Incorrect energy.")

    def test_cc2(self):
        coordinates = numpy.array([
            [-1.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = numpy.array([1, -1, 1], float)
        ff = molmod.pairff.CoulombFF(coordinates, charges)
        self.assertAlmostEqual(ff.energy(), -1.5, 5, "Incorrect energy.")

    def test_cc3(self):
        coordinates = numpy.array([
            [ 0.0,  1.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = numpy.array([1, -1, 1], float)
        ff = molmod.pairff.CoulombFF(coordinates, charges)
        self.assertAlmostEqual(ff.energy(), -2.0 + 1/math.sqrt(2), 5, "Incorrect energy.")

    def test_cd1(self):
        coordinates = numpy.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = numpy.array([ 0, 1], float)
        dipoles = numpy.array([
            [ 1.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], float)
        ff = molmod.pairff.CoulombFF(coordinates, charges, dipoles)
        self.assertAlmostEqual(ff.energy(), 1.0, 5, "Incorrect energy.")


    def test_dd1(self):
        coordinates = numpy.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = numpy.array([ 0, 0], float)
        dipoles = numpy.array([
            [ 1.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        ff = molmod.pairff.CoulombFF(coordinates, charges, dipoles)
        self.assertAlmostEqual(ff.energy(), -2.0, 5, "Incorrect energy.")






