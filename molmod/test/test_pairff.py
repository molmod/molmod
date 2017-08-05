# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
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


from __future__ import division

from builtins import range
import unittest

import numpy as np

from molmod.test.common import BaseTestCase
from molmod import *


__all__ = ["PairFFTestCase", "CoulombFFTestCase"]


class Debug1FF(PairFF):
    def yield_pair_energies(self, index1, index2):
        yield self.distances[index1,index2]**2, 1

    def yield_pair_gradients(self, index1, index2):
        yield 2*self.distances[index1,index2], np.zeros(3)

    def yield_pair_hessians(self, index1, index2):
        yield 2, np.zeros((3,3))


class Debug2FF(PairFF):
    def yield_pair_energies(self, index1, index2):
        yield 1, sum((self.coordinates[index1] - self.coordinates[index2])**2)

    def yield_pair_gradients(self, index1, index2):
        yield 0, 2*(self.coordinates[index1] - self.coordinates[index2])

    def yield_pair_hessians(self, index1, index2):
        yield 0, 2*np.identity(3, float)


class Debug3FF(PairFF):
    def yield_pair_energies(self, index1, index2):
        yield self.distances[index1,index2]**3, sum((self.coordinates[index1] - self.coordinates[index2])**2)

    def yield_pair_gradients(self, index1, index2):
        yield 3*self.distances[index1,index2]**2, 2*(self.coordinates[index1] - self.coordinates[index2])

    def yield_pair_hessians(self, index1, index2):
        yield 6*self.distances[index1,index2], 2*np.identity(3, float)


class Debug4FF(PairFF):
    def yield_pair_energies(self, index1, index2):
        yield self.distances[index1,index2]**2, (self.coordinates[index1] - self.coordinates[index2])[0]**2

    def yield_pair_gradients(self, index1, index2):
        yield 2*self.distances[index1,index2], 2*(self.coordinates[index1] - self.coordinates[index2])*np.array([1, 0, 0])

    def yield_pair_hessians(self, index1, index2):
        yield 2, np.array([[2, 0, 0], [0, 0, 0], [0, 0, 0]], float)


class PairFFTestCase(unittest.TestCase):
    def make_coulombff(self, do_charges, do_dipoles):
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        if do_charges:
            charges = np.array([0.3, 0.5, -0.8], float)
        else:
            charges = None

        if do_dipoles:
            dipoles = np.array([
                [ 0.2, -0.1, -0.7],
                [-0.8,  0.5,  0.3],
                [ 0.0, -0.6, -0.1]
            ], float)
        else:
            dipoles = None

        return CoulombFF(scaling, charges, dipoles, coordinates)

    def test_coulombff_c(self):
        self.check_ff(self.make_coulombff(do_charges=True,  do_dipoles=False))

    def test_coulombff_d(self):
        self.check_ff(self.make_coulombff(do_charges=False, do_dipoles=True))

    def test_coulombff_cd(self):
        self.check_ff(self.make_coulombff(do_charges=True,  do_dipoles=True))

    def make_dispersionff(self):
        atom_strengths = np.array([0.3, 0.5, 0.8], float)
        strengths = np.outer(atom_strengths, atom_strengths)
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return DispersionFF(scaling, strengths, coordinates)

    def test_dispersionff(self):
        self.check_ff(self.make_dispersionff())

    def make_pauliff(self):
        atom_strengths = np.array([0.3, 0.5, 0.8], float)
        strengths = np.outer(atom_strengths, atom_strengths)
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return PauliFF(scaling, strengths, coordinates)

    def test_pauliff(self):
        self.check_ff(self.make_pauliff())

    def make_exprepff(self):
        atom_As = np.array([0.3, 0.5, 0.8], float)
        As = np.sqrt(np.outer(atom_As, atom_As))
        atom_Bs = np.array([0.1, 0.2, 0.3], float)
        Bs = 0.5*np.add.outer(atom_Bs, atom_Bs)
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return ExpRepFF(scaling, As, Bs, coordinates)

    def test_exprepff(self):
        self.check_ff(self.make_exprepff())

    def make_debug1ff(self):
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return Debug1FF(scaling, coordinates)

    def test_debug1ff(self):
        self.check_ff(self.make_debug1ff())

    def make_debug2ff(self):
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return Debug2FF(scaling, coordinates)

    def test_debug2ff(self):
        self.check_ff(self.make_debug2ff())

    def make_debug3ff(self):
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return Debug3FF(scaling, coordinates)

    def test_debug3ff(self):
        self.check_ff(self.make_debug3ff())

    def make_debug4ff(self):
        coordinates = np.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]
        ], float)
        scaling = np.array([
            [0.0, 0.5, 0.9],
            [0.5, 0.0, 0.1],
            [0.9, 0.1, 0.0],
        ], float)

        return Debug4FF(scaling, coordinates)

    def test_debug4ff(self):
        self.check_ff(self.make_debug4ff())

    def check_ff(self, ff):
        coordinates = ff.coordinates
        numc = len(coordinates)

        energy = ff.energy()
        gradient = ff.gradient()
        hessian = ff.hessian()

        #print ff.hessian_flat()

        delta = 1e-5
        # test the yield_pair_gradient and yield_pair_hessian generators:
        for atom1 in range(len(coordinates)):
            for atom2 in range(atom1):
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

                num_vg = np.zeros(an_vg.shape, float)
                num_vh = np.zeros(an_vh.shape, float)
                for i in range(3):
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
        for atom in range(numc):
            error = sum((hessian[atom,:,atom,:] - hessian[atom,:,atom,:].transpose()).ravel()**2)
            reference = sum(hessian[atom,:,atom,:].ravel()**2)
            self.assertAlmostEqual(error, 0.0, 3, "1a) Diagonal hessian block %i is not symmetric: % 12.8f / % 12.8f" % (atom, error, reference))

        # 1b) test the off-diagonal hessian blocks
        for atom1 in range(numc):
            for atom2 in range(atom1):
                error = sum((hessian[atom1,:,atom2,:] - hessian[atom2,:,atom1,:].transpose()).ravel()**2)
                reference = sum(hessian[atom1,:,atom2,:].ravel()**2)
                self.assertAlmostEqual(error, 0.0, 3, "1a) Off-diagonal hessian block (%i,%i) is not symmetric: % 12.8f / % 12.8f" % (atom1, atom2, error, reference))

        # 2) test the cartesian gradient/hessian

        # 2a) test the analytical gradient
        numerical_gradient = np.zeros(gradient.shape, float)
        for atom in range(len(coordinates)):
            for index in range(3):
                delta_coordinates = coordinates.copy()
                delta_coordinates[atom,index] += delta
                ff.update_coordinates(delta_coordinates)
                numerical_gradient[atom,index] = (ff.energy() - energy) / delta
        error = sum((numerical_gradient - gradient).ravel()**2)
        reference = sum((numerical_gradient).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2a) The analytical gradient is incorrect: % 12.8f / % 12.8f" % (error, reference))

        # 2 pre_b) create a mask for the diagonal
        diagonal_mask = np.zeros(hessian.shape, float)
        for index in range(numc):
            diagonal_mask[index,:,index,:] = 1
        off_diagonal_mask = 1 - diagonal_mask

        # 2b) test the analytical hessian
        numerical_hessian = np.zeros(hessian.shape, float)
        for atom1 in range(len(coordinates)):
            for atom2 in range(len(coordinates)):
                for index1 in range(3):
                    for index2 in range(3):
                        delta_coordinates = coordinates.copy()
                        delta_coordinates[atom1,index1] += delta
                        delta_coordinates[atom2,index2] += delta
                        ff.update_coordinates(delta_coordinates)
                        numerical_hessian[atom1,index1,atom2,index2] = (ff.energy() - energy - delta*numerical_gradient[atom1,index1] - delta*numerical_gradient[atom2,index2])/(delta*delta)
        #print hessian
        #print ((numerical_hessian - hessian)*diagonal_mask)
        error = sum(((numerical_hessian - hessian)*diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))

        error = sum(((numerical_hessian - hessian)*off_diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*off_diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The off-diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))

        # 2c) test the analytical hessian in another way
        numerical_hessian = np.zeros(hessian.shape, float)
        for atom in range(len(coordinates)):
            for index in range(3):
                delta_coordinates = coordinates.copy()
                delta_coordinates[atom,index] += delta
                ff.update_coordinates(delta_coordinates)
                numerical_hessian[atom,index,:,:] = (ff.gradient() - gradient) / delta

        #print  numerical_hessian

        error = sum(((numerical_hessian - hessian)*diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))

        error = sum(((numerical_hessian - hessian)*off_diagonal_mask).ravel()**2)
        reference = sum((numerical_hessian*off_diagonal_mask).ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The off-diagonal blocks of the analytical hessian are incorrect: % 12.8f / %12.8f" % (error, reference))


class CoulombFFTestCase(BaseTestCase):
    def test_cc1(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        charges = np.array([-1, 1], float)
        ff = CoulombFF(scaling, charges, coordinates=coordinates)
        self.assertAlmostEqual(ff.energy(), -1.0, 5, "Incorrect energy.")

    def test_cc2(self):
        coordinates = np.array([
            [-1.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(3, float)
        charges = np.array([1, -1, 1], float)
        ff = CoulombFF(scaling, charges, coordinates=coordinates)
        self.assertAlmostEqual(ff.energy(), -1.5, 5, "Incorrect energy.")

    def test_cc3(self):
        coordinates = np.array([
            [ 0.0,  1.0,  0.0],
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(3, float)
        charges = np.array([1, -1, 1], float)
        ff = CoulombFF(scaling, charges, coordinates=coordinates)
        self.assertAlmostEqual(ff.energy(), -2.0 + 1/np.sqrt(2), 5, "Incorrect energy.")

    def test_cd1(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = np.array([ 0, 1], float)
        dipoles = np.array([
            [ 1.0,  0.0,  0.0],
            [ 0.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        ff = CoulombFF(scaling, charges, dipoles, coordinates)
        self.assertAlmostEqual(ff.energy(), 1.0, 5, "Incorrect energy.")

    def test_dd1(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = np.array([ 0, 0], float)
        dipoles = np.array([
            [ 1.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        ff = CoulombFF(scaling, charges, dipoles, coordinates=coordinates)
        self.assertAlmostEqual(ff.energy(), -2.0, 5, "Incorrect energy.")

    def test_esp_point_c(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        charges = np.array([-1, 1], float)
        ff = CoulombFF(scaling, charges, coordinates=coordinates)
        self.assertAlmostEqual(ff.esp_point(np.array([1.0, 0.0, 0.0])), 0.0, 5)
        self.assertAlmostEqual(ff.esp_point(np.array([0.0, 2.0, 0.0])), -0.5+0.125**0.5, 5)

    def test_esp_point_d(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        dipoles = np.array([
            [ 1.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        ff = CoulombFF(scaling, charges=None, dipoles=dipoles, coordinates=coordinates)
        self.assertAlmostEqual(ff.esp_point(np.array([0.5, 0.0, 0.0])), 0.0, 5)
        self.assertAlmostEqual(ff.esp_point(np.array([0.0, 1.0, 0.0])), -0.5*0.5**0.5, 5)

    def test_esp_c(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        charges = np.array([-1, 1], float)
        ff = CoulombFF(scaling, charges, coordinates=coordinates)
        self.assertArraysAlmostEqual(ff.esp(), np.array([0.5, -0.5]), 1e-5)

    def test_esp_d(self):
        coordinates = np.array([
            [ 0.0,  0.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        dipoles = np.array([
            [ 0.0,  1.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        ff = CoulombFF(scaling, charges=None, dipoles=dipoles, coordinates=coordinates)
        self.assertArraysAlmostEqual(ff.esp(), np.array([-1.0, 0.0]), 1e-5)

    def test_efield_point(self):
        coordinates = np.array([
            [-2.0,  0.0,  0.0],
            [ 2.0,  0.0,  0.0]
        ], float)
        scaling = 1 - np.identity(2, float)
        dipoles = np.array([
            [ 0.0,  1.0,  0.0],
            [ 1.0,  0.0,  0.0]
        ], float)
        charges = np.array([-1, 1], float)
        ff = CoulombFF(scaling, charges=charges, dipoles=dipoles, coordinates=coordinates)
        eps = 1e-5
        for i in range(10):
            center = np.random.uniform(-0.2, 0.2, 3)
            esp0 = ff.esp_point(center)
            efield0 = ff.efield_point(center)
            for j in range(3):
                ex_center = center.copy()
                ex_center[j] += eps
                self.assertAlmostEqual(-efield0[j], (ff.esp_point(ex_center) - esp0)/eps, 3)

    def test_efield(self):
        coordinates = np.random.uniform(-1, 1, (3,3))
        point = coordinates[0]
        scaling = 1 - np.identity(3, float)
        dipoles = np.random.uniform(-1, 1, (3,3))
        charges = np.random.uniform(-1, 1, 3)
        dipoles[0] = 0
        charges[0] = 1
        ff1 = CoulombFF(scaling, charges=charges, dipoles=dipoles, coordinates=coordinates)
        coordinates = coordinates[[1,2]]
        scaling = 1 - np.identity(2, float)
        dipoles = dipoles[[1,2]]
        charges = charges[[1,2]]
        ff2 = CoulombFF(scaling, charges=charges, dipoles=dipoles, coordinates=coordinates)
        self.assertArraysAlmostEqual(ff1.gradient()[0], -ff1.efield()[0])
        self.assertArraysAlmostEqual(ff1.gradient()[0], -ff2.efield_point(point))
