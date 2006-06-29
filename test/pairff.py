# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of MolMod.
# 
# MolMod is free software; you can redistribute it and/or
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

from molmod.pairff import CoulombFF

import unittest, numpy

__all__ = ["PairFF"]


class PairFF(unittest.TestCase):
    def test_coulombff(self):
        charges = numpy.array([0.3, 0.5, -0.8], float)
        coordinates = numpy.array([
            [ 0.5, 2.5, 0.1],
            [-1.2, 0.4, 0.3],
            [ 0.3, 0.9, 0.7]],
            float
        )
        numc = len(coordinates)

        ff = CoulombFF(coordinates, charges)
        energy = ff.energy()
        gradient = ff.gradient()
        hessian = ff.hessian()

        # 0) hessian should be symmetric
        hessian_flat = ff.hessian_flat()
        error = sum((hessian_flat - hessian_flat.transpose()).ravel()**2)
        reference = sum(hessian_flat.ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "0) The hessian is not symmetric: % 12.8f / % 12.8f" % (error, reference))

        delta = 1e-5

        # 1) test the individual pairs

        for index1 in xrange(numc):
            for index2 in xrange(numc):
                if index1 == index2:
                    continue
                # 1a) test the pair_gradient
                e1 = ff.pair_energy(ff.distances[index1,index2], index1, index2)
                e2 = ff.pair_energy(ff.distances[index1,index2]+delta, index1, index2)
                g1 = ff.pair_gradient(ff.distances[index1,index2], index1, index2)
                error = (g1 - (e2-e1)/delta)**2
                reference = ((e2-e1)/delta)**2
                self.assertAlmostEqual(error, 0.0, 3, "1a) The pair gradient is wrong: % 12.8f / % 12.8f" % (error, reference))
                
                # 1b) test the pair_hessian
                g2 = ff.pair_gradient(ff.distances[index1,index2]+delta, index1, index2)
                h1 = ff.pair_hessian(ff.distances[index1,index2], index1, index2)
                error = (h1 - (g2-g1)/delta)**2
                reference = ((g2-g1)/delta)**2
                self.assertAlmostEqual(error, 0.0, 3, "1b) The pair hessian is wrong: % 12.8f / % 12.8f" % (error, reference))
        
        
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
        error = sum((numerical_hessian - hessian).ravel()**2)
        reference = sum(numerical_hessian.ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2b) The analytical hessian is incorrect: % 12.8f / %12.8f" % (error, reference))
        
        # 2c) test the analytical hessian in another way
        numerical_hessian = numpy.zeros(hessian.shape, float)
        for atom in xrange(len(coordinates)):
            for index in xrange(3):
                delta_coordinates = coordinates.copy()
                delta_coordinates[atom,index] += delta
                ff.update_coordinates(delta_coordinates)
                numerical_hessian[atom,:,index,:] = (ff.gradient() - gradient) / delta
        error = sum((numerical_hessian - hessian).ravel()**2)
        reference = sum(numerical_hessian.ravel()**2)
        self.assertAlmostEqual(error, 0.0, 3, "2c) The analytical hessian is incorrect: % 12.8f / %12.8f" % (error, reference))

