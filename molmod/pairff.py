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
"""Reference implementations of a few conventional nonbonding force fields

   The main purpose of this implementation is reliability, not speed. These
   routines can be used to validate an efficient low level implementation.
"""


from __future__ import division

from builtins import range
import numpy as np

__all__ = [
    "PairFF", "CoulombFF", "DispersionFF", "PauliFF", "ExpRepFF",
]


class PairFF(object):
    """Evaluates the energy, gradient and Hessian of pairwise potential

       The expressions is of the energetic model has the form
       E = sum_{i!=j} sum_k s_k(r_ij)*v_(bar{r}_ij).
       In the derived classes one must provide functions that iterate over all
       the corresponding function values, derivatives and second derivatives of
       s and v for a given r_ij.
    """

    def __init__(self, scaling, coordinates=None):
        """Initialize a pair potential object

           Arguments:
             scaling  --  symmetric NxN array with pairwise scaling factors.
                          When an element is set to zero, it will be excluded.

           Optional argument:
             coordinates  --  the initial Cartesian coordinates of the system,
                              which can be updated with the update_coordinates
                              method
        """
        if coordinates is not None:
            self.update_coordinates(coordinates)
        self.scaling = scaling
        self.scaling.ravel()[::len(self.scaling)+1] = 0

    def update_coordinates(self, coordinates=None):
        """Update the coordinates (and derived quantities)

           Argument:
             coordinates  --  new Cartesian coordinates of the system
        """
        if coordinates is not None:
            self.coordinates = coordinates
        self.numc = len(self.coordinates)
        self.distances = np.zeros((self.numc, self.numc), float)
        self.deltas = np.zeros((self.numc, self.numc, 3), float)
        self.directions = np.zeros((self.numc, self.numc, 3), float)
        self.dirouters = np.zeros((self.numc, self.numc, 3, 3), float)
        for index1, coordinate1 in enumerate(self.coordinates):
            for index2, coordinate2 in enumerate(self.coordinates):
                delta = coordinate1 - coordinate2
                self.deltas[index1, index2] = delta
                distance = np.linalg.norm(delta)
                self.distances[index1, index2] = distance
                if index1 != index2:
                    tmp = delta/distance
                    self.directions[index1, index2] = tmp
                    self.dirouters[index1, index2] = np.outer(tmp, tmp)

    def yield_pair_energies(self, index1, index2):
        """Yields pairs ((s(r_ij), v(bar{r}_ij))"""
        raise NotImplementedError

    def yield_pair_gradients(self, index1, index2):
        """Yields pairs ((s'(r_ij), grad_i v(bar{r}_ij))"""
        raise NotImplementedError

    def yield_pair_hessians(self, index1, index2):
        """Yields pairs ((s''(r_ij), grad_i (x) grad_i v(bar{r}_ij))"""
        raise NotImplementedError

    def energy(self):
        """Compute the energy of the system"""
        result = 0.0
        for index1 in range(self.numc):
            for index2 in range(index1):
                if self.scaling[index1, index2] > 0:
                    for se, ve in self.yield_pair_energies(index1, index2):
                        result += se*ve*self.scaling[index1, index2]
        return result

    def gradient_component(self, index1):
        """Compute the gradient of the energy for one atom"""
        result = np.zeros(3, float)
        for index2 in range(self.numc):
            if self.scaling[index1, index2] > 0:
                for (se, ve), (sg, vg) in zip(self.yield_pair_energies(index1, index2), self.yield_pair_gradients(index1, index2)):
                    result += (sg*self.directions[index1, index2]*ve + se*vg)*self.scaling[index1, index2]
        return result

    def gradient(self):
        """Compute the gradient of the energy for all atoms"""
        result = np.zeros((self.numc, 3), float)
        for index1 in range(self.numc):
            result[index1] = self.gradient_component(index1)
        return result

    def hessian_component(self, index1, index2):
        """Compute the hessian of the energy for one atom pair"""
        result = np.zeros((3, 3), float)
        if index1 == index2:
            for index3 in range(self.numc):
                if self.scaling[index1, index3] > 0:
                    d_1 = 1/self.distances[index1, index3]
                    for (se, ve), (sg, vg), (sh, vh) in zip(
                        self.yield_pair_energies(index1, index3),
                        self.yield_pair_gradients(index1, index3),
                        self.yield_pair_hessians(index1, index3)
                    ):
                        result += (
                            +sh*self.dirouters[index1, index3]*ve
                            +sg*(np.identity(3, float) - self.dirouters[index1, index3])*ve*d_1
                            +sg*np.outer(self.directions[index1, index3],  vg)
                            +sg*np.outer(vg, self.directions[index1, index3])
                            +se*vh
                        )*self.scaling[index1, index3]
        elif self.scaling[index1, index2] > 0:
            d_1 = 1/self.distances[index1, index2]
            for (se, ve), (sg, vg), (sh, vh) in zip(
                self.yield_pair_energies(index1, index2),
                self.yield_pair_gradients(index1, index2),
                self.yield_pair_hessians(index1, index2)
            ):
                result -= (
                    +sh*self.dirouters[index1, index2]*ve
                    +sg*(np.identity(3, float) - self.dirouters[index1, index2])*ve*d_1
                    +sg*np.outer(self.directions[index1, index2],  vg)
                    +sg*np.outer(vg, self.directions[index1, index2])
                    +se*vh
                )*self.scaling[index1, index2]
        return result

    def hessian(self):
        """Compute the hessian of the energy"""
        result = np.zeros((self.numc, 3, self.numc, 3), float)
        for index1 in range(self.numc):
            for index2 in range(self.numc):
                result[index1, :, index2, :] = self.hessian_component(index1, index2)
        return result

    def gradient_flat(self):
        """Return the gradient a 3N array"""
        return self.gradient().ravel()

    def hessian_flat(self):
        """Return the hessian a 3N x 3N array"""
        return self.hessian().reshape((self.numc*3, self.numc*3))


class CoulombFF(PairFF):
    """Computes the electrostatic interactions using charges and point dipoles"""

    def __init__(self, scaling, charges=None, dipoles=None, coordinates=None):
        """Initialize a CoulombFF object

           Arguments:
             scaling  --  symmetric NxN array with pairwise scaling factors.
                          When an element is set to zero, it will be excluded.

           Optio1nal arguments:
             charges  --  the atomic partial charges
             dipoles  --  atomic dipole moments
             coordinates  --  the initial Cartesian coordinates of the system,
                              which can be updated with the update_coordinates
                              method
        """
        PairFF.__init__(self, scaling, coordinates)
        self.charges = charges
        self.dipoles = dipoles

    def yield_pair_energies(self, index1, index2):
        """Yields pairs ((s(r_ij), v(bar{r}_ij))"""
        d_1 = 1/self.distances[index1, index2]
        if self.charges is not None:
            c1 = self.charges[index1]
            c2 = self.charges[index2]
            yield c1*c2*d_1, 1
        if self.dipoles is not None:
            d_3 = d_1**3
            d_5 = d_1**5
            delta = self.deltas[index1, index2]
            p1 = self.dipoles[index1]
            p2 = self.dipoles[index2]
            yield d_3*np.dot(p1, p2), 1
            yield -3*d_5, np.dot(p1, delta)*np.dot(delta, p2)
            if self.charges is not None:
                yield c1*d_3, np.dot(p2, delta)
                yield c2*d_3, np.dot(p1, -delta)

    def yield_pair_gradients(self, index1, index2):
        """Yields pairs ((s'(r_ij), grad_i v(bar{r}_ij))"""
        d_2 = 1/self.distances[index1, index2]**2
        if self.charges is not None:
            c1 = self.charges[index1]
            c2 = self.charges[index2]
            yield -c1*c2*d_2, np.zeros(3)
        if self.dipoles is not None:
            d_4 = d_2**2
            d_6 = d_2**3
            delta = self.deltas[index1, index2]
            p1 = self.dipoles[index1]
            p2 = self.dipoles[index2]
            yield -3*d_4*np.dot(p1, p2), np.zeros(3)
            yield 15*d_6, p1*np.dot(p2, delta) + p2*np.dot(p1, delta)
            if self.charges is not None:
                yield -3*c1*d_4, p2
                yield -3*c2*d_4, -p1

    def yield_pair_hessians(self, index1, index2):
        """Yields pairs ((s''(r_ij), grad_i (x) grad_i v(bar{r}_ij))"""
        d_1 = 1/self.distances[index1, index2]
        d_3 = d_1**3
        if self.charges is not None:
            c1 = self.charges[index1]
            c2 = self.charges[index2]
            yield 2*c1*c2*d_3, np.zeros((3, 3))
        if self.dipoles is not None:
            d_5 = d_1**5
            d_7 = d_1**7
            p1 = self.dipoles[index1]
            p2 = self.dipoles[index2]
            yield 12*d_5*np.dot(p1, p2), np.zeros((3, 3))
            yield -90*d_7, np.outer(p1, p2) + np.outer(p2, p1)
            if self.charges is not None:
                yield 12*c1*d_5, np.zeros((3, 3))
                yield 12*c2*d_5, np.zeros((3, 3))

    def esp_point(self, point):
        result = 0.0
        for index2 in range(self.numc):
            delta = point - self.coordinates[index2]
            d = np.linalg.norm(delta)
            if self.charges is not None:
                result += self.charges[index2]/d
            if self.dipoles is not None:
                direction = delta/d
                result += np.dot(self.dipoles[index2], direction)/d**2
        return result

    def esp_component(self, index1):
        result = 0.0
        for index2 in range(self.numc):
            if self.scaling[index1, index2] > 0:
                d = self.distances[index1, index2]
                if self.charges is not None:
                    result += self.charges[index2]/d
                if self.dipoles is not None:
                    result += np.dot(self.dipoles[index2], self.directions[index1, index2])/d**2
        return result

    def esp(self):
        """Compute the electrostatic potential at each atom due to other atoms"""
        result = np.zeros(self.numc, float)
        for index1 in range(self.numc):
            result[index1] = self.esp_component(index1)
        return result

    def efield_point(self, point):
        result = 0.0
        for index2 in range(self.numc):
            delta = point - self.coordinates[index2]
            d = np.linalg.norm(delta)
            direction = delta/d
            if self.charges is not None:
                result += self.charges[index2]*direction/d**2
            if self.dipoles is not None:
                p = self.dipoles[index2]
                result += (3*np.dot(p, direction)*direction - p)/d**3
        return result

    def efield_component(self, index1):
        result = 0.0
        for index2 in range(self.numc):
            if self.scaling[index1, index2] > 0:
                d = self.distances[index1, index2]
                direction = self.directions[index1, index2]
                if self.charges is not None:
                    result += self.charges[index2]*direction/d**2
                if self.dipoles is not None:
                    p = self.dipoles[index2]
                    result += (3*np.dot(p, direction)*direction - p)/d**3
        return result

    def efield(self):
        """Compute the electrostatic potential at each atom due to other atoms"""
        result = np.zeros((self.numc,3), float)
        for index1 in range(self.numc):
            result[index1] = self.efield_component(index1)
        return result



class DispersionFF(PairFF):
    """Computes the London dispersion interaction"""

    def __init__(self, scaling, strengths, coordinates=None):
        """Initialize a DispersionFF object

           Arguments:
             scaling  --  symmetric NxN array with pairwise scaling factors.
                          When an element is set to zero, it will be excluded.
             strengths  --  a symmetric with linear coefficients in front of
                            r**-6 for each atom pair

           Optional arguments:
             coordinates  --  the initial Cartesian coordinates of the system,
                              which can be updated with the update_coordinates
                              method
        """
        PairFF.__init__(self, scaling, coordinates)
        self.strengths = strengths

    def yield_pair_energies(self, index1, index2):
        """Yields pairs ((s(r_ij), v(bar{r}_ij))"""
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield strength*distance**(-6), 1

    def yield_pair_gradients(self, index1, index2):
        """Yields pairs ((s'(r_ij), grad_i v(bar{r}_ij))"""
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield -6*strength*distance**(-7), np.zeros(3)

    def yield_pair_hessians(self, index1, index2):
        """Yields pairs ((s''(r_ij), grad_i (x) grad_i v(bar{r}_ij))"""
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield 42*strength*distance**(-8), np.zeros((3, 3))


class PauliFF(PairFF):
    """Computes the Pauli repulsion interaction"""

    def __init__(self, scaling, strengths, coordinates=None):
        """Initialize a PauliFF

           Arguments:
             scaling  --  symmetric NxN array with pairwise scaling factors.
                          When an element is set to zero, it will be excluded.
             strengths  --  a symmetric with linear coefficients in front of
                            r**-12 for each atom pair

           Optional arguments:
             coordinates  --  the initial Cartesian coordinates of the system,
                              which can be updated with the update_coordinates
                              method
        """
        PairFF.__init__(self, scaling, coordinates)
        self.strengths = strengths

    def yield_pair_energies(self, index1, index2):
        """Yields pairs ((s(r_ij), v(bar{r}_ij))"""
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield strength*distance**(-12), 1

    def yield_pair_gradients(self, index1, index2):
        """Yields pairs ((s'(r_ij), grad_i v(bar{r}_ij))"""
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield -12*strength*distance**(-13), np.zeros(3)

    def yield_pair_hessians(self, index1, index2):
        """Yields pairs ((s''(r_ij), grad_i (x) grad_i v(bar{r}_ij))"""
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield 12*13*strength*distance**(-14), np.zeros((3, 3))


class ExpRepFF(PairFF):
    """Computes the exponential repulsion interaction"""

    def __init__(self, scaling, As, Bs, coordinates=None):
        """Initialize a ExpRepFF

           Arguments:
             scaling  --  symmetric NxN array with pairwise scaling factors.
                          When an element is set to zero, it will be excluded.
             As  --  A matrix with pre-exponential factors
             Bs  --  A matrix with exponents

           The repulsion has the form A*exp(-B*r)

           Optional arguments:
             coordinates  --  the initial Cartesian coordinates of the system,
                              which can be updated with the update_coordinates
                              method
        """
        PairFF.__init__(self, scaling, coordinates)
        self.As = As
        self.Bs = Bs

    def yield_pair_energies(self, index1, index2):
        """Yields pairs ((s(r_ij), v(bar{r}_ij))"""
        A = self.As[index1, index2]
        B = self.Bs[index1, index2]
        distance = self.distances[index1, index2]
        yield A*np.exp(-B*distance), 1

    def yield_pair_gradients(self, index1, index2):
        """Yields pairs ((s'(r_ij), grad_i v(bar{r}_ij))"""
        A = self.As[index1, index2]
        B = self.Bs[index1, index2]
        distance = self.distances[index1, index2]
        yield -B*A*np.exp(-B*distance), np.zeros(3)

    def yield_pair_hessians(self, index1, index2):
        """Yields pairs ((s''(r_ij), grad_i (x) grad_i v(bar{r}_ij))"""
        A = self.As[index1, index2]
        B = self.Bs[index1, index2]
        distance = self.distances[index1, index2]
        yield B*B*A*np.exp(-B*distance), np.zeros((3, 3))
