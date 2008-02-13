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


import math, numpy


class PairFF(object):
    """
    Evaluates the energy, gradient and hessian of a energetic model consisting
    of pair interactions between point particles.

    The expressions is of the energetic model has the form
    E = sum_{i!=j} sum_k s_k(r_ij)*v_(bar{r}_ij).
    In the derived classes one must provide functions that yield all the
    corresponding function values, derivatives and second derivatives of s and
    v for a given r_ij.

    The class is initialized with the coordinates of the points, pairs of
    points whos interactions have to be excluded, and userdefined data that
    specifies the interactions in the derived classes.
    """
    def __init__(self, coordinates, exclude_pairs=[]):
        self.update_coordinates(coordinates)
        self.exclude_pairs = exclude_pairs

    def update_coordinates(self, coordinates=None):
        if coordinates is not None:
            self.coordinates = coordinates
        self.numc = len(self.coordinates)
        self.distances = numpy.zeros((self.numc, self.numc), float)
        self.deltas = numpy.zeros((self.numc, self.numc, 3), float)
        self.directions = numpy.zeros((self.numc, self.numc, 3), float)
        self.dirouters = numpy.zeros((self.numc, self.numc, 3, 3), float)
        for index1, coordinate1 in enumerate(self.coordinates):
            for index2, coordinate2 in enumerate(self.coordinates):
                delta = coordinate1 - coordinate2
                self.deltas[index1,index2] = delta
                distance = math.sqrt(numpy.dot(delta, delta))
                self.distances[index1, index2] = distance
                if index1 != index2:
                    tmp = delta/distance
                    self.directions[index1, index2] = tmp
                    self.dirouters[index1, index2] = numpy.outer(tmp, tmp)

    def yield_pair_energies(self, index1, index2):
        """
        Yields pairs ((s(r_ij), v(bar{r}_ij)). (implemented in derived classes)
        """
        raise NotImplementedError

    def yield_pair_gradients(self, index1, index2):
        """
        Yields pairs ((s'(r_ij), grad_i v(bar{r}_ij)). (implemented in derived classes)
        """
        raise NotImplementedError

    def yield_pair_hessians(self, index1, index2):
        """
        Yields pairs ((s''(r_ij), grad_i (x) grad_i v(bar{r}_ij)). (implemented in derived classes)
        """
        raise NotImplementedError

    def energy(self):
        result = 0.0
        for index1 in xrange(self.numc):
            for index2 in xrange(index1):
                if not frozenset([index1, index2]) in self.exclude_pairs:
                    for se, ve in self.yield_pair_energies(index1, index2):
                        result += se*ve
        return result

    def gradient_component(self, index1):
        result = numpy.zeros(3, float)
        for index2 in xrange(self.numc):
            if index2 != index1 and not frozenset([index1, index2]) in self.exclude_pairs:
                for (se, ve), (sg, vg) in zip(self.yield_pair_energies(index1, index2), self.yield_pair_gradients(index1, index2)):
                    result += sg*self.directions[index1, index2]*ve + se*vg
        return result

    def gradient(self):
        result = numpy.zeros((self.numc, 3), float)
        for index1 in xrange(self.numc):
            result[index1] = self.gradient_component(index1)
        return result

    def hessian_component(self, index1, index2):
        result = numpy.zeros((3, 3), float)
        if index1 == index2:
            for index3 in xrange(self.numc):
                if index3 != index1 and not frozenset([index1, index3]) in self.exclude_pairs:
                    d_1 = 1/self.distances[index1,index3]
                    for (se, ve), (sg, vg), (sh, vh) in zip(
                        self.yield_pair_energies(index1, index3),
                        self.yield_pair_gradients(index1, index3),
                        self.yield_pair_hessians(index1, index3)
                    ):
                        result += (
                            +sh*self.dirouters[index1,index3]*ve
                            +sg*(numpy.identity(3, float) - self.dirouters[index1, index3])*ve*d_1
                            +sg*numpy.outer(self.directions[index1, index3],  vg)
                            +sg*numpy.outer(vg, self.directions[index1, index3])
                            +se*vh
                        )
        elif not frozenset([index1, index2]) in self.exclude_pairs:
            d_1 = 1/self.distances[index1,index2]
            for (se, ve), (sg, vg), (sh, vh) in zip(
                self.yield_pair_energies(index1, index2),
                self.yield_pair_gradients(index1, index2),
                self.yield_pair_hessians(index1, index2)
            ):
                result -= (
                    +sh*self.dirouters[index1,index2]*ve
                    +sg*(numpy.identity(3, float) - self.dirouters[index1, index2])*ve*d_1
                    +sg*numpy.outer(self.directions[index1, index2],  vg)
                    +sg*numpy.outer(vg, self.directions[index1, index2])
                    +se*vh
                )
        return result

    def hessian(self):
        result = numpy.zeros((self.numc, self.numc, 3, 3), float)
        for index1 in xrange(self.numc):
            for index2 in xrange(self.numc):
                result[index1, index2] = self.hessian_component(index1, index2)
        return result

    def gradient_flat(self):
        return self.gradient().ravel()

    def hessian_flat(self):
        tmp = self.hessian()
        result = numpy.zeros((self.numc*3,self.numc*3), float)
        for index1 in xrange(self.numc):
            for index2 in xrange(self.numc):
                result[index1*3:(index1+1)*3,index2*3:(index2+1)*3] = tmp[index1,index2]
        return result


class CoulombFF(PairFF):
    def __init__(self, coordinates, charges=None, dipoles=None, exclude_pairs=[]):
        PairFF.__init__(self, coordinates, exclude_pairs)
        self.charges = charges
        self.dipoles = dipoles

    def yield_pair_energies(self, index1, index2):
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
            yield d_3*numpy.dot(p1, p2), 1
            yield -3*d_5, numpy.dot(p1, delta)*numpy.dot(delta, p2)
            if self.charges is not None:
                yield c1*d_3, numpy.dot(p2, delta)
                yield c2*d_3, numpy.dot(p1, -delta)

    def yield_pair_gradients(self, index1, index2):
        d_2 = 1/self.distances[index1, index2]**2
        if self.charges is not None:
            c1 = self.charges[index1]
            c2 = self.charges[index2]
            yield -c1*c2*d_2, numpy.zeros(3)
        if self.dipoles is not None:
            d_4 = d_2**2
            d_6 = d_2**3
            delta = self.deltas[index1, index2]
            p1 = self.dipoles[index1]
            p2 = self.dipoles[index2]
            yield -3*d_4*numpy.dot(p1, p2), numpy.zeros(3)
            yield 15*d_6, p1*numpy.dot(p2, delta) + p2*numpy.dot(p1, delta)
            if self.charges is not None:
                yield -3*c1*d_4, p2
                yield -3*c2*d_4, -p1

    def yield_pair_hessians(self, index1, index2):
        d_1 = 1/self.distances[index1, index2]
        d_3 = d_1**3
        if self.charges is not None:
            c1 = self.charges[index1]
            c2 = self.charges[index2]
            yield 2*c1*c2*d_3, numpy.zeros((3,3))
        if self.dipoles is not None:
            d_5 = d_1**5
            d_7 = d_1**7
            delta = self.deltas[index1, index2]
            p1 = self.dipoles[index1]
            p2 = self.dipoles[index2]
            yield 12*d_5*numpy.dot(p1, p2), numpy.zeros((3,3))
            yield -90*d_7, numpy.outer(p1, p2) + numpy.outer(p2, p1)
            if self.charges is not None:
                yield 12*c1*d_5, numpy.zeros((3,3))
                yield 12*c2*d_5, numpy.zeros((3,3))


class DispersionFF(PairFF):
    def __init__(self, coordinates, strengths, exclude_pairs=[]):
        PairFF.__init__(self, coordinates, exclude_pairs)
        self.strengths = strengths

    def yield_pair_energies(self, index1, index2):
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield strength*distance**(-6), 1

    def yield_pair_gradients(self, index1, index2):
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield -6*strength*distance**(-7), numpy.zeros(3)

    def yield_pair_hessians(self, index1, index2):
        strength = self.strengths[index1, index2]
        distance = self.distances[index1, index2]
        yield 42*strength*distance**(-8), numpy.zeros((3,3))




