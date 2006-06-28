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


import math, numpy


class SpringFF(object):
    def __init__(self, coordinates):
        self.update_coordinates(coordinates)

    def update_coordinates(self, coordinates):
        self.coordinates = coordinates
        self.numc = len(self.coordinates)
        self.distances = numpy.zeros((self.numc, self.numc), float)
        self.directions = numpy.zeros((self.numc, self.numc, 3), float)
        for index1, coordinate1 in enumerate(self.coordinates):
            for index2, coordinate2 in enumerate(self.coordinates):
                delta = coordinate1 - coordinate2
                distance = math.sqrt(numpy.dot(delta, delta))
                self.distances[index1, index2] = distance
                if index1 != index2:
                    self.directions[index1, index2] = delta/distance

    def spring_energy(self, distance, index1, index2):
        raise NotImplementedError

    def spring_gradient(self, distance, index1, index2):
        raise NotImplementedError

    def spring_hessian(self, distance, index2):
        raise NotImplementedError

    def energy(self):
        result = 0.0
        for index1 in xrange(self.numc):
            for index2 in xrange(index1):
                result += self.spring_energy(self.distances[index1, index2], index1, index2)
        return result

    def gradient(self):
        result = numpy.zeros((self.numc, 3), float)
        for index1 in xrange(self.numc):
            for index2 in xrange(self.numc):
                if index2 != index1:
                    result[index1] += self.spring_gradient(self.distances[index1, index2], index1, index2)*self.directions[index1,index2]
        return result

    def hessian(self):
        result = numpy.zeros((self.numc, self.numc, 3, 3), float)
        for index1 in xrange(self.numc):
            for index2 in xrange(self.numc):
                if index1 == index2:
                    for index3 in xrange(self.numc):
                        if index3 != index1:
                            distance = self.distances[index1, index3]
                            result[index1, index1] += (
                                (
                                    self.spring_gradient(distance, index1, index3)
                                    *(numpy.identity(3, float) - numpy.outer(self.directions[index1,index3], self.directions[index1,index3]))
                                    /distance
                                ) + (
                                    self.spring_hessian(distance, index1, index3)
                                    *numpy.outer(self.directions[index1,index3], self.directions[index1,index3])
                                )
                            )
                else:
                    distance = self.distances[index1, index2]
                    result[index1, index2] = -(
                        (
                            self.spring_gradient(distance, index1, index2)
                            *(numpy.identity(3, float) - numpy.outer(self.directions[index1,index2], self.directions[index1,index2]))
                            /distance
                        ) + (
                            self.spring_hessian(distance, index1, index2)
                            *numpy.outer(self.directions[index1,index2], self.directions[index1,index2])
                        )
                    )
        return result

    def gradient_float(self):
        return self.gradient().ravel()

    def hessian_flat(self):
        tmp = self.hessian()
        result = numpy.zeros((self.numc*3,self.numc*3), float)
        for index1 in xrange(self.numc):
            for index2 in xrange(self.numc):
                result[index1*3:(index1+1)*3,index2*3:(index2+1)*3] = tmp[index1,index2]
        return result
                

class CoulombFF(SpringFF):
    def __init__(self, coordinates, charges):
        SpringFF.__init__(self, coordinates)
        self.charges = charges

    def spring_energy(self, distance, index1, index2):
        return self.charges[index1]*self.charges[index2]/distance

    def spring_gradient(self, distance, index1, index2):
        return -self.charges[index1]*self.charges[index2]*distance**(-2)

    def spring_hessian(self, distance, index1, index2):
        return 2*self.charges[index1]*self.charges[index2]*distance**(-3)


class DispersionFF(SpringFF):
    def __init__(self, coordinates, strength_matrix):
        SpringFF.__init__(self, coordinates)
        self.strength_matrix = strength_matrix

    def spring_energy(self, distance, index1, index2):
        return self.strength_matrix[index1, index2]*distance**(-6)

    def spring_gradient(self, distance, index1, index2):
        return (-7)*self.strength_matrix[index1, index2]*distance**(-7)

    def spring_hessian(self, distance, index1, index2):
        return (42)*self.strength_matrix[index1, index2]*distance**(-8)


