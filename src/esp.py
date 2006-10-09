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


import numpy


class ESPCostFunction(object):
    def __init__(self, potential_data, charges, dipoles, weights=None, external_potential=None, epsilons=None):
        if potential_data is not None:
            self.charges = charges
            self.dipoles = dipoles

            design_matrix = numpy.zeros((len(potential_data), len(charges) + len(dipoles)*3), float)
            expected_values = potential_data[:,-1].copy()

            if external_potential is not None:
                expected_values -= external_potential(potential_data[:,:-1])

            for row, record in enumerate(potential_data):
                point = record[:-1]
                if len(charges) > 0:
                    design_matrix[row,:len(charges)] = 1/numpy.sqrt(sum((point - charges).transpose()**2))
                if len(dipoles) > 0:
                    design_matrix[row,len(charges):] = ((point - dipoles).transpose()/numpy.sqrt(sum((point - dipoles).transpose()**2)**3)).transpose().ravel()

            if weights is not None:
                expected_values *= weights
                tmp = design_matrix.transpose()
                tmp *= weights

            self.A = numpy.dot(design_matrix.transpose(), design_matrix)
            self.B = numpy.dot(design_matrix.transpose(), expected_values)
            self.C = numpy.dot(expected_values, expected_values)

            if epsilons is not None:
                self.A[::len(self.A)+1] += epsilons



    def evaluate(self, unknowns):
        return numpy.dot(numpy.dot(unknowns, self.A) - 2*self.B, unknowns) + self.C

    def gradient(self, unknowns):
        return 2*(numpy.dot(unknowns, self.A) - self.B)

    def solve(self):
        return numpy.linalg.solve(self.A, self.B)

    def solve_constrained(self):
        projection = numpy.array((len(self.charges) + len(self.dipoles) - 1, len(self.charges) + len(self.dipoles)), float)
        projection[0:len(charges),0] = -1
        for index in xrange(len(charges) + len(dipoles) - 1):
            projection[index+1, index] = 1

        A = numpy.dot(projection, numpy.dot(self.A, projection.transpose()))
        B = numpy.dot(projection, B)
        return numpy.dot(projection.transpose(), numpy.linalg.solve(A, B))
