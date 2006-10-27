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


import numpy, math

from molmod.helpers import potential


class ESPCostFunction(object):
    def __init__(self, charges, dipoles, design_matrix, expected_values, weights=None):
        self.charges = charges
        self.dipoles = dipoles
        self.design_matrix = design_matrix
        self.expected_values = expected_values
        self.weights = weights
        self.num_equations = len(expected_values)
        self.recompute()

    def recompute(self):
        if self.weights is not None:
            tmp_ev = self.weights*self.expected_values
            tmp_dm = (self.design_matrix.transpose()*self.weights).transpose()
        else:
            tmp_ev = self.expected_values
            tmp_dm = self.design_matrix
        self.A = numpy.dot(tmp_dm.transpose(), tmp_dm)
        self.B = numpy.dot(tmp_dm.transpose(), tmp_ev)
        self.C = numpy.dot(tmp_ev, tmp_ev)

    def evaluate(self, unknowns):
        return numpy.dot(numpy.dot(unknowns, self.A) - 2*self.B, unknowns) + self.C

    def gradient(self, unknowns):
        return 2*(numpy.dot(unknowns, self.A) - self.B)

    def solve(self):
        return numpy.linalg.solve(self.A, self.B)

    def solve_constrained(self):
        projection = numpy.zeros((len(self.charges) + 3*len(self.dipoles) - 1, len(self.charges) + 3*len(self.dipoles)), float)
        projection[0:len(self.charges)-1,0] = -1
        for index in xrange(len(self.charges) + 3*len(self.dipoles) - 1):
            projection[index, index+1] = 1
        A = numpy.dot(projection, numpy.dot(self.A, projection.transpose()))
        B = numpy.dot(projection, self.B)
        return numpy.dot(projection.transpose(), numpy.linalg.solve(A, B))

    def cost_to_error(self, cost):
        return math.sqrt(cost/self.num_equations)


class StandardESPCostFunction(ESPCostFunction):
    def __init__(self, points, potential, charges, dipoles, weights=None, external_potential=None, epsilons=None):
        if points is not None:
            design_matrix = numpy.zeros((len(potential), len(charges) + len(dipoles)*3), float)
            expected_values = potential.copy()

            if external_potential is not None:
                expected_values -= external_potential(points)

            for row, point in enumerate(points):
                if len(charges) > 0:
                    design_matrix[row,:len(charges)] = 1/numpy.sqrt(sum((point - charges).transpose()**2))
                if len(dipoles) > 0:
                    design_matrix[row,len(charges):] = ((point - dipoles).transpose()/numpy.sqrt(sum((point - dipoles).transpose()**2)**3)).transpose().ravel()

            ESPCostFunction.__init__(self, charges, dipoles, design_matrix, expected_values, weights)

            if epsilons is not None:
                self.A[::len(self.A)+1] += epsilons


class SphericallySwitchedESPCostFunction(ESPCostFunction):
    def __init__(self, cube, samples, charges, dipoles):
        assert cube.volumes is not None
        ion_coordinates = cube.molecule.coordinates
        ion_charges = cube.molecule.numbers

        self.num_equations = len(samples)

        design_matrix = numpy.zeros((self.num_equations, len(charges) + len(dipoles)*3), float)
        expected_values = numpy.zeros(self.num_equations, float)
        self.electron_depletion = numpy.zeros(self.num_equations, float)
        self.ion_depletion = numpy.zeros(self.num_equations, float)
        self.closest_distances = numpy.zeros(self.num_equations, float)

        self.all_distances = numpy.zeros((len(cube.grid_data), len(cube.molecule.numbers)), float)
        #potential.all_distances(
        for index, coordinate in enumerate(cube.molecule.coordinates):
            self.all_distances[:,index] = numpy.sqrt(((cube.grid_data[:,:3] - coordinate)**2).sum(axis=1))
        closest_indices = self.all_distances.argmin(axis=1)


        for row, (point, switch_low, switch_high) in enumerate(samples):
            if row % 100 == 0: print "%06i/%06i" % (row, len(samples))
            distances = numpy.sqrt(((cube.molecule.coordinates - point)**2).sum(axis=1))
            self.closest_distances[row] = distances.min()
            #print "  point:", point
            v_ai, qe, qi, ionic_s = potential.spherical_switchv(
                point, cube.grid_data, cube.volumes,
                ion_charges, ion_coordinates,
                switch_low, switch_high,
                closest_indices
            )
            #print " ".join("% 4i" % int(s*100) for s in ionic_s)
            self.electron_depletion[row] = qe
            self.ion_depletion[row] = qi
            #print "Potential:", v_ai
            #print "Net charge:", q_net
            counter = 0
            #print " charges"
            for charge, s in zip(charges, ionic_s):
                distance = numpy.linalg.norm(point - charge)
                #s = potential.switch_cos(distance, switch_low, switch_high)
                design_matrix[row, counter] = s/distance
                counter += 1
            #print " dipoles"
            for dipole, s in zip(dipoles, ionic_s):
                delta = point - dipole
                distance = numpy.linalg.norm(delta)
                #s = potential.switch_cos(distance, switch_low, switch_high)
                design_matrix[row, counter:counter+3] = delta*s/distance**3
                counter += 3
            expected_values[row] = v_ai

        ESPCostFunction.__init__(self, charges, dipoles, design_matrix, expected_values)


class LinearlySwitchedESPCostFunction(ESPCostFunction):
    def __init__(self, cube, samples, charges, dipoles):
        assert cube.volumes is not None
        ion_coordinates = cube.molecule.coordinates
        ion_charges = cube.molecule.numbers

        self.num_equations = len(samples)
        self.num_coeffs = len(charges) + len(dipoles)*3

        design_matrix = numpy.zeros((self.num_equations, self.num_coeffs), float)
        expected_values = numpy.zeros(self.num_equations, float)

        def equation(observer, direction, length, switch_low, switch_high):
            coefficients = numpy.zeros(self.num_coeffs, float)
            v_ai, qe, qi = potential.linear_switchv(
                observer, direction, length,
                cube.grid_data, cube.volumes,
                ion_charges, ion_coordinates,
                switch_low, switch_high,
            )
            #print "Potential:", v_ai
            #print "Net charge:", q_net
            counter = 0
            #ssum = 0
            for charge in charges:
                distance = numpy.linalg.norm(observer - charge)
                s = potential.switch_cos(potential.measure_ortho(charge, observer, direction, length), switch_low, switch_high)
                if s > 0:
                    coefficients[counter] = s/distance
                counter += 1
                #ssum += s
            for dipole in dipoles:
                distance = numpy.linalg.norm(observer - dipole)
                s = potential.switch_cos(potential.measure_ortho(dipole, observer, direction, length), switch_low, switch_high)
                if s > 0:
                    coefficients[counter:counter+3] = delta*s/distance**3
                counter += 3
                #ssum += s
            #print "ssum:", ssum
            #print
            return coefficients, v_ai

        for row, (a, b, switch_low, switch_high) in enumerate(samples):
            if row % 100 == 0: print "%06i/%06i" % (row, len(samples))
            distance = numpy.linalg.norm(b - a)
            direction = (b - a)/distance
            coefficients1, expected1 = equation(a, direction, distance, switch_low, switch_high)
            coefficients2, expected2 = equation(b, -direction, distance, switch_low, switch_high)
            design_matrix[row] = coefficients1 + coefficients2
            expected_values[row] = expected1 + expected2

        ESPCostFunction.__init__(self, charges, dipoles, design_matrix, expected_values)


