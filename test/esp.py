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


import molmod.esp

import unittest, numpy


__all__ = ["ESPExample"]


class ESPExample(unittest.TestCase):
    def test_esp(self):
        charges = numpy.array([
            [0.4, 1.2, 1.8],
            [1.5, 0.1, 1.0],
            [0.5, 1.1, 0.6],
            [1.3, 1.8, 0.2],
        ], float)
        dipoles = numpy.array([
            [1.3, 1.1, 0.8],
            [0.4, 0.9, 0.1],
            [1.9, 0.1, 1.8],
            [1.5, 1.2, 1.1],
        ], float)

        grid_spacing = 0.2
        grid_origin = -0.05
        grid_size = 11

        potential_data = numpy.zeros((grid_size*grid_size*grid_size, 4), float)
        esp_data = numpy.array([
             0.8, -0.5,  1.3,  0.2,
             0.2,  1.3, -0.5,
            -0.9, -0.1,  1.2,
            -0.3,  0.4, -1.3,
             1.2, -0.8,  0.1,
        ])
        counter = 0
        for x in xrange(grid_size):
            for y in xrange(grid_size):
                for z in xrange(grid_size):
                    point = potential_data[counter, :-1]
                    point[:] = grid_origin + numpy.array([x, y, z])*grid_spacing
                    potential = 0.0
                    for charge, value in zip(charges, esp_data[:4]):
                        distance = numpy.linalg.norm(point - charge)
                        potential += value/distance
                    for dipole, vector in zip(dipoles, esp_data[4:].reshape((4, 3))):
                        delta = point - dipole
                        distance = numpy.linalg.norm(delta)
                        potential += numpy.dot(delta, vector)/distance**3
                    potential_data[counter,-1] = potential
                    counter += 1

        cost_fn = molmod.esp.ESPCostFunction(potential_data, charges, dipoles)
        self.assert_(sum((esp_data -cost_fn.solve())**2) / sum(esp_data**2) < 1e-6, "Wrong ESP solution")
        cost_fn.gradient(esp_data)
        cost_fn.evaluate(esp_data)

