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


from molmod.io.gaussian03 import *
from molmod.io.output_parsers import OutputParser

import numpy

import unittest


__all__ = ["Gaussian03"]


class Gaussian03(unittest.TestCase):
    def test_parsers(self):
        output_parser = OutputParser([
            HessianParser(),
            MassParser(),
            StandardOrientationCoordinatesParser(),
            InputOrientationCoordinatesParser(),
            EnergyParser(),
            LowFrequenciesParser(),
            SelectedFrequenciesParser(),
            StandardOrientationGradientParser(),
            InputOrientationGradientParser(),
            IsOptimizedParser(),
            OptimizedCoordinatesParser(),
        ])

        result = output_parser.parse("input/g98_1")

        expected_keys = [
            'energies',
            'masses',
            'so_coordinates_list',
            'io_coordinates_list',
            'low_frequencies',
            'selected_frequencies',
            'hessian',
            'so_gradient_list',
            'io_gradient_list',
            'optimized',
            'optimized_coordinates'
        ]

        for expected_key in expected_keys:
            self.assert_(expected_key in result)
            #print "'%s': %s" % (expected_key, result[expected_key])

        #hessian = result["hessian"]
        #masses = result["masses"]
        #coordinates_list = result["coordinates_list"]

        #print "symmetry hessian:", sum(numpy.ravel(hessian - numpy.transpose(hessian))**2)

        #print "DOF hess:   % 3i, % 3i" % (hessian.shape[0], hessian.shape[1])
        #print "DOF masses: % 3i" % (len(masses)*3)
        #for coordinates in coordinates_list:
        #    print "DOF coords: % 3i" % (len(coordinates)*3)

    def test_fchk(self):
        fchk = FCHKFile("input/1TOH.b3lyp.fchk")
        fchk = FCHKFile("input/1TOH.b3lyp.trim.fchk", ignore_errors=True)
        #print fchk.molecule.numbers
        #print fchk.molecule.coordinates
        #print fchk.optimization_coordinates()
        #print fchk.optimized_molecule()







