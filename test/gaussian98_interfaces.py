 # PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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

from pychem.interfaces.gaussian98.file_parsers import *
from pychem.interfaces.output_parsers import OutputParser


import math, numpy, numpy.linalg
import unittest


__all__ = ["Gaussian98Interface"]


class Gaussian98Interface(unittest.TestCase):
    def test_parser(self):
        output_parser = OutputParser([
            HessianParser(),
            MassParser(),
            CoordinatesParser(),
            EnergyParser(),
            LowFrequenciesParser(),
            SelectedFrequenciesParser(),
            GradientParser(),
        ])

        result = output_parser.parse("input", "g98_1")

        expected_keys = ['energies', 'masses', 'coordinates_list', 'low_frequencies', 'selected_frequencies', 'hessian', 'gradient_list']
        for expected_key in expected_keys:
            self.assert_(expected_key in result)
            #print "'%s': %s" % (expected_key, results[expected_key])

        #hessian = result["hessian"]
        #masses = result["masses"]
        #coordinates_list = result["coordinates_list"]

        #print "symmetry hessian:", sum(numpy.ravel(hessian - numpy.transpose(hessian))**2)
        
        #print "DOF hess:   % 3i, % 3i" % (hessian.shape[0], hessian.shape[1])
        #print "DOF masses: % 3i" % (len(masses)*3)
        #for coordinates in coordinates_list:
        #    print "DOF coords: % 3i" % (len(coordinates)*3)
        
