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

from pychem.interfaces.cpmd.file_parsers import *
from pychem.interfaces.output_parsers import OutputParser


import math, Numeric, LinearAlgebra
import unittest


__all__ = ["CpmdInterface"]


class CpmdInterface(unittest.TestCase):
    def test_md_parser(self):
        num_steps_parser = NumStepsParser()
        num_every_parser = NumEveryParser()
        elements_parser = ElementsParser()
        
        output_parser = OutputParser([
            num_steps_parser,
            num_every_parser,
            elements_parser,
            TimeStepsParser(),
            CoordinatesGradientsParser(num_steps_parser, num_every_parser, elements_parser),
            EnergiesParser(num_steps_parser, num_every_parser)
        ])

        result = output_parser.parse("input", "cpmd_md1")
        
        expected_keys = ['coor_grad', 'elements', 'energies', 'time_steps', 'num_steps', 'num_every']
        for expected_key in expected_keys:
            self.assert_(expected_key in result)

        #print "num_steps:", result["num_steps"]
        #print "num_every:", result["num_every"]
        #coordinates, gradients = result["coor_grad"]
        #print "coordinates:", coordinates[0]
        #print "gradient:", gradients[0]
        #print "energy:", result["energies"][0]
        #print "elements:", result["elements"]
        #print "time_steps", result["time_steps"]
        
