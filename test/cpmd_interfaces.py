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

from molmod.interfaces.cpmd.file_parsers import *
from molmod.interfaces.output_parsers import OutputParser


import math, numpy, numpy.linalg
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
            MassesParser(),
            CoordinatesGradientsParser(num_steps_parser, num_every_parser, elements_parser),
            EnergiesParser(num_steps_parser, num_every_parser),
            EnergiesFileParser(num_steps_parser),
            TrajectoryFileParser(num_steps_parser, elements_parser),
            CellDimensionParser()
        ])

        result = output_parser.parse("input/cpmd", "2TOH.md")
        
        expected_keys = ['coor_grad', 'elements', 'energies', 'time_steps', 'num_steps', 'num_every', 'energies_table', 'coor_velo', 'masses']
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
        #print result["energies_table"][:10]
        #coordinates, velocities = result["coor_velo"]
        #print "coordinates:", coordinates[0:2]
        #print "velocities:", velocities[0:2]
        #print "masses", result["masses"]
        #print "cell_dimension", result["cell_dimension"]
        
