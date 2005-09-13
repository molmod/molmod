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

from pychem.interfaces.mpqc import SimpleMpqcJobSinglePoint, SimpleMpqcJobOptimize
from pychem.molecules import molecule_from_xyz_filename

import math, Numeric
import unittest


suite = unittest.TestSuite()

__all__ = ["suite"]


class TestMpqcInterface(unittest.TestCase):
    gridsize = 1.0
    
    def test_single_point(self):
        def validate():
            self.assert_(job.completed)
            self.assertAlmostEqual(job.energy, -75.9734488121, 8)
            self.assertAlmostEqual(job.gradient[0,0],  0.01174361, 6)
            self.assertAlmostEqual(job.gradient[0,1],  0.0,        6)
            self.assertAlmostEqual(job.gradient[0,2],  0.0,        6)
            self.assertAlmostEqual(job.gradient[1,0], -0.0058718,  6)
            self.assertAlmostEqual(job.gradient[1,1],  0.0,        6)
            self.assertAlmostEqual(job.gradient[1,2], -0.01381411, 6)
            self.assertAlmostEqual(job.gradient[2,0], -0.0058718,  6)
            self.assertAlmostEqual(job.gradient[2,1],  0.0,        6)
            self.assertAlmostEqual(job.gradient[2,2],  0.01381411, 6)
            
        water = molecule_from_xyz_filename("input/water.xyz")
        job = SimpleMpqcJobSinglePoint(
            "output/water_sp",
            "Water single point berekening", 
            water,
            "KS (xc = B3LYP)",
            "3-21G*",
            do_gradient=True
        )
        job.run(user_overwrite=True)
        validate()
        job.run()
        validate()
                  
    def test_optimize(self):
        def validate():
            self.assert_(job.completed)
            self.assertAlmostEqual(job.energies[-1], -75.973963163199997, 8)
            coordinates = job.output_molecule.coordinates
            delta = coordinates[0]-coordinates[1]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 1.88335259871, 6)
            delta = coordinates[0]-coordinates[2]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 1.88335259871, 6)
            delta = coordinates[1]-coordinates[2]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 2.96668446577, 6)
        
        water = molecule_from_xyz_filename("input/water.xyz")
        job = SimpleMpqcJobOptimize(
            "output/water_opt",
            "Water single point berekening", 
            water,
            "KS (xc = B3LYP)",
            "3-21G*"
        )
        job.run(user_overwrite=True)
        validate()
        job.run()
        validate()
        

suite.addTest(unittest.makeSuite(TestMpqcInterface))
