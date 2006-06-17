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


from molmod.lone_pair import all_lone_pairs
from molmod.molecules import molecule_xyz_from_filename

import unittest

__all__ = ["LonePairs"]


class LonePairs(unittest.TestCase):
    gridsize = 1.0
    
    def molecule_test(self, filename):
        m = molecule_xyz_from_filename("input/"+filename)
        print all_lone_pairs(m)
                
    def test_water(self):
        self.molecule_test("water.xyz")
