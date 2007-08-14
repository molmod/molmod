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

import unittest, numpy

from molmod.quaternions import *


__all__ = ["QuaternionTestCase"]


class QuaternionTestCase(unittest.TestCase):
    def test_rotation_matrix(self):
        q1 = numpy.random.normal(0,1,4)
        q1 /= numpy.linalg.norm(q1)
        r1 = quaternion_to_rotation_matrix(q1)
        factor, q2 = quaternion_from_rotation_matrix(r1)
        r2 = quaternion_to_rotation_matrix(q2)
        self.assert_(abs(q1-q2).max() < 1e-10)
        self.assert_(abs(r1-r2).max() < 1e-10)

