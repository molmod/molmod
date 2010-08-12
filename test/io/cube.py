# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
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


from common import BaseTestCase

from molmod.io.cube import *
from molmod.units import *

import numpy


__all__ = ["CubeTestCase"]


class CubeTestCase(BaseTestCase):
    def test_cube_reader(self):
        cr = CubeReader("input/alanine.cube")
        self.assertEqual(cr.numbers[0], 7)
        self.assertEqual(cr.numbers[8], 1)
        self.assertEqual(cr.numbers[-1], 8)
        self.assertAlmostEqual(cr.coordinates[0,0], -1.020557*angstrom, 5)
        self.assertAlmostEqual(cr.coordinates[7,2], 0.256423*angstrom, 5)
        vector, value = cr.next()
        self.assertAlmostEqual(value, 1.16792E-14)
        self.assertArraysAlmostEqual(vector, numpy.array([-9.375592, -8.571340, -6.768095]))
        vector, value = cr.next()
        self.assertAlmostEqual(value, 2.77030E-13)
        self.assertArraysAlmostEqual(vector, numpy.array([-9.375592, -8.571340, -4.994197]))
