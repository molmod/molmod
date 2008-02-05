# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

from molmod.io.array_state import *

import numpy, unittest


__all__ = ["ArrayStateTestCase"]


class TestObject(object):
    def __init__(self):
        self.a = numpy.random.normal(0, 1, (5, 3))
        self.b = numpy.random.randint(0, 40, (4, 7, 9))
        self.state = ArrayState()
        self.state.set_field("a", self.a)
        self.state.set_field("b", self.b)


class ArrayStateTestCase(BaseTestCase):
    def test_consistency(self):
        test1 = TestObject()
        test1.state.dump("output/test")
        test2 = TestObject()
        test2.state.load("output/test")
        self.assertArraysAlmostEqual(test1.a, test2.a, 1e-10)
        self.assertArraysAlmostEqual(test1.b, test2.b, 1e-10)
        test2.b[:] = 0.0
        test2.state.load("output/test", ["a"])
        self.assertAlmostEqual(abs(test2.b).max(), 0.0)
