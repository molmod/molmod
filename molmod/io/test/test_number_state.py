# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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


import numpy as np

from molmod.test.common import *
from molmod.io import *


__all__ = ["NumberStateTestCase"]


class TestObject(object):
    def __init__(self):
        self.a = np.random.normal(0, 1, (5, 3))
        self.b = np.random.randint(0, 40, (4, 7, 9))
        self.c = np.random.normal(0, 2)
        self.d = np.random.randint(0, 10)
        self.state = NumberState(self, ["a", "b", "c", "d"])


class NumberStateTestCase(BaseTestCase):
    def test_consistency(self):
        with tmpdir(__name__, 'test_consistency') as dn:
            test1 = TestObject()
            test1.state.dump("%s/test" % dn)
            test2 = TestObject()
            test2.state.load("%s/test" % dn)
            self.assertArraysAlmostEqual(test1.a, test2.a, 1e-10)
            self.assertArraysAlmostEqual(test1.b, test2.b, 1e-10)
            self.assertAlmostEqual(test1.c, test2.c)
            self.assertAlmostEqual(test1.d, test2.d)
            test2.a[:] = 0.0
            test2.b[:] = 0
            test2.c = 0.0
            test2.d = 0
            test2.state.load("%s/test" % dn, ["a", "d"])
            self.assertArraysAlmostEqual(test1.a, test2.a, 1e-10)
            self.assertAlmostEqual(abs(test2.b).max(), 0.0)
            self.assertAlmostEqual(test2.c, 0.0)
            self.assertAlmostEqual(test1.d, test2.d)

        test1 = TestObject()
        test2 = TestObject()
        subset = set(["a", "d"])
        state = test1.state.get(subset)
        self.assertEqual(len(state), 2)
        self.assert_("a" in state)
        self.assert_("d" in state)
        test2.state.set(state, subset)
        self.assertArraysAlmostEqual(test1.a, test2.a, 1e-10)
        self.assertAlmostEqual(test1.d, test2.d)
