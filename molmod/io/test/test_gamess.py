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


from __future__ import division

import pkg_resources

from molmod.test.common import BaseTestCase
from molmod.io import *
from molmod import *


__all__ = ["GamessTestCase"]


class GamessTestCase(BaseTestCase):
    def test_punch(self):
        punch = PunchFile(pkg_resources.resource_filename(__name__, "../../data/test/PCGamess_PUNCH"))
        self.assertEqual(punch.title, "Simple example sample optimization with Hessian output for Toon")
        self.assertEqual(punch.symmetry, "C1")
        self.assertEqual(punch.symbols, ["CL", "H", "H", "H", "H", "F", "F", "F", "F", "H", "F"])
        N = len(punch.symbols)
        self.assertEqual(punch.numbers.shape, (N,))
        self.assertEqual(punch.numbers[0], 17)
        self.assertEqual(punch.numbers[1], 1)
        self.assertEqual(punch.numbers[-1], 9)
        self.assertEqual(punch.coordinates.shape, (N,3))
        self.assertAlmostEqual(punch.coordinates[0,1]/angstrom, -0.1843157808)
        self.assertAlmostEqual(punch.coordinates[3,-1]/angstrom, 1.2926708150)
        self.assertAlmostEqual(punch.coordinates[-1,0]/angstrom, 3.8608437748)
        self.assertAlmostEqual(punch.energy, -959.9675629527)
        self.assertEqual(punch.gradient.shape, (N,3))
        self.assert_(abs(punch.gradient[0,1] - 1.5314677838E-05) < 1e-10)
        self.assert_(abs(punch.gradient[3,-1] - 8.5221217336E-06) < 1e-10)
        self.assert_(abs(punch.gradient[-1,0] - 2.1211421041E-05) < 1e-10)
        self.assertEqual(punch.hessian.shape, (3*N,3*N))
        self.assert_(abs(punch.hessian - punch.hessian.transpose()).max() < 1e-10)
        self.assert_(abs(punch.hessian[0,0] - 2.51645239E-02) < 1e-10)
        self.assert_(abs(punch.hessian[0,-1] - -1.27201108E-04) < 1e-10)
        self.assert_(abs(punch.hessian[-1,0] - -1.27201108E-04) < 1e-10)
        self.assert_(abs(punch.hessian[-1,-1] - 7.34538698E-03) < 1e-10)
        self.assertEqual(punch.masses.shape, (N,))
        self.assertAlmostEqual(punch.masses[0]/amu, 34.96885)
        self.assertAlmostEqual(punch.masses[3]/amu, 1.00782)
        self.assertAlmostEqual(punch.masses[-1]/amu, 18.99840)
