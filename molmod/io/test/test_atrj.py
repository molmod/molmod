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

import unittest

import pkg_resources

from molmod.test.common import BaseTestCase
from molmod.io import *
from molmod import *


__all__ = ["ATRJTestCase"]


class ATRJTestCase(BaseTestCase):
    def test_load(self):
        # A) normal
        atrj_reader = ATRJReader(pkg_resources.resource_filename(__name__, "../../data/test/bartek.atrj"))
        self.assertEqual(atrj_reader.num_atoms, 1293)
        frames = list(atrj_reader)
        self.assertEqual(len(frames), 3)
        # check time
        self.assertAlmostEqual(frames[0].time/picosecond, 1.0)
        self.assertAlmostEqual(frames[1].time/picosecond, 2.0)
        self.assertAlmostEqual(frames[2].time/picosecond, 3.0)
        # check step
        self.assertAlmostEqual(frames[0].step, 1000)
        self.assertAlmostEqual(frames[1].step, 2000)
        self.assertAlmostEqual(frames[2].step, 3000)
        # check total energy
        self.assertAlmostEqual(frames[0].total_energy/kcalmol, 3.4186035768405162e2)
        self.assertAlmostEqual(frames[1].total_energy/kcalmol, 3.3443356630787252e2)
        self.assertAlmostEqual(frames[2].total_energy/kcalmol, 3.3613629561285467e2)
        # check (some of) the coordinates
        self.assertAlmostEqual(frames[0].coordinates[0,0]/angstrom, 1.1953453341349823e1)
        self.assertAlmostEqual(frames[1].coordinates[5,1]/angstrom, 1.2284911431298470e1)
        self.assertAlmostEqual(frames[-1].coordinates[-5,-1]/angstrom, 2.1392983758428979e1)

        # B) sliced
        atrj_reader = ATRJReader(pkg_resources.resource_filename(__name__, "../../data/test/bartek.atrj"), slice(None, None, 2))
        self.assertEqual(atrj_reader.num_atoms, 1293)
        frames = list(atrj_reader)
        self.assertEqual(len(frames), 2)
        # check time
        self.assertAlmostEqual(frames[0].time/picosecond, 1.0)
        self.assertAlmostEqual(frames[1].time/picosecond, 3.0)
