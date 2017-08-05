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


import pkg_resources

from molmod.test.common import BaseTestCase
from molmod.io import *
from molmod import *


__all__ = ["CPMDTestCase"]


class CPMDTestCase(BaseTestCase):
    def test_trajectory_reader(self):
        ctr = CPMDTrajectoryReader(pkg_resources.resource_filename(__name__, "../../data/test/TRAJECTORY_H2_CPMD"))
        self.assertEqual(ctr.num_atoms, 2)
        for pos, vel in ctr:
            self.assertAlmostEqual(pos[0,2], 7.55931800897628)
            self.assertAlmostEqual(vel[-1,0], 0.00025263238305)
            break

        ctr = CPMDTrajectoryReader(pkg_resources.resource_filename(__name__, "../../data/test/TRAJECTORY_H2_CPMD"), slice(2,10,2))
        self.assertEqual(ctr.num_atoms, 2)
        pos, vel = next(ctr)
        self.assertAlmostEqual(pos[0,0], 8.28262598957809)
        self.assertAlmostEqual(vel[0,2], 0.00010143872376)
        pos, vel = next(ctr)
        self.assertAlmostEqual(pos[1,1], 7.55433910653197)
        self.assertAlmostEqual(vel[-1,-1], -0.00010042196797)
