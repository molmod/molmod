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


__all__ = ["LAMMPSTestCase"]


class LAMMPSTestCase(BaseTestCase):
    def test_dump_reader(self):
        ldr = LAMMPSDumpReader(pkg_resources.resource_filename(__name__, "../../data/test/lammps_dump.txt"), [angstrom]*3 + [angstrom/femtosecond]*3)
        for fields in ldr:
            self.assertAlmostEqual(fields[0], 0)
            self.assertAlmostEqual(fields[1][11]/angstrom, -1.68253)
            self.assertAlmostEqual(fields[2][3]/angstrom, -2.71309)
            self.assertAlmostEqual(fields[3][9]/angstrom, -1.81368)
            self.assertAlmostEqual(fields[4][0]/angstrom*femtosecond, -0.00598226)
            self.assertAlmostEqual(fields[5][2]/angstrom*femtosecond, 0.00113906)
            self.assertAlmostEqual(fields[6][8]/angstrom*femtosecond, -0.00390626)
            break

        ldr = LAMMPSDumpReader(pkg_resources.resource_filename(__name__, "../../data/test/lammps_dump.txt"), [angstrom]*3 + [angstrom/femtosecond]*3)
        self.assertEqual(len(list(ldr)), 7)

        ldr = LAMMPSDumpReader(pkg_resources.resource_filename(__name__, "../../data/test/lammps_dump.txt"), [angstrom]*3 + [angstrom/femtosecond]*3, sub=slice(1,5,2))
        self.assertEqual(len(list(ldr)), 2)
