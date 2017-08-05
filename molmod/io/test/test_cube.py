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
import pkg_resources

from molmod.test.common import *
from molmod.io import *
from molmod import *


__all__ = ["CubeTestCase"]


class CubeTestCase(BaseTestCase):
    def test_cube_reader(self):
        cr = CubeReader(pkg_resources.resource_filename(__name__, "../../data/test/alanine.cube"))
        self.assertEqual(cr.molecule.numbers[0], 7)
        self.assertEqual(cr.molecule.numbers[8], 1)
        self.assertEqual(cr.molecule.numbers[-1], 8)
        self.assertEqual(cr.nuclear_charges[0], 7.0)
        self.assertEqual(cr.nuclear_charges[8], 1.0)
        self.assertEqual(cr.nuclear_charges[-1], 8.0)
        self.assertAlmostEqual(cr.molecule.coordinates[0,0], -1.020557*angstrom, 5)
        self.assertAlmostEqual(cr.molecule.coordinates[7,2], 0.256423*angstrom, 5)
        vector, value = next(cr)
        self.assertAlmostEqual(value, 1.16792E-14)
        self.assertArraysAlmostEqual(vector, np.array([-9.375592, -8.571340, -6.768095]))
        vector, value = next(cr)
        self.assertAlmostEqual(value, 2.77030E-13)
        self.assertArraysAlmostEqual(vector, np.array([-9.375592, -8.571340, -4.994197]))

    def test_cube_reader_size(self):
        cr = CubeReader(pkg_resources.resource_filename(__name__, "../../data/test/alanine.cube"))
        self.assertEqual(cr.nrep[0], 11)
        self.assertEqual(cr.nrep[1], 10)
        self.assertEqual(cr.nrep[2], 9)
        l = list(cr)
        self.assertEqual(len(cr.f.read().strip()), 0)
        self.assertEqual(len(l), 11*10*9)

    def test_cube(self):
        # Make sure the file is read properly
        cf = Cube.from_file(pkg_resources.resource_filename(__name__, "../../data/test/alanine.cube"))
        self.assertEqual(cf.molecule.numbers[0], 7)
        self.assertEqual(cf.molecule.numbers[8], 1)
        self.assertEqual(cf.molecule.numbers[-1], 8)
        self.assertEqual(cf.nuclear_charges[0], 7.0)
        self.assertEqual(cf.nuclear_charges[8], 1.0)
        self.assertEqual(cf.nuclear_charges[-1], 8.0)
        self.assertAlmostEqual(cf.molecule.coordinates[0,0], -1.020557*angstrom, 5)
        self.assertAlmostEqual(cf.molecule.coordinates[7,2], 0.256423*angstrom, 5)
        self.assertEqual(cf.nrep[0], 11)
        self.assertEqual(cf.nrep[1], 10)
        self.assertEqual(cf.nrep[2], 9)
        self.assertAlmostEqual(cf.data[0,0,0], 1.16792E-14, 17)
        self.assertAlmostEqual(cf.data[0,0,1], 2.77030E-13, 16)
        self.assertAlmostEqual(cf.data[0,0,2], 4.85996E-12, 15)
        self.assertAlmostEqual(cf.data[0,1,0], 1.42126E-12, 15)
        self.assertAlmostEqual(cf.data[0,2,0], 3.72924E-11, 14)
        self.assertAlmostEqual(cf.data[-1,-1,-1], 1.93947E-16, 19)
        # Make sure the file is written properly
        with tmpdir(__name__, 'test_cube') as dn:
            cf.write_to_file('%s/alanine.cube' % dn)
            f1 = open(pkg_resources.resource_filename(__name__, "../../data/test/alanine.cube"))
            lines1 = f1.readlines()
            f1.close()
            f2 = open('%s/alanine.cube' % dn)
            lines2 = f2.readlines()
            f2.close()
        self.assertEqual(len(lines1), len(lines2))
        self.assertEqual(lines1, lines2)
        # Test copying
        cf2 = cf.copy(np.sqrt(cf.data))
        self.assertEqual(cf2.molecule.numbers[0], 7)
        self.assertEqual(cf2.molecule.numbers[8], 1)
        self.assertEqual(cf2.molecule.numbers[-1], 8)
        self.assertEqual(cf2.nuclear_charges[0], 7.0)
        self.assertEqual(cf2.nuclear_charges[8], 1.0)
        self.assertEqual(cf2.nuclear_charges[-1], 8.0)
        self.assertAlmostEqual(cf2.molecule.coordinates[0,0], -1.020557*angstrom, 5)
        self.assertAlmostEqual(cf2.molecule.coordinates[7,2], 0.256423*angstrom, 5)
        self.assertEqual(cf2.nrep[0], 11)
        self.assertEqual(cf2.nrep[1], 10)
        self.assertEqual(cf2.nrep[2], 9)
        self.assertAlmostEqual(cf2.data[0,0,0], np.sqrt(1.16792E-14), 17)
        self.assertAlmostEqual(cf2.data[0,0,1], np.sqrt(2.77030E-13), 16)
        self.assertAlmostEqual(cf2.data[0,0,2], np.sqrt(4.85996E-12), 15)
        self.assertAlmostEqual(cf2.data[0,1,0], np.sqrt(1.42126E-12), 15)
        self.assertAlmostEqual(cf2.data[0,2,0], np.sqrt(3.72924E-11), 14)
        self.assertAlmostEqual(cf2.data[-1,-1,-1], np.sqrt(1.93947E-16), 19)
        # Test the points
        points = cf.get_points()
        self.assertArraysAlmostEqual(points[0,0,0], cf.origin)
        self.assertArraysAlmostEqual(points[1,0,0], cf.origin + cf.axes[0])
        self.assertArraysAlmostEqual(points[0,1,0], cf.origin + cf.axes[1])
        self.assertArraysAlmostEqual(points[0,0,1], cf.origin + cf.axes[2])
        self.assertArraysAlmostEqual(points[5,3,2], cf.origin + 5*cf.axes[0] + 3*cf.axes[1] + 2*cf.axes[2])
