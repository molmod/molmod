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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
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
        self.assert_((abs(q1-q2).max() < 1e-10) or (abs(q1+q2).max() < 1e-10))
        self.assert_(abs(r1-r2).max() < 1e-10)

    def test_quaternion_rotation(self):
        q = numpy.random.normal(0,1,4)
        q /= numpy.linalg.norm(q)
        r = quaternion_to_rotation_matrix(q)

        p = numpy.random.normal(0,1,3)
        pa = numpy.dot(r, p)
        pb = quaternion_rotation(q, p)
        self.assert_(abs(pa-pb).max() < 1e-10)

    def test_quaternion_product(self):
        q1 = numpy.random.normal(0,1,4)
        q1 /= numpy.linalg.norm(q1)
        q2 = numpy.random.normal(0,1,4)
        q2 /= numpy.linalg.norm(q2)
        q3 = quaternion_product(q1, q2)
        r1 = quaternion_to_rotation_matrix(q1)
        r2 = quaternion_to_rotation_matrix(q2)
        r3 = numpy.dot(r1, r2)
        foo, q3_check = quaternion_from_rotation_matrix(r3)
        self.assert_((abs(q3-q3_check).max() < 1e-10) or (abs(q3+q3_check).max() < 1e-10))



