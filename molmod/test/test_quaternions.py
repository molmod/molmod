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

import unittest

import numpy as np

from molmod import *


__all__ = ["QuaternionTestCase"]


class QuaternionTestCase(unittest.TestCase):
    def test_rotation_matrix(self):
        q1 = np.random.normal(0,1,4)
        q1 /= np.linalg.norm(q1)
        r1 = quaternion_to_rotation_matrix(q1)
        factor, q2 = rotation_matrix_to_quaternion(r1)
        r2 = quaternion_to_rotation_matrix(q2)
        self.assert_((abs(q1-q2).max() < 1e-10) or (abs(q1+q2).max() < 1e-10))
        self.assert_(abs(r1-r2).max() < 1e-10)

    def test_quaternion_rotation(self):
        q = np.random.normal(0,1,4)
        q /= np.linalg.norm(q)
        r = quaternion_to_rotation_matrix(q)

        p = np.random.normal(0,1,3)
        pa = np.dot(r, p)
        pb = quaternion_rotation(q, p)
        self.assert_(abs(pa-pb).max() < 1e-10)

    def test_quaternion_product(self):
        q1 = np.random.normal(0,1,4)
        q1 /= np.linalg.norm(q1)
        q2 = np.random.normal(0,1,4)
        q2 /= np.linalg.norm(q2)
        q3 = quaternion_product(q1, q2)
        r1 = quaternion_to_rotation_matrix(q1)
        r2 = quaternion_to_rotation_matrix(q2)
        r3 = np.dot(r1, r2)
        foo, q3_check = rotation_matrix_to_quaternion(r3)
        self.assert_((abs(q3-q3_check).max() < 1e-10) or (abs(q3+q3_check).max() < 1e-10))
