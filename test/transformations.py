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


from molmod.transformations import Complete

import numpy, copy, random, math
import unittest


__all__ = ["Apply"]


class Apply(unittest.TestCase):
    def setUp(self):
        self.test_transformations = []
        for i in xrange(20):
            test_transformation = Complete()
            test_transformation.set_rotation_properties(random.random()*math.pi*2, numpy.random.uniform(-3, 3, 3), random.sample([True, False], 1)[0])
            test_transformation.t = numpy.random.uniform(-3, 3, 3)
            self.test_transformations.append(test_transformation)

    def test_apply_after(self):
        for tt1 in self.test_transformations:
            for tt2 in self.test_transformations:
                temp = copy.deepcopy(tt1)
                temp.apply_after(tt2)
                r = numpy.dot(tt2.r, tt1.r)
                rotation_error = numpy.sum(numpy.ravel((r - temp.r)**2))/9.0
                self.assertAlmostEqual(rotation_error, 0)
                t = numpy.dot(tt2.r, tt1.t) + tt2.t
                translation_error = numpy.sum((t - temp.t)**2)/3.0
                self.assertAlmostEqual(translation_error, 0)

    def test_apply_before(self):
        for tt1 in self.test_transformations:
            for tt2 in self.test_transformations:
                temp = copy.deepcopy(tt1)
                temp.apply_before(tt2)
                r = numpy.dot(tt1.r, tt2.r)
                rotation_error = numpy.sum(numpy.ravel((r - temp.r)**2))/9.0
                self.assertAlmostEqual(rotation_error, 0)
                t = numpy.dot(tt1.r, tt2.t) + tt1.t
                translation_error = numpy.sum((t - temp.t)**2)/3.0
                self.assertAlmostEqual(translation_error, 0)

    def test_apply_inverse_after(self):
        for tt1 in self.test_transformations:
            for tt2 in self.test_transformations:
                temp = copy.deepcopy(tt1)
                temp.apply_inverse_after(tt2)
                r = numpy.dot(numpy.transpose(tt2.r), tt1.r)
                rotation_error = numpy.sum(numpy.ravel((r - temp.r)**2))/9.0
                self.assertAlmostEqual(rotation_error, 0)
                t = numpy.dot(numpy.transpose(tt2.r), tt1.t - tt2.t)
                translation_error = numpy.sum((t - temp.t)**2)/3.0
                self.assertAlmostEqual(translation_error, 0)

    def test_apply_inverse_before(self):
        for tt1 in self.test_transformations:
            for tt2 in self.test_transformations:
                temp = copy.deepcopy(tt1)
                temp.apply_inverse_before(tt2)
                r = numpy.dot(tt1.r, numpy.transpose(tt2.r))
                rotation_error = numpy.sum(numpy.ravel((r - temp.r)**2))/9.0
                self.assertAlmostEqual(rotation_error, 0)
                t = numpy.dot(tt1.r, -numpy.dot(numpy.transpose(tt2.r), tt2.t)) + tt1.t
                translation_error = numpy.sum((t - temp.t)**2)/3.0
                self.assertAlmostEqual(translation_error, 0)


