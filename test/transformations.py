# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

from molmod.transformations import Complete, coincide
from molmod.vectors import random_unit

import numpy, copy, random, math
import unittest


__all__ = ["TransformationsTestCase"]


class TransformationsTestCase(BaseTestCase):
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

    def test_coincide(self):
        # create a few test sets with random data points, including degenerate
        # situations. (e.g. one point, two points, linear points)
        references = [ # a list of 2-tuples: (points, degenerate)
            (numpy.random.normal(0, 5, (n, 3)), False) for n in xrange(4, 40)
        ] + [
            (numpy.array([[0,0,1]], float), True),
            (numpy.array([[0,0,0],[0,0,1]], float), True),
            (numpy.array([[0,0,0],[0,0,1],[0,0,2]], float), True),
            (numpy.array([[0,0,0],[0,0,1],[0,0,2],[0,0,4]], float), True),
            (numpy.random.normal(0, 5, (3, 3)), True)
        ]

        # Do a random transformation on the points
        randomized = []
        for points, degenerate in references:
            axis = random_unit(3)
            angle = numpy.random.uniform(0, numpy.pi*2)
            transformation = Complete()
            transformation.set_rotation_properties(angle, axis, False)
            transformation.t[:] = 0#numpy.random.normal(0, 5, 3)
            randomized.append((
                numpy.array([transformation.vector_apply(p) for p in points]),
                transformation
            ))

        for (ref_points, degenerate), (tr_points, transf) in zip(references, randomized):
            check_transf = coincide(ref_points, tr_points)
            # first check whether check_transf brings the tr_points back to the ref_points
            check_points = numpy.dot(tr_points, check_transf.r.transpose()) + check_transf.t
            self.assertArraysAlmostEqual(ref_points, check_points, 1e-5)
            # check whether the rotation matrix is orthogonal
            self.assertArraysAlmostEqual(numpy.dot(check_transf.r, check_transf.r.transpose()), numpy.identity(3, float), 1e-5)
            if not degenerate:
                # if the set of points is not degenerate, check whether transf and check_transf
                # are each other inverses
                tmp = Complete()
                tmp.apply_before(transf)
                tmp.apply_before(check_transf)
                self.assertArraysAlmostEqual(numpy.dot(transf.r, check_transf.r), numpy.identity(3, float), 1e-5)
                self.assertArrayAlmostZero(tmp.t, 1e-5)


        # Add some distortion to the transformed points
        randomized = []
        for tr_points, transf in randomized:
            tr_points[:] += numpy.random.normal(0, 1.0, len(tr_points))

        # Do a blind test
        for (ref_points, degenerate), (tr_points, transf) in zip(references, randomized):
            coincide(ref_points, tr_points)



