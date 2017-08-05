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


from builtins import range
import numpy as np

from molmod.test.common import BaseTestCase
from molmod import *


__all__ = ["TransformationsTestCase"]


class TransformationsTestCase(BaseTestCase):
    def iter_random_translations(self, n=20):
        for i in range(n):
            yield Translation(
                np.random.uniform(-3, 3, 3)
            )

    def iter_random_rotations(self, n=20):
        for i in range(n):
            yield Rotation.random()

    def iter_random_completes(self, n=20):
        for i in range(n):
            yield Complete.from_properties(
                np.random.uniform(0,np.pi*2),
                random_unit(),
                np.random.randint(0,1),
                np.random.uniform(-3, 3, 3)
            )

    def test_translation_apply_to_vector(self):
        for t in self.iter_random_translations():
            x = np.random.uniform(-3,3,3)
            self.assertArraysAlmostEqual(t.inv*x, x-t.t)
            self.assertArraysAlmostEqual(t.inv*(t*x), x)
            self.assertArraysAlmostEqual(t*(t.inv*x), x)
            self.assertArraysAlmostEqual((t*t.inv)*x, x)
            self.assertArraysAlmostEqual((t.inv*t)*x, x)
            self.assertArraysAlmostEqual(t.apply_to(x), x+t.t)
            self.assertArraysAlmostEqual(t*x, x+t.t)

    def test_translation_apply_to_vectors(self):
        for t in self.iter_random_translations():
            x = np.random.uniform(-3,3,(10,3))
            self.assertArraysAlmostEqual(t.inv*x, x-t.t)
            self.assertArraysAlmostEqual(t.inv*(t*x), x)
            self.assertArraysAlmostEqual(t*(t.inv*x), x)
            self.assertArraysAlmostEqual((t*t.inv)*x, x)
            self.assertArraysAlmostEqual((t.inv*t)*x, x)
            self.assertArraysAlmostEqual(t.apply_to(x), x+t.t)
            self.assertArraysAlmostEqual(t*x, x+t.t)
            self.assertArraysAlmostEqual(t.apply_to(x).transpose(), t.apply_to(x.transpose(), columns=True))

    def test_translation_apply_to_translation(self):
        for t1, t2 in zip(self.iter_random_translations(), self.iter_random_translations()):
            self.assertArraysAlmostEqual((t1*t2).matrix, np.dot(t1.matrix,t2.matrix))
            self.assertArraysAlmostEqual((t1*t2).t, (t2*t1).t)
            self.assertArraysAlmostEqual((t1.inv*t2.inv).t, (t2*t1).inv.t)
            self.assertArraysAlmostEqual((t1*t2).t, t1.t+t2.t)
            self.assertArraysAlmostEqual((t1*(t2*t1)).t, ((t1*t2)*t1).t)

    def test_translation_apply_to_rotation(self):
        for t, r in zip(self.iter_random_translations(), self.iter_random_rotations()):
            self.assertArraysAlmostEqual((t*r).matrix, np.dot(t.matrix,r.matrix))
            self.assertArraysAlmostEqual((t*r).t, t.t)
            self.assertArraysAlmostEqual((t*r).r, r.r)
            self.assertArraysAlmostEqual((t*r.inv).t, t.t)
            self.assertArraysAlmostEqual((t*r.inv).r, r.r.transpose())
            self.assertArraysAlmostEqual((t.inv*r).t, -t.t)
            self.assertArraysAlmostEqual((t.inv*r).r, r.r)
            self.assertArraysAlmostEqual((t.inv*r.inv).t, -t.t)
            self.assertArraysAlmostEqual((t.inv*r.inv).r, r.r.transpose())

    def test_translation_apply_to_complete(self):
        for t, c in zip(self.iter_random_translations(), self.iter_random_completes()):
            self.assertArraysAlmostEqual((t*c).matrix, np.dot(t.matrix,c.matrix))
            self.assertArraysAlmostEqual((t*c).t, t.t+c.t)
            self.assertArraysAlmostEqual((t*c).r, c.r)
            self.assertArraysAlmostEqual((t*c.inv).t, t.t-np.dot(c.r.transpose(),c.t))
            self.assertArraysAlmostEqual((t*c.inv).r, c.r.transpose())
            self.assertArraysAlmostEqual((t.inv*c).t, -t.t+c.t)
            self.assertArraysAlmostEqual((t.inv*c).r, c.r)
            self.assertArraysAlmostEqual((t.inv*c.inv).t, -t.t-np.dot(c.r.transpose(),c.t))
            self.assertArraysAlmostEqual((t.inv*c.inv).r, c.r.transpose())

    def test_rotation_apply_to_vector(self):
        for r in self.iter_random_rotations():
            x = np.random.uniform(-3,3,3)
            self.assertArraysAlmostEqual(r.inv*(r*x), x)
            self.assertArraysAlmostEqual(r*(r.inv*x), x)
            self.assertArraysAlmostEqual((r*r.inv)*x, x)
            self.assertArraysAlmostEqual((r.inv*r)*x, x)
            self.assertArraysAlmostEqual(r.apply_to(x), np.dot(r.r,x))
            self.assertArraysAlmostEqual(r*x, np.dot(r.r,x))

    def test_rotation_apply_to_vectors(self):
        for r in self.iter_random_rotations():
            x = np.random.uniform(-3,3,(10,3))
            self.assertArraysAlmostEqual(r.inv*(r*x), x)
            self.assertArraysAlmostEqual(r*(r.inv*x), x)
            self.assertArraysAlmostEqual((r*r.inv)*x, x)
            self.assertArraysAlmostEqual((r.inv*r)*x, x)
            self.assertArraysAlmostEqual(r.apply_to(x), np.dot(r.r,x.transpose()).transpose())
            self.assertArraysAlmostEqual(r*x, np.dot(r.r,x.transpose()).transpose())
            self.assertArraysAlmostEqual(r.apply_to(x).transpose(), r.apply_to(x.transpose(), columns=True))

    def test_rotation_apply_to_translation(self):
        for r, t in zip(self.iter_random_rotations(), self.iter_random_translations()):
            self.assertArraysAlmostEqual((r*t).matrix, np.dot(r.matrix,t.matrix))
            self.assertArraysAlmostEqual((r*t).t, np.dot(r.r,t.t))
            self.assertArraysAlmostEqual((r*t).r, r.r)
            self.assertArraysAlmostEqual((r*t.inv).t, -np.dot(r.r,t.t))
            self.assertArraysAlmostEqual((r*t.inv).r, r.r)
            self.assertArraysAlmostEqual((r.inv*t).t, np.dot(r.r.transpose(),t.t))
            self.assertArraysAlmostEqual((r.inv*t).r, r.r.transpose())
            self.assertArraysAlmostEqual((r.inv*t.inv).t, -np.dot(r.r.transpose(),t.t))
            self.assertArraysAlmostEqual((r.inv*t.inv).r, r.r.transpose())

    def test_rotation_apply_to_rotation(self):
        for r1, r2 in zip(self.iter_random_rotations(), self.iter_random_rotations()):
            self.assertArraysAlmostEqual((r1*r2).matrix, np.dot(r1.matrix,r2.matrix))
            self.assertArraysAlmostEqual((r1.inv*r2.inv).r, (r2*r1).inv.r)
            self.assertArraysAlmostEqual((r1*r2).r, np.dot(r1.r,r2.r))
            self.assertArraysAlmostEqual((r1*(r2*r1)).r, ((r1*r2)*r1).r)

    def test_rotation_apply_to_complete(self):
        for r, c in zip(self.iter_random_rotations(), self.iter_random_completes()):
            self.assertArraysAlmostEqual((r*c).matrix, np.dot(r.matrix,c.matrix))
            self.assertArraysAlmostEqual((r*c).t, np.dot(r.r,c.t))
            self.assertArraysAlmostEqual((r*c).r, np.dot(r.r,c.r))
            self.assertArraysAlmostEqual((r*c.inv).t, -np.dot(r.r,np.dot(c.r.transpose(),c.t)))
            self.assertArraysAlmostEqual((r*c.inv).r, np.dot(r.r,c.r.transpose()))
            self.assertArraysAlmostEqual((r.inv*c).t, np.dot(r.r.transpose(),c.t))
            self.assertArraysAlmostEqual((r.inv*c).r, np.dot(r.r.transpose(),c.r))
            self.assertArraysAlmostEqual((r.inv*c.inv).t, -np.dot(r.r.transpose(),np.dot(c.r.transpose(),c.t)))
            self.assertArraysAlmostEqual((r.inv*c.inv).r, np.dot(r.r.transpose(),c.r.transpose()))

    def test_complete_apply_to_vector(self):
        for c in self.iter_random_completes():
            x = np.random.uniform(-3,3,3)
            self.assertArraysAlmostEqual(c.inv*(c*x), x)
            self.assertArraysAlmostEqual(c*(c.inv*x), x)
            self.assertArraysAlmostEqual((c*c.inv)*x, x)
            self.assertArraysAlmostEqual((c.inv*c)*x, x)
            self.assertArraysAlmostEqual(c.apply_to(x), np.dot(c.r,x) + c.t)
            self.assertArraysAlmostEqual(c*x, np.dot(c.r,x) + c.t)

    def test_complete_apply_to_vectors(self):
        for c in self.iter_random_completes():
            x = np.random.uniform(-3,3,(10,3))
            self.assertArraysAlmostEqual(c.inv*(c*x), x)
            self.assertArraysAlmostEqual(c*(c.inv*x), x)
            self.assertArraysAlmostEqual((c*c.inv)*x, x)
            self.assertArraysAlmostEqual((c.inv*c)*x, x)
            self.assertArraysAlmostEqual(c.apply_to(x), np.dot(c.r,x.transpose()).transpose() + c.t)
            self.assertArraysAlmostEqual(c*x, np.dot(c.r,x.transpose()).transpose() + c.t)
            self.assertArraysAlmostEqual(c.apply_to(x).transpose(), c.apply_to(x.transpose(), columns=True))

    def test_complete_apply_to_translation(self):
        for c, t in zip(self.iter_random_completes(), self.iter_random_translations()):
            self.assertArraysAlmostEqual((c*t).matrix, np.dot(c.matrix,t.matrix))
            self.assertArraysAlmostEqual((c*t).t, np.dot(c.r,t.t) + c.t)
            self.assertArraysAlmostEqual((c*t).r, c.r)
            self.assertArraysAlmostEqual((c*t.inv).t, -np.dot(c.r,t.t) + c.t)
            self.assertArraysAlmostEqual((c*t.inv).r, c.r)
            self.assertArraysAlmostEqual((c.inv*t).t, np.dot(c.r.transpose(),t.t-c.t))
            self.assertArraysAlmostEqual((c.inv*t).r, c.r.transpose())
            self.assertArraysAlmostEqual((c.inv*t.inv).t, -np.dot(c.r.transpose(),t.t+c.t))
            self.assertArraysAlmostEqual((c.inv*t.inv).r, c.r.transpose())

    def test_complete_apply_to_rotation(self):
        for c, r in zip(self.iter_random_completes(), self.iter_random_rotations()):
            self.assertArraysAlmostEqual((c*r).matrix, np.dot(c.matrix,r.matrix))
            self.assertArraysAlmostEqual((c*r).t, c.t)
            self.assertArraysAlmostEqual((c*r).r, np.dot(c.r,r.r))
            self.assertArraysAlmostEqual((c*r.inv).t, c.t)
            self.assertArraysAlmostEqual((c*r.inv).r, np.dot(c.r,r.r.transpose()))
            self.assertArraysAlmostEqual((c.inv*r).t, -np.dot(c.r.transpose(),c.t))
            self.assertArraysAlmostEqual((c.inv*r).r, np.dot(c.r.transpose(),r.r))
            self.assertArraysAlmostEqual((c.inv*r.inv).t, -np.dot(c.r.transpose(),c.t))
            self.assertArraysAlmostEqual((c.inv*r.inv).r, np.dot(c.r.transpose(),r.r.transpose()))

    def test_complete_apply_to_complete(self):
        for c1, c2 in zip(self.iter_random_completes(), self.iter_random_completes()):
            self.assertArraysAlmostEqual((c1*c2).matrix, np.dot(c1.matrix,c2.matrix))
            self.assertArraysAlmostEqual((c1.inv*c2.inv).t, (c2*c1).inv.t)
            self.assertArraysAlmostEqual((c1.inv*c2.inv).r, (c2*c1).inv.r)
            self.assertArraysAlmostEqual((c1*(c2*c1)).t, ((c1*c2)*c1).t)
            self.assertArraysAlmostEqual((c1*(c2*c1)).r, ((c1*c2)*c1).r)
            self.assertArraysAlmostEqual((c1*c2).t, c1.t + np.dot(c1.r,c2.t))
            self.assertArraysAlmostEqual((c1*c2).r, np.dot(c1.r,c2.r))
            self.assertArraysAlmostEqual((c1*c2.inv).t, c1.t - np.dot(c1.r,np.dot(c2.r.transpose(),c2.t)))
            self.assertArraysAlmostEqual((c1*c2.inv).r, np.dot(c1.r,c2.r.transpose()))
            self.assertArraysAlmostEqual((c1.inv*c2).t, np.dot(c1.r.transpose(),c2.t-c1.t))
            self.assertArraysAlmostEqual((c1.inv*c2).r, np.dot(c1.r.transpose(),c2.r))
            self.assertArraysAlmostEqual((c1.inv*c2.inv).t, np.dot(c1.r.transpose(),-np.dot(c2.r.transpose(),c2.t)-c1.t))
            self.assertArraysAlmostEqual((c1.inv*c2.inv).r, np.dot(c1.r.transpose(),c2.r.transpose()))

    def test_inv_translation(self):
        for t in self.iter_random_translations():
            self.assertArraysAlmostEqual(t.inv.matrix, np.linalg.inv(t.matrix))
            self.assertArraysAlmostEqual((t*t.inv).t, np.zeros(3,float), doabs=True)
            self.assertArraysAlmostEqual((t.inv*t).t, np.zeros(3,float), doabs=True)
            self.assertTrue((t*t.inv).compare(Translation.identity()))
            self.assertTrue((t.inv*t).compare(Translation.identity()))
        self.assertEqual(id(t), id(t.inv.inv))

    def test_inv_rotation(self):
        for r in self.iter_random_rotations():
            self.assertArraysAlmostEqual(r.inv.matrix, np.linalg.inv(r.matrix))
            self.assertArraysAlmostEqual((r*r.inv).r, np.identity(3,float))
            self.assertArraysAlmostEqual((r.inv*r).r, np.identity(3,float))
            self.assertTrue((r*r.inv).compare(Rotation.identity()))
            self.assertTrue((r.inv*r).compare(Rotation.identity()))
        self.assertEqual(id(r), id(r.inv.inv))

    def test_inv_complete(self):
        for c in self.iter_random_completes():
            self.assertArraysAlmostEqual(c.inv.matrix, np.linalg.inv(c.matrix))
            self.assertArraysAlmostEqual((c*c.inv).t, np.zeros(3,float), doabs=True)
            self.assertArraysAlmostEqual((c.inv*c).t, np.zeros(3,float), doabs=True)
            self.assertArraysAlmostEqual((c*c.inv).r, np.identity(3,float))
            self.assertArraysAlmostEqual((c.inv*c).r, np.identity(3,float))
            self.assertTrue((c*c.inv).compare(Complete.identity()))
            self.assertTrue((c.inv*c).compare(Complete.identity()))
        self.assertEqual(id(c), id(c.inv.inv))

    def test_translation_from_matrix(self):
        for t in self.iter_random_translations():
            self.assert_(t.compare(Translation.from_matrix(t.matrix)))

    def test_rotation_from_matrix(self):
        for r in self.iter_random_rotations():
            self.assert_(r.compare(Rotation.from_matrix(r.matrix)))

    def test_complete_from_matrix(self):
        for c in self.iter_random_completes():
            self.assert_(c.compare(Complete.from_matrix(c.matrix)))

    def test_rotation_from_properties(self):
        for r in self.iter_random_rotations():
            self.assert_(r.compare(Rotation.from_properties(*r.properties)))

    def test_complete_from_properties(self):
        for c in self.iter_random_completes():
            self.assert_(c.compare(Complete.from_properties(*c.properties)))

    def test_compare_translations(self):
        for t1, t2 in zip(self.iter_random_translations(), self.iter_random_translations()):
            self.assertFalse(t1.compare(t2))

    def test_compare_rotations(self):
        for r1, r2 in zip(self.iter_random_rotations(), self.iter_random_rotations()):
            self.assertFalse(r1.compare(r2))

    def test_compare_completes(self):
        for c1, c2 in zip(self.iter_random_completes(), self.iter_random_completes()):
            self.assertFalse(c1.compare(c2))

    def test_cast_translation(self):
        for t in self.iter_random_translations():
            self.assertArraysAlmostEqual(t.matrix, Complete.cast(t).matrix)

    def test_cast_rotation(self):
        for r in self.iter_random_rotations():
            self.assertArraysAlmostEqual(r.matrix, Complete.cast(r).matrix)

    def test_cast_complete(self):
        for c in self.iter_random_completes():
            self.assertArraysAlmostEqual(c.matrix, Complete.cast(c).matrix)


    def test_rotation_about_axis(self):
        x = np.array([2,0,0])
        c = Complete.about_axis(np.array([1,0,0]), np.pi, np.array([1,1,0]), False)
        self.assertArraysAlmostEqual(c*x, np.array([1,1,0]))

    def test_superpose(self):
        # create a few test sets with random data points, including degenerate
        # situations. (e.g. one point, two points, linear points)
        references = [ # a list of 2-tuples: (points, degenerate)
            (np.random.normal(0, 5, (n, 3)), False) for n in range(4, 40)
        ] + [
            (np.array([[0,0,1]], float), True),
            (np.array([[0,0,0],[0,0,1]], float), True),
            (np.array([[0,0,0],[0,0,1],[0,0,2]], float), True),
            (np.array([[0,0,0],[0,0,1],[0,0,2],[0,0,4]], float), True),
            (np.random.normal(0, 5, (3, 3)), True)
        ]

        # Do a random transformation on the points
        randomized = []
        for points, degenerate in references:
            #points[:] -= points.mean(axis=0)
            axis = random_unit()
            angle = np.random.uniform(0, np.pi*2)
            transformation = Complete.from_properties(angle, axis, False, np.random.normal(0, 5, 3))
            randomized.append((transformation*points, transformation))

        for (ref_points, degenerate), (tr_points, transf) in zip(references, randomized):
            check_transf = superpose(ref_points, tr_points)
            check_transf2 = fit_rmsd(ref_points, tr_points)[0]
            self.assertArraysAlmostEqual(check_transf.r, check_transf2.r)
            self.assertArraysAlmostEqual(check_transf.t, check_transf2.t)
            # check whether the rotation matrix is orthogonal
            self.assertArraysAlmostEqual(np.dot(check_transf.r, check_transf.r.transpose()), np.identity(3, float), 1e-5)
            # first check whether check_transf brings the tr_points back to the ref_points
            check_points = np.dot(tr_points, check_transf.r.transpose()) + check_transf.t
            self.assertArraysAlmostEqual(ref_points, check_points, 1e-5, doabs=True)
            if not degenerate:
                # if the set of points is not degenerate, check whether transf and check_transf
                # are each other inverses
                self.assert_((transf*check_transf).compare(Complete.identity()))


        # Add some distortion to the transformed points
        for tr_points, transf in randomized:
            tr_points[:] += np.random.normal(0, 1.0, tr_points.shape)

        # Test wether a small modification to the optimal transformation leads
        # to a high rmsd.
        for (ref_points, degenerate), (tr_points, transf) in zip(references, randomized):
            transf = superpose(ref_points, tr_points)
            lowest_rmsd = compute_rmsd(ref_points, transf*tr_points)
            for i in range(10):
                transf_small = Complete.from_properties(0.01, random_unit(), True, np.random.normal(0,0.001,3))
                transf_bis = transf*transf_small
                higher_rmsd = compute_rmsd(ref_points, transf_bis*tr_points)
                self.assert_(lowest_rmsd < higher_rmsd)

    def test_copy(self):
        import copy
        t = Translation([0.8, 3, -0.1])
        self.assert_(t is copy.copy(t))
        self.assert_(t.t.ctypes.data == copy.copy(t).t.ctypes.data)
        self.assert_(t.t.ctypes.data == copy.deepcopy(t).t.ctypes.data)
        r = Rotation.from_properties(0.1, [3.0, 0.0, 1.0], False)
        self.assert_(r is copy.copy(r))
        self.assert_(r.r.ctypes.data == copy.copy(r).r.ctypes.data)
        self.assert_(r.r.ctypes.data == copy.deepcopy(r).r.ctypes.data)
        c = Complete.from_properties(0.1, [3.0, 0.0, 1.0], False, [0.8, 3, -0.1])
        self.assert_(c is copy.copy(c))
        self.assert_(c.r.ctypes.data == copy.copy(c).r.ctypes.data)
        self.assert_(c.r.ctypes.data == copy.deepcopy(c).r.ctypes.data)
        self.assert_(c.t.ctypes.data == copy.copy(c).t.ctypes.data)
        self.assert_(c.t.ctypes.data == copy.deepcopy(c).t.ctypes.data)

    def test_fit_rmsd(self):
        a = np.random.normal(0,1,(20,3))
        trans, a_trans, rmsd = fit_rmsd(a,a)
        self.assertArraysAlmostEqual(trans.r, np.identity(3, float))
        self.assertArraysAlmostEqual(trans.t, np.zeros(3, float), doabs=True)
        self.assertArraysAlmostEqual(a, a_trans)
        self.assertAlmostEqual(rmsd, 0.0)
