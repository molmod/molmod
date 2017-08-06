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

from builtins import range
import numpy as np
from nose.plugins.skip import SkipTest

from molmod.test.common import BaseTestCase
from molmod import *


__all__ = ["UnitCellTestCase"]


def get_random_uc(scale=2, num_active=3, min_spacing=0.1):
    while True:
        active = np.random.randint(0, 2, 3).astype(bool)
        if active.sum() == num_active:
            break
    while True:
        try:
            result = UnitCell(np.random.uniform(-scale, scale, (3, 3)), active)
            if not active.any() or result.spacings[active].min() >= min_spacing:
                return result
        except ValueError:
            pass


class UnitCellTestCase(BaseTestCase):
    def test_parameters(self):
        for counter in range(100):
            in_lengths = np.random.uniform(0.5, 1, (3,))
            in_angles = np.random.uniform(0.3, np.pi/2, (3,))
            try:
                uc = UnitCell.from_parameters3(in_lengths, in_angles)
            except ValueError as e:
                continue
            out_lengths, out_angles = uc.parameters
            self.assertArraysAlmostEqual(in_lengths, out_lengths)
            self.assertArraysAlmostEqual(in_angles, out_angles)

    def test_reciprocal(self):
        for counter in range(100):
            uc = get_random_uc(num_active=2)
            if uc.active[0]:
                self.assertAlmostEqual(np.dot(uc.matrix[:,0], uc.reciprocal[:,0]), 1.0)
                self.assertAlmostEqual(np.dot(uc.matrix[:,0], uc.reciprocal[:,1]), 0.0)
                self.assertAlmostEqual(np.dot(uc.matrix[:,0], uc.reciprocal[:,2]), 0.0)
            if uc.active[1]:
                self.assertAlmostEqual(np.dot(uc.matrix[:,1], uc.reciprocal[:,0]), 0.0)
                self.assertAlmostEqual(np.dot(uc.matrix[:,1], uc.reciprocal[:,1]), 1.0)
                self.assertAlmostEqual(np.dot(uc.matrix[:,1], uc.reciprocal[:,2]), 0.0)
            if uc.active[2]:
                self.assertAlmostEqual(np.dot(uc.matrix[:,2], uc.reciprocal[:,0]), 0.0)
                self.assertAlmostEqual(np.dot(uc.matrix[:,2], uc.reciprocal[:,1]), 0.0)
                self.assertAlmostEqual(np.dot(uc.matrix[:,2], uc.reciprocal[:,2]), 1.0)

    def test_reciprocal_bis(self):
        for i in range(100):
            uc = get_random_uc(num_active=2)
            active, inactive = uc.active_inactive
            for i in inactive:
                self.assertEqual(abs(uc.reciprocal[:,i]).max(), 0.0)
            if len(active) == 1:
                unit = uc.reciprocal[:,active[0]].copy()
                r = np.linalg.norm(unit)
                unit /= r
                cell_vector = unit/r
                self.assertArraysAlmostEqual(cell_vector, uc.matrix[:,active[0]])
            elif len(active) == 2:
                # construct an auxiliary normal vector
                normal = np.cross(uc.reciprocal[:,active[0]], uc.reciprocal[:,active[1]])
                norm = np.linalg.norm(normal)
                # reconstruct the cell vectors
                cell_vector0 = np.cross(uc.reciprocal[:,active[1]], normal)/norm**2
                cell_vector1 = -np.cross(uc.reciprocal[:,active[0]], normal)/norm**2
                self.assertArraysAlmostEqual(cell_vector0, uc.matrix[:,active[0]])
                self.assertArraysAlmostEqual(cell_vector1, uc.matrix[:,active[1]])
            elif len(active) == 3:
                self.assertArraysEqual(uc.reciprocal, uc.reciprocal)

    def test_to_fractional(self):
        for i in range(100):
            uc = get_random_uc()
            fractional = np.random.uniform(-0.5, 0.5, 3)
            cartesian = fractional[0]*uc.matrix[:,0] + fractional[1]*uc.matrix[:,1] + fractional[2]*uc.matrix[:,2]
            fractional_bis = uc.to_fractional(cartesian)
            self.assertArraysAlmostEqual(fractional, fractional_bis)
        for i in range(100):
            uc = get_random_uc()
            cartesian = np.random.uniform(-3, 3, (10,3))
            fractional = uc.to_fractional(cartesian)
            for i in range(10):
                fractional_bis = uc.to_fractional(cartesian[i])
                self.assertArraysAlmostEqual(fractional[i], fractional_bis)

    def test_to_cartesian(self):
        for i in range(100):
            uc = get_random_uc()
            cartesian = np.random.uniform(-3, 3, 3)
            fractional = cartesian[0]*uc.reciprocal[0] + cartesian[1]*uc.reciprocal[1] + cartesian[2]*uc.reciprocal[2]
            cartesian_bis = uc.to_cartesian(fractional)
            self.assertArraysAlmostEqual(cartesian, cartesian_bis)
        for i in range(100):
            uc = get_random_uc()
            fractional = np.random.uniform(-0.5, 0.5, (10,3))
            cartesian = uc.to_cartesian(fractional)
            for i in range(10):
                cartesian_bis = uc.to_cartesian(fractional[i])
                self.assertArraysAlmostEqual(cartesian[i], cartesian_bis)

    def test_consistency(self):
        for i in range(100):
            uc = get_random_uc()
            cartesian = np.random.uniform(-3, 3, 3)
            fractional = uc.to_fractional(cartesian)
            cartesian_bis = uc.to_cartesian(fractional)
            self.assertArraysAlmostEqual(cartesian, cartesian_bis)
        for i in range(100):
            uc = get_random_uc()
            fractional = np.random.uniform(-0.5, 0.5, (10,3))
            cartesian = uc.to_cartesian(fractional)
            fractional_bis = uc.to_fractional(cartesian)
            self.assertArraysAlmostEqual(fractional, fractional_bis)

    def test_add_periodicities(self):
        for counter in range(100):
            uc0 = UnitCell(np.identity(3, float), np.zeros(3,bool))
            uc1 = uc0.add_cell_vector(np.random.uniform(-2,2,3))
            uc2 = uc1.add_cell_vector(np.random.uniform(-2,2,3))
            uc3 = uc2.add_cell_vector(np.random.uniform(-2,2,3))

    def test_shortest_vector(self):
        raise SkipTest
        # simple cases
        uc = UnitCell(np.identity(3,float)*3)
        self.assertArraysAlmostEqual(uc.shortest_vector([3, 0, 1]), np.array([0, 0, 1]))
        self.assertArraysAlmostEqual(uc.shortest_vector([-3, 0, 1]), np.array([0, 0, 1]))
        self.assertArraysAlmostEqual(uc.shortest_vector([-2, 0, 1]), np.array([1, 0, 1]))
        self.assertArraysAlmostEqual(uc.shortest_vector([-1.6, 1, 1]), np.array([1.4, 1, 1]))
        self.assertArraysAlmostEqual(uc.shortest_vector([-1.4, 1, 1]), np.array([-1.4, 1, 1]))
        # simple cases
        uc = UnitCell(np.identity(3,float)*3, np.array([True, False, False]))
        self.assertArraysAlmostEqual(uc.shortest_vector([3, 0, 1]), np.array([0, 0, 1]))
        self.assertArraysAlmostEqual(uc.shortest_vector([3, 0, 3]), np.array([0, 0, 3]))
        # random tests
        for uc_counter in range(1000):
            uc = get_random_uc(num_active=2)
            for r_counter in range(10):
                r0 = np.random.normal(0, 10, 3)
                r1 = uc.shortest_vector(r0)
                change = r1 - r0
                self.assert_(np.dot(change, r0) <= 0)
                #self.assert_(np.linalg.norm(r0) >= np.linalg.norm(r1))
                index = uc.to_fractional(r0-r1)
                self.assertArraysAlmostEqual(index, np.round(index), doabs=True)
                index = uc.to_fractional(r1)
                self.assert_(index.max()<0.5)
                self.assert_(index.max()>=-0.5)
            r0 = np.random.normal(0, 10, (10,3))
            r1 = uc.shortest_vector(r0)
            for i in range(10):
                r1_row_bis = uc.shortest_vector(r0[i])
                self.assertArraysAlmostEqual(r1_row_bis, r1[i], doabs=True)

    def test_shortest_vector_trivial(self):
        uc = UnitCell(np.identity(3, float))
        half = np.array([0.5,0.5,0.5])
        self.assertArraysEqual(uc.shortest_vector(half), -half)
        self.assertArraysEqual(uc.shortest_vector(-half), -half)

    def test_spacings(self):
        uc = UnitCell(np.identity(3,float)*3)
        self.assertArraysAlmostEqual(uc.spacings, np.ones(3, float)*3.0)
        for i in range(100):
            uc = get_random_uc()
            a = uc.matrix[:,0]
            b = uc.matrix[:,1]
            c = uc.matrix[:,2]

            ap = np.cross(b,c)
            ap /= np.linalg.norm(ap)
            spacing = abs(np.dot(a, ap))
            self.assertAlmostEqual(uc.spacings[0], spacing)

    def test_radius_ranges(self):
        for i in range(20):
            uc = get_random_uc()
            radius = np.random.uniform(1,5)
            ranges = uc.get_radius_ranges(radius)
            for j in range(100):
                c0 = uc.to_cartesian(np.random.uniform(-0.5, 0.5, 3))
                c1 = c0 + radius*random_unit()
                f1 = uc.to_fractional(c1)
                self.assert_((abs(f1) <= ranges+0.5).all(), "f1=%s  ranges=%s" % (f1, ranges))

    def test_radius_ranges_2d(self):
        uc = UnitCell(np.identity(3, float), np.array([True, True, False]))
        self.assertArraysEqual(uc.get_radius_ranges(3), np.array([3,3,0]))

    def test_radius_indexes(self):
        for i in range(20):
            uc = get_random_uc()
            radius = np.random.uniform(1,2)*abs(uc.volume)**(0.333)

            #uc = UnitCell.from_parameters3(
            #    np.array([3.73800243,  2.35503196,  3.25130153]),
            #    np.array([29.78777448, 64.81228452, 86.7641093])*deg,
            #)
            #lengths, angles = uc.parameters
            #radius = 1.98042040465

            #matrix = np.array([[2,0,0],[0,2,0],[0,0,2]], float)
            #uc = UnitCell(matrix)
            #radius = 5.3

            #lengths, angles = uc.parameters
            #print lengths
            #print angles/deg
            #print radius

            ranges = uc.get_radius_ranges(radius)
            #print "ranges", ranges
            if np.product(ranges) > 100:
                continue
            indexes = uc.get_radius_indexes(radius)
            self.assert_(len(indexes)*abs(uc.volume) > 4.0/3.0*np.pi*radius**3)
            #print "testing distances"
            for j in range(20):
                c0 = uc.to_cartesian(np.random.uniform(-0.5, 0.5, 3))
                c1 = uc.to_cartesian(np.random.uniform(-0.5, 0.5, 3))
                relative = c1 - c0
                # compute all distances between c0 and c1 based on radius
                # ranges
                distances_slow = []
                for i0 in range(-ranges[0], ranges[0]+1):
                    for i1 in range(-ranges[1], ranges[1]+1):
                        for i2 in range(-ranges[2], ranges[2]+1):
                            delta = uc.to_cartesian([i0,i1,i2])
                            distance = np.linalg.norm(relative + delta)
                            if distance <= radius:
                                distances_slow.append(distance)
                distances_slow.sort()
                distances_fast = []
                for index in indexes:
                    delta = uc.to_cartesian(index)
                    distance = np.linalg.norm(relative + delta)
                    if distance <= radius:
                        distances_fast.append(distance)
                distances_fast.sort()
                #print distances_slow
                #print distances_fast
                #print
                self.assertArraysAlmostEqual(np.array(distances_slow), np.array(distances_fast))

    def test_radius_indexes_1d(self):
        uc = UnitCell(np.identity(3, float), np.array([True, False, False]))
        indexes = uc.get_radius_indexes(0.5)
        expected_indexes = np.array([
            [-1,  0,  0],
            [ 0,  0,  0],
            [ 1,  0,  0],
        ])
        self.assertArraysEqual(indexes, expected_indexes)
        uc = UnitCell(np.identity(3, float), np.array([True, False, False]))
        indexes = uc.get_radius_indexes(1.8, np.array([4,-1,-1]))
        expected_indexes = np.array([
            [-2,  0,  0],
            [-1,  0,  0],
            [ 0,  0,  0],
            [ 1,  0,  0],
        ])
        self.assertArraysEqual(indexes, expected_indexes)

    def test_radius_indexes_2d(self):
        uc = UnitCell(np.identity(3, float), np.array([True, True, False]))
        indexes = uc.get_radius_indexes(0.5)
        expected_indexes = np.array([
            [-1, -1,  0],
            [-1,  0,  0],
            [-1,  1,  0],
            [ 0, -1,  0],
            [ 0,  0,  0],
            [ 0,  1,  0],
            [ 1, -1,  0],
            [ 1,  0,  0],
            [ 1,  1,  0],
        ])
        self.assertArraysEqual(indexes, expected_indexes)

    def test_radius_indexes_2d_graphical(self):
        raise SkipTest
        #uc = UnitCell(np.array([
        #    [2.0, 1.0, 0.0],
        #    [0.0, 0.2, 0.0],
        #    [0.0, 0.0, 10.0],
        #]))
        #radius = 0.8
        uc = UnitCell(np.array([
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 10.0],
        ]))
        radius = 5.3
        #uc = UnitCell(np.array([
        #    [1.0, 1.0, 0.0],
        #    [0.0, 1.0, 0.0],
        #    [0.0, 0.0, 1.0],
        #]))
        #radius = 0.9

        fracs = np.arange(-0.5, 0.55, 0.1)
        import pylab
        from matplotlib.patches import Circle, Polygon
        from matplotlib.lines import Line2D
        pylab.clf()
        for i0 in fracs:
            for i1 in fracs:
                center = uc.to_cartesian([i0,i1,0.0])
                pylab.gca().add_artist(Circle((center[0], center[1]), radius, fill=True, fc='#7777AA', ec='none'))
        pylab.gca().add_artist(Circle((0, 0), radius, fill=True, fc='#0000AA', ec='none'))

        ranges = uc.get_radius_ranges(radius)
        indexes = uc.get_radius_indexes(radius)

        for i in range(-ranges[0]-1, ranges[0]+1):
            start = uc.to_cartesian([i+0.5, -ranges[1]-0.5, 0])
            end = uc.to_cartesian([i+0.5, ranges[1]+0.5, 0])
            pylab.gca().add_artist(Line2D([start[0], end[0]], [start[1], end[1]], color="k", linewidth=1))

        for i in range(-ranges[1]-1, ranges[1]+1):
            start = uc.to_cartesian([-ranges[0]-0.5, i+0.5, 0])
            end = uc.to_cartesian([ranges[0]+0.5, i+0.5, 0])
            pylab.gca().add_artist(Line2D([start[0], end[0]], [start[1], end[1]], color="k", linewidth=1))

        for i in range(-ranges[0], ranges[0]+1):
            start = uc.to_cartesian([i, -ranges[1]-0.5, 0])
            end = uc.to_cartesian([i, ranges[1]+0.5, 0])
            pylab.gca().add_artist(Line2D([start[0], end[0]], [start[1], end[1]], color="k", linewidth=0.5, linestyle="--"))

        for i in range(-ranges[1], ranges[1]+1):
            start = uc.to_cartesian([-ranges[0]-0.5, i, 0])
            end = uc.to_cartesian([ranges[0]+0.5, i, 0])
            pylab.gca().add_artist(Line2D([start[0], end[0]], [start[1], end[1]], color="k", linewidth=0.5, linestyle="--"))

        for i0,i1,i2 in indexes:
            if i2 != 0:
                continue
            corners = uc.to_cartesian(np.array([
                [i0-0.5, i1-0.5, 0.0],
                [i0-0.5, i1+0.5, 0.0],
                [i0+0.5, i1+0.5, 0.0],
                [i0+0.5, i1-0.5, 0.0],
            ]))
            pylab.gca().add_artist(Polygon(corners[:,:2], fill=True, ec='none', fc='r', alpha=0.5))

        corners = uc.to_cartesian(np.array([
            [-ranges[0]-0.5, -ranges[1]-0.5, 0.0],
            [-ranges[0]-0.5, +ranges[1]+0.5, 0.0],
            [+ranges[0]+0.5, +ranges[1]+0.5, 0.0],
            [+ranges[0]+0.5, -ranges[1]-0.5, 0.0],
        ]))
        pylab.xlim(1.1*corners[:,:2].min(), 1.1*corners[:,:2].max())
        pylab.ylim(1.1*corners[:,:2].min(), 1.1*corners[:,:2].max())
        #pylab.xlim(-1.5*radius, 1.5*radius)
        #pylab.ylim(-1.5*radius, 1.5*radius)
        #pylab.savefig("radius_indexes_2d.png")

    def test_div(self):
        for i in range(20):
            uc0 = get_random_uc(num_active=2)
            x = np.random.uniform(1,2,3)
            uc1 = uc0/x
            self.assertArraysAlmostEqual(uc0.matrix/x, uc1.matrix)
            self.assertArraysEqual(uc0.active, uc1.active)
            self.assertAlmostEqual(
                np.dot(uc0.matrix[:,0], uc1.matrix[:,0])
                /np.linalg.norm(uc0.matrix[:,0])
                /np.linalg.norm(uc1.matrix[:,0]), 1
            )
            self.assertAlmostEqual(
                np.dot(uc0.matrix[:,1], uc1.matrix[:,1])
                /np.linalg.norm(uc0.matrix[:,1])
                /np.linalg.norm(uc1.matrix[:,1]), 1
            )
            self.assertAlmostEqual(
                np.dot(uc0.matrix[:,2], uc1.matrix[:,2])
                /np.linalg.norm(uc0.matrix[:,2])
                /np.linalg.norm(uc1.matrix[:,2]), 1
            )

    def test_mul(self):
        for i in range(20):
            uc0 = get_random_uc(num_active=2)
            x = np.random.uniform(1,2,3)
            uc1 = uc0*x
            self.assertArraysAlmostEqual(uc0.matrix*x, uc1.matrix)
            self.assertArraysEqual(uc0.active, uc1.active)
            self.assertAlmostEqual(
                np.dot(uc0.matrix[:,0], uc1.matrix[:,0])
                /np.linalg.norm(uc0.matrix[:,0])
                /np.linalg.norm(uc1.matrix[:,0]), 1
            )
            self.assertAlmostEqual(
                np.dot(uc0.matrix[:,1], uc1.matrix[:,1])
                /np.linalg.norm(uc0.matrix[:,1])
                /np.linalg.norm(uc1.matrix[:,1]), 1
            )
            self.assertAlmostEqual(
                np.dot(uc0.matrix[:,2], uc1.matrix[:,2])
                /np.linalg.norm(uc0.matrix[:,2])
                /np.linalg.norm(uc1.matrix[:,2]), 1
            )

    def test_volume(self):
        matrix = np.array([[10,0,0],[0,10,0],[0,0,10]], float)
        self.assertAlmostEqual(UnitCell(matrix, np.array([False,False,False])).volume, -1)
        self.assertAlmostEqual(UnitCell(matrix, np.array([True,False,False])).volume, 10)
        self.assertAlmostEqual(UnitCell(matrix, np.array([True,True,False])).volume, 100)
        self.assertAlmostEqual(UnitCell(matrix, np.array([True,True,True])).volume, 1000)

    def test_alignment_a(self):
        matrix = np.array([[10,10,0],[-10,10,0],[0,0,10]], float)
        uc = UnitCell(matrix)
        r = uc.alignment_a
        sh = np.sqrt(0.5)
        expected_r = np.array([
            [sh, -sh, 0],
            [sh, sh, 0],
            [0, 0, 1],
        ], float)
        self.assertArraysAlmostEqual(r.r, expected_r)
        uc = r*uc
        expected_matrix = np.array([
            [10/sh, 0, 0],
            [0, 10/sh, 0],
            [0, 0, 10],
        ])
        self.assertArraysAlmostEqual(uc.matrix, expected_matrix)

    def test_alignment_c(self):
        matrix = np.array([[10,0,0],[0,10,-10],[0,10,10]], float)
        uc = UnitCell(matrix)
        r = uc.alignment_c
        uc = r*uc
        sh = np.sqrt(0.5)
        expected_matrix = np.array([
            [10, 0, 0],
            [0, 10/sh, 0],
            [0, 0, 10/sh],
        ])
        self.assertArraysAlmostEqual(uc.matrix, expected_matrix)

    def test_shortest_vector_aperiodic(self):
        unit_cell = UnitCell(np.identity(3, float), np.zeros(3, bool))
        shortest = unit_cell.shortest_vector(np.ones(3, float))
        expected = np.ones(3, float)
        self.assertArraysAlmostEqual(shortest, expected)

    def test_ordered(self):
        for i in range(10):
            uc0 = get_random_uc(num_active=2)
            uc1 = uc0.ordered
            self.assertAlmostEqual(uc0.volume, uc1.volume)
            self.assertEqual(uc0.active.sum(), uc1.active.sum())
            self.assert_(uc1.active[0] >= uc1.active[1])
            self.assert_(uc1.active[1] >= uc1.active[2])
            for i1, i0 in enumerate(uc0.active_inactive[0]):
                self.assertArraysAlmostEqual(uc0.matrix[:,i0], uc1.matrix[:,i1])
                self.assertEqual(uc0.active[i0], uc1.active[i1])
