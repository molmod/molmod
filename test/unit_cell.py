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


from molmod.unit_cell import UnitCell
from molmod.vectors import random_normal
from molmod.units import deg

import numpy, unittest


__all__ = ["UnitCellTestCase"]


class UnitCellTestCase(unittest.TestCase):
    def get_random_uc(self, r=3, full=True):
        if full:
            return UnitCell(numpy.random.uniform(-r, r, (3, 3)))
        else:
            return UnitCell(
                numpy.random.uniform(-r, r, (3, 3)),
                numpy.random.randint(1,size=3).astype(bool)
            )

    def test_parameters(self):
        for counter in xrange(100):
            in_lengths = numpy.random.uniform(0.5, 1, (3,))
            in_angles = numpy.random.uniform(0.3, numpy.pi/2, (3,))
            try:
                uc = UnitCell.from_parameters3(in_lengths, in_angles)
            except ValueError, e:
                continue
            out_lengths, out_angles = uc.parameters
            self.assertAlmostEqual(sum((in_lengths - out_lengths)**2), 0.0, 5, "Lengths mismatch.")
            self.assertAlmostEqual(sum((in_angles - out_angles)**2), 0.0, 5, "Angles mismatch: %s and %s" % (in_angles, out_angles))

    def test_reciprocal(self):
        for counter in xrange(100):
            uc = UnitCell(numpy.random.uniform(-1, 1, (3, 3)))
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,0], uc.reciprocal[:,0]), 1.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,0], uc.reciprocal[:,1]), 0.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,0], uc.reciprocal[:,2]), 0.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,1], uc.reciprocal[:,0]), 0.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,1], uc.reciprocal[:,1]), 1.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,1], uc.reciprocal[:,2]), 0.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,2], uc.reciprocal[:,0]), 0.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,2], uc.reciprocal[:,1]), 0.0)
            self.assertAlmostEqual(numpy.dot(uc.matrix[:,2], uc.reciprocal[:,2]), 1.0)

    def test_to_fractional(self):
        for i in xrange(100):
            uc = self.get_random_uc()
            fractional = numpy.random.uniform(-0.5, 0.5, 3)
            cartesian = fractional[0]*uc.matrix[:,0] + fractional[1]*uc.matrix[:,1] + fractional[2]*uc.matrix[:,2]
            fractional_bis = uc.to_fractional(cartesian)
            self.assertAlmostEqual(fractional[0], fractional_bis[0])
            self.assertAlmostEqual(fractional[1], fractional_bis[1])
            self.assertAlmostEqual(fractional[2], fractional_bis[2])
        for i in xrange(100):
            uc = self.get_random_uc()
            cartesian = numpy.random.uniform(-3, 3, (10,3))
            fractional = uc.to_fractional(cartesian)
            for i in xrange(10):
                fractional_bis = uc.to_fractional(cartesian[i])
                self.assertAlmostEqual(fractional[i,0], fractional_bis[0])
                self.assertAlmostEqual(fractional[i,1], fractional_bis[1])
                self.assertAlmostEqual(fractional[i,2], fractional_bis[2])

    def test_to_cartesian(self):
        for i in xrange(100):
            uc = self.get_random_uc()
            cartesian = numpy.random.uniform(-3, 3, 3)
            fractional = cartesian[0]*uc.reciprocal[0] + cartesian[1]*uc.reciprocal[1] + cartesian[2]*uc.reciprocal[2]
            cartesian_bis = uc.to_cartesian(fractional)
            self.assertAlmostEqual(cartesian[0], cartesian_bis[0])
            self.assertAlmostEqual(cartesian[1], cartesian_bis[1])
            self.assertAlmostEqual(cartesian[2], cartesian_bis[2])
        for i in xrange(100):
            uc = self.get_random_uc()
            fractional = numpy.random.uniform(-0.5, 0.5, (10,3))
            cartesian = uc.to_cartesian(fractional)
            for i in xrange(10):
                cartesian_bis = uc.to_cartesian(fractional[i])
                self.assertAlmostEqual(cartesian[i,0], cartesian_bis[0])
                self.assertAlmostEqual(cartesian[i,1], cartesian_bis[1])
                self.assertAlmostEqual(cartesian[i,2], cartesian_bis[2])

    def test_consistency(self):
        for i in xrange(100):
            uc = self.get_random_uc()
            cartesian = numpy.random.uniform(-3, 3, 3)
            fractional = uc.to_fractional(cartesian)
            cartesian_bis = uc.to_cartesian(fractional)
            self.assertAlmostEqual(cartesian[0], cartesian_bis[0])
            self.assertAlmostEqual(cartesian[1], cartesian_bis[1])
            self.assertAlmostEqual(cartesian[2], cartesian_bis[2])
        for i in xrange(100):
            uc = self.get_random_uc()
            fractional = numpy.random.uniform(-0.5, 0.5, (10,3))
            cartesian = uc.to_cartesian(fractional)
            fractional_bis = uc.to_fractional(cartesian)
            for i in xrange(10):
                self.assertAlmostEqual(fractional[i,0], fractional_bis[i,0])
                self.assertAlmostEqual(fractional[i,1], fractional_bis[i,1])
                self.assertAlmostEqual(fractional[i,2], fractional_bis[i,2])

    def test_add_periodicities(self):
        for counter in xrange(100):
            uc0 = UnitCell(active=numpy.zeros(3,bool))
            uc1 = uc0.add_cell_vector(numpy.random.uniform(-2,2,3))
            uc2 = uc1.add_cell_vector(numpy.random.uniform(-2,2,3))
            uc3 = uc2.add_cell_vector(numpy.random.uniform(-2,2,3))

    def test_shortest_vector(self):
        # simple case
        uc = UnitCell(numpy.identity(3,float)*3)
        self.assertAlmostEqual(uc.spacings[0], 3.0)
        self.assertAlmostEqual(uc.spacings[1], 3.0)
        self.assertAlmostEqual(uc.spacings[2], 3.0)
        self.assert_(abs(uc.shortest_vector([3, 0, 1]) - numpy.array([0, 0, 1])).max() < 1e-10)
        self.assert_(abs(uc.shortest_vector([-3, 0, 1]) - numpy.array([0, 0, 1])).max() < 1e-10)
        self.assert_(abs(uc.shortest_vector([-2, 0, 1]) - numpy.array([1, 0, 1])).max() < 1e-10)
        self.assert_(abs(uc.shortest_vector([-1.6, 1, 1]) - numpy.array([1.4, 1, 1])).max() < 1e-10)
        self.assert_(abs(uc.shortest_vector([-1.4, 1, 1]) - numpy.array([-1.4, 1, 1])).max() < 1e-10)
        # random tests
        for uc_counter in xrange(10):
            uc = UnitCell(numpy.random.uniform(-1, 1, (3, 3)))
            for r_counter in xrange(10):
                r0 = numpy.random.normal(0, 10, 3)
                r1 = uc.shortest_vector(r0)
                self.assert_(numpy.linalg.norm(r0) >= numpy.linalg.norm(r1))
                index = uc.to_fractional(r0-r1)
                self.assertAlmostEqual(index[0], round(index[0]))
                self.assertAlmostEqual(index[1], round(index[1]))
                self.assertAlmostEqual(index[2], round(index[2]))
            r0 = numpy.random.normal(0, 10, (10,3))
            r1 = uc.shortest_vector(r0)
            for i in xrange(10):
                r1_row_bis = uc.shortest_vector(r0[i])
                self.assertAlmostEqual(r1_row_bis[0], r1[i,0])
                self.assertAlmostEqual(r1_row_bis[1], r1[i,1])
                self.assertAlmostEqual(r1_row_bis[2], r1[i,2])

    def test_spacings(self):
        uc = UnitCell(numpy.identity(3,float)*3)
        self.assertAlmostEqual(uc.spacings[0], 3.0)
        self.assertAlmostEqual(uc.spacings[1], 3.0)
        self.assertAlmostEqual(uc.spacings[2], 3.0)
        for i in xrange(100):
            uc = self.get_random_uc()
            a = uc.matrix[:,0]
            b = uc.matrix[:,1]
            c = uc.matrix[:,2]

            ap = numpy.cross(b,c)
            ap /= numpy.linalg.norm(ap)
            spacing = abs(numpy.dot(a, ap))
            self.assertAlmostEqual(uc.spacings[0], spacing)

    def test_radius_ranges(self):
        for i in xrange(20):
            uc = self.get_random_uc()
            lengths, angles = uc.parameters
            radius = numpy.random.uniform(1,5)
            ranges = uc.get_radius_ranges(radius)
            for j in xrange(100):
                c0 = uc.to_cartesian(numpy.random.uniform(-0.5, 0.5, 3))
                c1 = c0 + radius*random_normal()
                f1 = uc.to_fractional(c1)
                self.assert_((abs(f1) <= ranges+0.5).all(), "f1=%s  ranges=%s" % (f1, ranges))

    def test_div(self):
        for i in xrange(20):
            uc0 = self.get_random_uc(full=False)
            uc1 = uc0/4
            self.assert_(abs(uc0.matrix/4 - uc1.matrix).max() < 1e-10)
            self.assert_((uc0.active == uc1.active).all())

    def test_generalized_volume(self):
        matrix = numpy.array([[10,0,0],[0,10,0],[0,0,10]], float)
        self.assertEqual(UnitCell(matrix, numpy.array([False,False,False])).generalized_volume, -1)
        self.assertEqual(UnitCell(matrix, numpy.array([True,False,False])).generalized_volume, 10)
        self.assertEqual(UnitCell(matrix, numpy.array([True,True,False])).generalized_volume, 100)
        self.assertEqual(UnitCell(matrix, numpy.array([True,True,True])).generalized_volume, 1000)

    def test_alignment(self):
        matrix = numpy.array([[10,10,0],[-10,10,0],[0,0,10]], float)
        uc = UnitCell(matrix)
        r = uc.alignment
        sh = numpy.sqrt(0.5)
        self.assertAlmostEqual(r[0,0], sh)
        self.assertAlmostEqual(r[0,1], -sh)
        self.assertAlmostEqual(r[1,0], sh)
        self.assertAlmostEqual(r[1,1], sh)
        self.assertAlmostEqual(r[2,0], 0)
        self.assertAlmostEqual(r[2,1], 0)
        self.assertAlmostEqual(r[0,2], 0)
        self.assertAlmostEqual(r[1,2], 0)
        self.assertAlmostEqual(r[2,2], 1)

