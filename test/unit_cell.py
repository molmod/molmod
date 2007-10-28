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


from molmod.unit_cell import UnitCell

import numpy

import unittest, math


__all__ = ["UnitCellTestCase"]


class UnitCellTestCase(unittest.TestCase):
    def test_set(self):
        #print
        for counter in xrange(100):
            #print counter
            uc = UnitCell()
            uc.cell = numpy.random.uniform(-1, 1, (3, 3))
            in_lengths = numpy.random.uniform(0.5, 1, (3,))
            in_angles = numpy.random.uniform(0.3, math.pi/2, (3,))
            #print
            #print " === IN === "
            #print in_lengths
            #print in_angles
            det_before = numpy.linalg.det(uc.cell)
            try:
                uc.set_parameters(in_lengths, in_angles)
            except ValueError, e:
                #print
                #print e.__class__
                #print e
                #print in_lengths
                #print in_angles/math.pi*180
                #print "-"*20
                continue
            det_after = numpy.linalg.det(uc.cell)
            #print " === OUT === "
            out_lengths, out_angles = uc.get_parameters()
            #print out_lengths
            #print out_angles
            self.assertAlmostEqual(sum((in_lengths - out_lengths)**2), 0.0, 5, "Lengths mismatch.")
            self.assertAlmostEqual(sum((in_angles - out_angles)**2), 0.0, 5, "Angles mismatch: %s and %s" % (in_angles, out_angles))
            self.assert_(det_before * det_after > 0, "Handedness has changed.")

    def test_add_periodicities(self):
        for counter in xrange(100):
            #print counter
            uc = UnitCell()
            uc.cell = numpy.random.uniform(-1, 1, (3, 3))
            uc.cell_active[:] = False
            uc.add_cell_vector(numpy.random.uniform(-2,2,3))
            uc.add_cell_vector(numpy.random.uniform(-2,2,3))
            uc.add_cell_vector(numpy.random.uniform(-2,2,3))

    def test_radius_indexes(self):
        cell = numpy.array([
            [1.5, 0, 0],
            [0, 1.5, 0],
            [0, 0, 1.5],
        ], float)
        uc = UnitCell(cell, numpy.array([True, True, True]))
        self.assert_((uc.get_radius_ranges(3.1)==3).all())

        uc.cell = numpy.random.uniform(-1, 1, (3, 3))
        radius = 5.2
        for n in uc.get_radius_indexes(radius):
            self.assert_(numpy.linalg.norm(numpy.dot(uc.cell, n)) < radius)

