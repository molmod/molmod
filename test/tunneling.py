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


from molmod.tunneling import eckart
from molmod.units import kjmol, cm, unified
from molmod.constants import lightspeed
from molmod.data.periodic import periodic

import unittest, numpy


__all__ = ["TunnelingTestCase"]


invcm = lightspeed/cm

class TunnelingTestCase(unittest.TestCase):
    def test_leo(self):
        # Data is taken from:
        #   Molecular Physics, May 2003, Vol 1.1, No 9, p1329-1338
        #   M.L. Coote, M.A. Collins, L. Radom
        T = 298
        forward = numpy.array([71.7, 49.7, 62.0, 44.5,  6.8], float)*kjmol
        reverse = numpy.array([71.7, 86.8, 76.0, 78.2, 26.4], float)*kjmol

        frequencies = numpy.array([
            [3489, 3376, 3592, 3338, 3418],
            [3411, 3315, 3529, 3274, 3376],
            [3398, 3309, 3528, 3257, 3389],
            [2558, 2578, 2611, 2421, 2096],
            [2518, 2537, 2572, 2382, 2070],
            [2530, 2545, 2585, 2397, 2110],
            [1741, 1658, 1760, 1368,  796],
            [1713, 1640, 1739, 1358,  832],
            [1725, 1656, 1770, 1392,  926],
            [1657, 1592, 1673, 1311,  737],
            [1638, 1568, 1649, 1286,  674],
            [1638, 1579, 1659, 1306,  790],
            [1657, 1602, 1697, 1352,  981],
            [1773, 1739, 1824, 1456,  827],
            [2134, 2061, 2197, 1849, 1358],
            [1959, 1865, 2050, 1600, 1221],
            [1941, 1829, 2038, 1569, 1185],
            [2029, 2010, 2073, 1765, 1126],
            [1866, 1859, 1924, 1574, 1046],
        ], float).transpose()*invcm

        eckart_corrections = numpy.array([
            [7.8E+5, 3.2E+4, 4.0E+5, 1.1E+4, 7.2E+0],
            [5.9E+5, 2.8E+4, 3.3E+5, 1.0E+4, 7.1E+0],
            [5.6E+5, 2.7E+4, 3.3E+5, 9.8E+3, 7.1E+0],
            [1.2E+4, 2.6E+3, 9.5E+3, 8.0E+2, 5.0E+0],
            [1.0E+4, 2.2E+3, 7.9E+3, 6.9E+2, 4.9E+0],
            [1.1E+4, 2.3E+3, 8.4E+3, 7.3E+2, 5.0E+0],
            [9.7E+1, 4.1E+1, 9.4E+1, 1.0E+1, 1.8E+0],
            [8.2E+1, 3.7E+1, 8.3E+1, 9.8E+0, 1.8E+0],
            [8.8E+1, 4.0E+1, 9.9E+1, 1.1E+1, 2.0E+0],
            [5.9E+1, 3.0E+1, 5.7E+1, 8.2E+0, 1.7E+0],
            [5.2E+1, 2.7E+1, 5.0E+1, 7.5E+0, 1.5E+0],
            [5.2E+1, 2.8E+1, 5.3E+1, 8.1E+0, 1.8E+0],
            [5.9E+1, 3.1E+1, 6.5E+1, 9.6E+0, 2.0E+0],
            [1.2E+2, 6.0E+1, 1.4E+2, 1.4E+1, 1.8E+0],
            [1.1E+3, 2.9E+2, 1.1E+3, 7.8E+1, 3.1E+0],
            [3.7E+2, 1.1E+2, 5.0E+2, 2.6E+1, 2.8E+0],
            [3.3E+2, 9.3E+1, 4.7E+2, 2.3E+1, 2.7E+0],
            [5.8E+2, 2.2E+2, 5.7E+2, 5.4E+1, 2.5E+0],
            [2.1E+2, 1.1E+2, 2.4E+2, 2.4E+1, 2.3E+0],
        ], float).transpose()

        for Ef, Er, fs, cs in zip(forward, reverse, frequencies, eckart_corrections):
            for f, c in zip(fs, cs):
                our = eckart(T, nu=f, Ef=Ef, Er=Er)[0]
                self.assert_(abs(our-c)/c < 0.09)
                #if not (("%.1e" % our)==("%.1e" % c)):
                #    print "%.1e" % our, "%.1e" % c






