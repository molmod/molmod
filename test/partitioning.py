# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from molmod.helpers import partitioning
from molmod.descriptors import sph_harmonics

import unittest, numpy

__all__ = ["TestPartitioning"]


class TestPartitioning(unittest.TestCase):
    def test_rotation_invariance(self):
        for i in xrange(10):
            print " i", i
            distance = numpy.random.uniform(1, 10, 1)[0]
            counter = 0
            multipoles = numpy.zeros((10, 5), float)
            while counter < 10:
                print "  c", counter
                delta = numpy.random.uniform(-1, 1, 3)
                length = numpy.linalg.norm(delta)
                if length < 1e-4 or length > 1:
                    continue
                delta *= distance/length
                tmp = partitioning.mfn(delta[0],delta[1],delta[2])
                multipoles[counter,0] = tmp[0]
                multipoles[counter,1] = numpy.linalg.norm(tmp[1:4])
                multipoles[counter,2] = numpy.linalg.norm(tmp[4:9])
                multipoles[counter,3] = numpy.linalg.norm(tmp[9:16])
                multipoles[counter,4] = numpy.linalg.norm(tmp[16:25])
                counter += 1
            print multipoles
            self.assertAlmostEqual(((multipoles[:,0]-multipoles[:,0].mean())**2).sum(), 0, 5, "Problem with the monopole")
            self.assertAlmostEqual(((multipoles[:,1]-multipoles[:,1].mean())**2).sum(), 0, 5, "Problem with the dipole")
            self.assertAlmostEqual(((multipoles[:,2]-multipoles[:,2].mean())**2).sum(), 0, 5, "Problem with the quadrupole")
            self.assertAlmostEqual(((multipoles[:,3]-multipoles[:,3].mean())**2).sum(), 0, 5, "Problem with the octopole")
            self.assertAlmostEqual(((multipoles[:,4]-multipoles[:,4].mean())**2).sum(), 0, 5, "Problem with the hexadecapole")

    def tst_slow(self):
        def slow_mfn(x,y,z):
            r = numpy.sqrt(x*x + y*y + z*z)
            theta = numpy.arccos(z/r)
            phi = numpy.arctan2(y,z)
            s = sph_harmonics(4, phi, theta)
            print s
            counter = 0
            result = numpy.zeros(25, float)
            for row in s:
                l = len(row)/2
                factor = r**l*numpy.sqrt(4*numpy.pi/(2*l+1))
                for m in xrange(l+1):
                    if m == 0:
                        result[counter] = factor*row[l].real
                        counter += 1
                    else:
                        result[counter] = factor*row[l+m].real
                        counter += 1
                        result[counter] = factor*row[l-m].imag
                        counter += 1
            return result

        for i in xrange(10):
            delta = numpy.random.uniform(-1, 1, 3)
            v1 = slow_mfn(delta[0],delta[1],delta[2])
            v2 = partitioning.mfn(delta[0],delta[1],delta[2])
            print v1
            print v2
            print
            print



