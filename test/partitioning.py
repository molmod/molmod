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
        def test_mfn(mfn):
            for i in xrange(10):
                distance = numpy.random.uniform(1, 10, 1)[0]
                counter = 0
                multipoles = numpy.zeros((10, 5), float)
                while counter < 10:
                    delta = numpy.random.uniform(-1, 1, 3)
                    length = numpy.linalg.norm(delta)
                    if length < 1e-4 or length > 1:
                        continue
                    delta *= distance/length
                    tmp = mfn(delta[0],delta[1],delta[2])
                    multipoles[counter,0] = tmp[0]
                    multipoles[counter,1] = numpy.linalg.norm(tmp[1:4])
                    multipoles[counter,2] = numpy.linalg.norm(tmp[4:9])
                    multipoles[counter,3] = numpy.linalg.norm(tmp[9:16])
                    multipoles[counter,4] = numpy.linalg.norm(tmp[16:25])
                    counter += 1
                self.assertAlmostEqual(((multipoles[:,0]-multipoles[:,0].mean())**2).sum(), 0, 5, "Problem with the monopole")
                self.assertAlmostEqual(((multipoles[:,1]-multipoles[:,1].mean())**2).sum(), 0, 5, "Problem with the dipole")
                self.assertAlmostEqual(((multipoles[:,2]-multipoles[:,2].mean())**2).sum(), 0, 5, "Problem with the quadrupole")
                self.assertAlmostEqual(((multipoles[:,3]-multipoles[:,3].mean())**2).sum(), 0, 5, "Problem with the octopole")
                self.assertAlmostEqual(((multipoles[:,4]-multipoles[:,4].mean())**2).sum(), 0, 5, "Problem with the hexadecapole")

        test_mfn(partitioning.ext_mfn)
        test_mfn(partitioning.int_mfn)

    def test_slow_ext(self):
        sqrt2 = numpy.sqrt(2)
        def slow_mfn(x,y,z):
            r = numpy.sqrt(x*x + y*y + z*z)
            theta = numpy.arccos(z/r)
            phi = numpy.arctan2(y,x)
            s = sph_harmonics(4, phi, theta)
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
                        result[counter] = (-1)**m*factor*row[l+m].real*sqrt2
                        counter += 1
                        result[counter] = (-1)**m*factor*row[l+m].imag*sqrt2
                        counter += 1
            return result

        for i in xrange(10):
            delta = numpy.random.uniform(-1, 1, 3)
            v1 = slow_mfn(delta[0],delta[1],delta[2])
            v2 = partitioning.ext_mfn(delta[0],delta[1],delta[2])
            self.assertAlmostEqual(((v1/v2 - 1)**2).sum(), 0.0, 5, "Wrong multipole function %s" % str(v1/v2))

    def test_slow_int(self):
        sqrt2 = numpy.sqrt(2)
        def slow_mfn(x,y,z):
            r = numpy.sqrt(x*x + y*y + z*z)
            theta = numpy.arccos(z/r)
            phi = numpy.arctan2(y,x)
            s = sph_harmonics(4, phi, theta)
            counter = 0
            result = numpy.zeros(25, float)
            for row in s:
                l = len(row)/2
                factor = r**(-l-1)*numpy.sqrt(4*numpy.pi/(2*l+1))
                for m in xrange(l+1):
                    if m == 0:
                        result[counter] = factor*row[l].real
                        counter += 1
                    else:
                        result[counter] = (-1)**m*factor*row[l+m].real*sqrt2
                        counter += 1
                        result[counter] = (-1)**m*factor*row[l+m].imag*sqrt2
                        counter += 1
            return result

        for i in xrange(10):
            delta = numpy.random.uniform(-1, 1, 3)
            v1 = slow_mfn(delta[0],delta[1],delta[2])
            v2 = partitioning.int_mfn(delta[0],delta[1],delta[2])
            self.assertAlmostEqual(((v1/v2 - 1)**2).sum(), 0.0, 5, "Wrong multipole function %s" % str(v1/v2))

    def test_int_ext(self):
        c1 = numpy.array([0,0,0], float)
        c2 = numpy.array([0,0,10],float)
        d1 = numpy.array([0.1,0.2,0],float)
        d2 = numpy.array([0,0.2,0.1],float)
        q1 = 1
        q2 = -1
#        delta = d2-d1
#        deistance = numpy.linalg.norm(delta)

        print
        print "exact energy:", q1*q2/numpy.linalg.norm((c1+d1)-(c2+d2))

        extm1 = q1*partitioning.ext_mfn(d1[0], d1[1], d1[2])
        extm2 = q2*partitioning.ext_mfn(d2[0], d2[1], d2[2])

        e1 = c2+d2-c1
        intm1 = q2*partitioning.int_mfn(e1[0], e1[1], e1[2])
        e2 = c1+d1-c2
        intm2 = q1*partitioning.int_mfn(e2[0], e2[1], e2[2])

        charge1 = extm1[0]
        charge2 = extm2[0]
        dipole1 = numpy.array([extm1[2], extm1[3], extm1[1]], float)
        dipole2 = numpy.array([extm2[2], extm2[3], extm2[1]], float)
        print "charge at 1:", charge1
        print "charge at 2:", charge2
        print "dipole at 1:", dipole1
        print "dipole at 2:", dipole2

        potential1 = intm1[0]
        potential2 = intm2[0]
        field1 = -numpy.array([intm1[2], intm1[3], intm1[1]], float)
        field2 = -numpy.array([intm2[2], intm2[3], intm2[1]], float)
        print "potential at 1:", potential1
        print "potential at 2:", potential2
        print "field at 1:", field1
        print "field at 2:", field2

        print "energy1:", charge1*potential1 - numpy.dot(dipole1,field1)
        print "energy2:", charge2*potential2 - numpy.dot(dipole2,field2)


        print
        print "by hand"
        # at 1
        delta = c1-(c2+d2)
        distance = numpy.linalg.norm(delta)
        potential1 = q2/distance
        print "potential at 1:", potential1
        # at 2
        delta = c2-(c1+d1)
        distance = numpy.linalg.norm(delta)
        potential1 = q1/distance
        print "potential at 2:", potential2
        # at 1
        delta = c1-(c2+d2)
        distance = numpy.linalg.norm(delta)
        field1 = q2*delta/distance**3
        print "field at 1:", field1
        # at 2
        delta = c2-(c1+d1)
        distance = numpy.linalg.norm(delta)
        field2 = q1*delta/distance**3
        print "field at 2:", field2
