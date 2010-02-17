# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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


from molmod.volume import *
from molmod.molecules import Molecule
from molmod.units import angstrom
from molmod.transformations import Complete
from molmod.periodic import periodic


import unittest, numpy


__all__ = ["VolumeTestCase"]


class VolumeTestCase(unittest.TestCase):
    def test_volume_tpa(self):
        mol = Molecule.from_file("input/tpa.xyz")
        vdw_volume, sas_volume, ses_volume, error = estimate_volumes(mol)
        self.assert_(ses_volume > vdw_volume)
        self.assert_(sas_volume > ses_volume)
        #unit = angstrom**3
        #print "vdw_volume [A**3] = %.3f" % (vdw_volume/unit)
        #print "sas_volume [A**3] = %.3f" % (sas_volume/unit)
        #print "ses_volume [A**3] = %.3f" % (ses_volume/unit)
        #print "error [A**3] = %.3f" % (error/unit)

    def test_volume_water(self):
        mol = Molecule.from_file("input/water.xyz")
        vdw_volume, sas_volume, ses_volume, error = estimate_volumes(mol, 1.0)
        self.assert_(ses_volume > vdw_volume)
        self.assert_(sas_volume > ses_volume)

    def test_volume_argon(self):
        mol = Molecule.from_file("input/argon.xyz")
        vdw_volume, sas_volume, ses_volume, error = estimate_volumes(mol)
        self.assertEqual(vdw_volume, ses_volume)
        self.assert_(sas_volume > vdw_volume)

    def test_volume_dinitrogen(self):
        mol = Molecule.from_file("input/dinitrogen.xyz")
        vdw_volume, sas_volume, ses_volume, error = estimate_volumes(mol)
        self.assert_(ses_volume > vdw_volume)
        self.assert_(sas_volume > ses_volume)
        r = periodic["N"].vdw_radius
        ses_upper_limit = 4.0/3.0*numpy.pi*r**3 + numpy.pi*r**2*numpy.linalg.norm(mol.coordinates[0]-mol.coordinates[1])
        self.assert_(vdw_volume < ses_upper_limit)
        self.assert_(ses_volume < ses_upper_limit)

    def test_in_spheres(self):
        from molmodext import in_spheres
        for i in xrange(1000):
            probe = numpy.random.uniform(0, 5, 3)
            probe_radius = numpy.random.uniform(0, 3)
            spheres = numpy.random.uniform(0, 5, (4,3))
            sphere_radii = numpy.random.uniform(0, 3, 4)
            result = in_spheres(probe, probe_radius, spheres, sphere_radii)
            dists = numpy.sqrt(((spheres - probe)**2).sum(axis=1)) - (sphere_radii+probe_radius)
            if (dists > 0).all():
                self.assertEqual(result,0)
            else:
                self.assertEqual(result,-1)

    def test_in_spheres_all(self):
        from molmodext import in_spheres_all
        num = 4
        hits = numpy.zeros(num, numpy.int32)
        for i in xrange(1000):
            probe = numpy.random.uniform(0, 5, 3)
            probe_radius = numpy.random.uniform(0, 3)
            spheres = numpy.random.uniform(0, 5, (num,3))
            sphere_radii = numpy.random.uniform(0, 3, num)
            counter = in_spheres_all(probe, probe_radius, spheres, sphere_radii, hits)
            dists = numpy.sqrt(((spheres - probe)**2).sum(axis=1)) - (sphere_radii+probe_radius)
            self.assertEqual(counter, sum(dists<0))
            self.assertEqual(set(hits[:counter]), set((dists<0).nonzero()[0]))

    def test_center_ses1(self):
        from molmodext import center_ses1
        for i in xrange(1000):
            close1 = numpy.random.uniform(0, 5, 3)
            close1_radius = numpy.random.uniform(0, 5)
            probe = numpy.random.uniform(0, 5, 3)
            center = numpy.zeros(3, float)
            center_ses1(probe, close1, close1_radius, center)
            self.assertAlmostEqual(numpy.linalg.norm(center-close1), close1_radius)

    def test_center_ses2(self):
        from molmodext import center_ses2
        for i in xrange(1000):
            close1 = numpy.random.uniform(0, 5, 3)
            close1_radius = numpy.random.uniform(0, 5)
            close2 = numpy.random.uniform(0, 5, 3)
            close2_radius = numpy.random.uniform(0, 5)
            probe = numpy.random.uniform(0, 5, 3)
            center = numpy.zeros(3, float)
            result = bool(center_ses2(probe, close1, close1_radius, close2, close2_radius, center))
            check_result = (
                (numpy.linalg.norm(close1 - close2) > close1_radius + close2_radius) or
                (numpy.linalg.norm(close1 - close2) < abs(close1_radius - close2_radius))
            )
            self.assertEqual(result, check_result)
            if result == False:
                # test distances
                self.assertAlmostEqual(numpy.linalg.norm(center-close1), close1_radius)
                self.assertAlmostEqual(numpy.linalg.norm(center-close2), close2_radius)
                # test coplanarity
                mat = numpy.array([
                    probe - center,
                    probe - close1,
                    probe - close2,
                ])
                self.assertAlmostEqual(numpy.linalg.det(mat), 0)
                # define a rotation of 180 degrees about the close1-close2 axis
                rotation = Complete.about_axis(close1, numpy.pi, close2-close1)
                new_center = rotation*center
                # test if the rotated center also satisfies the distances and the coplanarity
                self.assertAlmostEqual(numpy.linalg.norm(new_center-close1), close1_radius)
                self.assertAlmostEqual(numpy.linalg.norm(new_center-close2), close2_radius)
                # test coplanarity
                mat = numpy.array([
                    probe - new_center,
                    probe - close1,
                    probe - close2,
                ])
                self.assertAlmostEqual(numpy.linalg.det(mat), 0)
                # test if the center is closes to the probe than new_center
                self.assert_(numpy.linalg.norm(center - probe) < numpy.linalg.norm(new_center - probe))

    def test_center_ses3(self):
        from molmodext import center_ses3
        for i in xrange(1000):
            close1 = numpy.random.uniform(0, 5, 3)
            close1_radius = numpy.random.uniform(0, 5)
            close2 = numpy.random.uniform(0, 5, 3)
            close2_radius = numpy.random.uniform(0, 5)
            close3 = numpy.random.uniform(0, 5, 3)
            close3_radius = numpy.random.uniform(0, 5)
            probe = numpy.random.uniform(0, 5, 3)
            #a = 0.71
            #close1 = numpy.array([1.0, 0.0, 0.0])
            #close1_radius = a
            #close2 = numpy.array([0.0, 1.0, 0.0])
            #close2_radius = a
            #close3 = numpy.array([0.0, 0.0, 1.0])
            #close3_radius = a
            #probe = numpy.array([0.0, 0.0, 0.0])

            center = numpy.zeros(3, float)
            result = bool(center_ses3(probe, close1, close1_radius, close2, close2_radius, close3, close3_radius, center))
            check_result = ( # necessary test, not sufficient
                (numpy.linalg.norm(close1 - close2) > close1_radius + close2_radius) or
                (numpy.linalg.norm(close2 - close3) > close2_radius + close3_radius) or
                (numpy.linalg.norm(close3 - close1) > close3_radius + close1_radius) or
                (numpy.linalg.norm(close1 - close2) < abs(close1_radius - close2_radius)) or
                (numpy.linalg.norm(close2 - close3) < abs(close2_radius - close3_radius)) or
                (numpy.linalg.norm(close3 - close1) < abs(close3_radius - close1_radius))
            )
            if check_result:
                self.assertEqual(result, True)
            if result == False:
                self.assertAlmostEqual(numpy.linalg.norm(center-close1), close1_radius)
                self.assertAlmostEqual(numpy.linalg.norm(center-close2), close2_radius)
                self.assertAlmostEqual(numpy.linalg.norm(center-close3), close3_radius)
                # reflect the center with respect to the plane formed by close1, close2 and close3
                normal = numpy.cross(close1-close2,close3-close2)
                normal /= numpy.linalg.norm(normal)
                delta = normal*numpy.dot(close1 - center, normal)
                new_center = center + 2*delta
                self.assertAlmostEqual(numpy.linalg.norm(new_center-close1), close1_radius)
                self.assertAlmostEqual(numpy.linalg.norm(new_center-close2), close2_radius)
                self.assertAlmostEqual(numpy.linalg.norm(new_center-close3), close3_radius)
                # the new_center must not be closer to the probe than the old center
                self.assert_(numpy.linalg.norm(center - probe) < numpy.linalg.norm(new_center - probe))

    def test_radius_water(self):
        mol = Molecule.from_file("input/water.xyz")
        radius, error = estimate_radius(mol)

    def test_volume_error_water(self):
        mol = Molecule.from_file("input/water.xyz")
        vdw_volumes = []
        sas_volumes = []
        ses_volumes = []
        errors = []
        for i in xrange(100):
            vdw_volume, sas_volume, ses_volume, error = estimate_volumes(mol, 2.0*angstrom)
            vdw_volumes.append(vdw_volume)
            sas_volumes.append(sas_volume)
            ses_volumes.append(ses_volume)
            errors.append(error)
        #print "vdw", numpy.mean(vdw_volumes), numpy.std(vdw_volumes), numpy.mean(errors)
        #print "sas", numpy.mean(sas_volumes), numpy.std(sas_volumes), numpy.mean(errors)
        #print "ses", numpy.mean(ses_volumes), numpy.std(ses_volumes), numpy.mean(errors)
        self.assert_(numpy.std(vdw_volumes) < numpy.mean(errors)*1.5)
        self.assert_(numpy.std(sas_volumes) < numpy.mean(errors)*1.5)
        self.assert_(numpy.std(ses_volumes) < numpy.mean(errors)*1.5)

    def check_volume_order(self, mol):
        radius1 = 0.5*angstrom
        radius2 = 1.5*angstrom
        vdw_volume1, sas_volume1, ses_volume1, error1 = estimate_volumes(mol, radius1, 0.5)
        vdw_volume2, sas_volume2, ses_volume2, error2 = estimate_volumes(mol, radius2, 0.5)
        self.assert_(sas_volume1 < sas_volume2)
        self.assert_(ses_volume1 < ses_volume2)

    def test_volume_order_water(self):
        mol = Molecule.from_file("input/water.xyz")
        self.check_volume_order(mol)

    def test_volume_order_dinitrogen(self):
        mol = Molecule.from_file("input/dinitrogen.xyz")
        self.check_volume_order(mol)

    def check_bigbox(self, mol):
        vdw_volume1, sas_volume1, ses_volume1, error1 = estimate_volumes(mol, 1.0, 2.0, bigbox=True)
        vdw_volume2, sas_volume2, ses_volume2, error2 = estimate_volumes(mol, 1.0, 2.0, bigbox=False)
        self.assert_(abs(vdw_volume1 - vdw_volume2) < error1+error2)
        self.assert_(sas_volume1 > ses_volume2)
        self.assert_(abs(ses_volume1 - ses_volume2) < error1+error2)

    def test_bigbox_water(self):
        mol = Molecule.from_file("input/water.xyz")
        self.check_bigbox(mol)

    def test_bigbox_dinitrogen(self):
        mol = Molecule.from_file("input/dinitrogen.xyz")
        self.check_bigbox(mol)

    def test_bigbox_argon(self):
        mol = Molecule.from_file("input/argon.xyz")
        self.check_bigbox(mol)

    def test_bigbox_tpa(self):
        mol = Molecule.from_file("input/tpa.xyz")
        self.check_bigbox(mol)


