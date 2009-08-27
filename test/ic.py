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

import molmod.ic as ic
from molmod.molecular_graphs import MolecularGraph, atom_criteria
from molmod.molecules import Molecule
from molmod.graphs import CriteriaSet
from molmod.units import angstrom, deg
from molmod.io.xyz import XYZFile

import unittest, numpy, sys


__all__ = ["ICTestCase"]

num = 100
big = 5
small = 1
eps = 1e-5

def yield_bonds():
    counter = 0
    while counter < num:
        p1 = numpy.random.normal(0,big,3)
        p2 = numpy.random.normal(0,big,3)
        if numpy.linalg.norm(p1-p2) < small:
            continue
        yield p1, p2
        counter += 1

def yield_bends():
    counter = 0
    while counter < num:
        p1 = numpy.random.normal(0,big,3)
        p2 = numpy.random.normal(0,big,3)
        if numpy.linalg.norm(p1-p2) < small:
            continue
        p3 = numpy.random.normal(0,big,3)
        if numpy.linalg.norm(p1-p3) < small:
            continue
        if numpy.linalg.norm(p2-p3) < small:
            continue
        yield p1, p2, p3
        counter += 1

def yield_diheds():
    counter = 0
    while counter < num:
        p1 = numpy.random.normal(0,big,3)
        p2 = numpy.random.normal(0,big,3)
        if numpy.linalg.norm(p1-p2) < small:
            continue
        d = (p1-p2)/numpy.linalg.norm(p1-p2)
        p3 = numpy.random.normal(0,big,3)
        p3_ortho = p3 - d*numpy.dot(p3-p1,d)
        if numpy.linalg.norm(p3_ortho) < small:
            continue
        p4 = numpy.random.normal(0,big,3)
        p4_ortho = p4 - d*numpy.dot(p4-p2,d)
        if numpy.linalg.norm(p4_ortho) < small:
            continue
        yield p1, p2, p3, p4
        counter += 1


class ICTestCase(BaseTestCase):
    def test_diff_bond(self):
        for a1, a2 in yield_bonds():
            d = numpy.random.normal(0,eps,6) # random displacements
            b1 = a1+d[:3] # displaced
            b2 = a2+d[3:] # displaced
            qa, ga, ha = ic.bond_length(a1,a2,2)
            qb, gb, hb = ic.bond_length(b1,b2,2)
            delta1 = qb - qa
            delta2 = numpy.dot(0.5*(ga+gb).ravel(), d)
            error = abs(delta1-delta2)
            oom = abs(delta1)
            self.assert_(error*1e5 < oom)
            delta1 = (gb - ga).ravel()
            delta2 = numpy.dot(0.5*(ha+hb).reshape((6,6)), d)
            error = ((delta1-delta2)**2).mean()
            oom = (delta1**2).mean()
            self.assert_(error*1e5 < oom)

    def test_diff_dot(self):
        def my_dot(x,y,deriv):
            x = ic.Vector3(6,deriv,x,(0,1,2))
            y = ic.Vector3(6,deriv,y,(3,4,5))
            return ic.dot(x,y).results()
        for a1, a2 in yield_bonds():
            d = numpy.random.normal(0,eps,6) # random displacements
            b1 = a1+d[:3] # displaced
            b2 = a2+d[3:] # displaced
            qa, ga, ha = my_dot(a1,a2,2)
            qb, gb, hb = my_dot(b1,b2,2)
            delta1 = qb - qa
            delta2 = numpy.dot(0.5*(ga+gb).ravel(), d)
            error = abs(delta1-delta2)
            oom = abs(delta1)
            self.assert_(error*1e5 < oom)
            delta1 = (gb - ga).ravel()
            delta2 = numpy.dot(0.5*(ha+hb).reshape((6,6)), d)
            error = ((delta1-delta2)**2).mean()
            oom = (delta1**2).mean()
            self.assert_(error*1e5 < oom)

    def test_diff_idiv(self):
        def my_div(t,deriv):
            t0 = ic.Scalar(2,deriv,t[0],0)
            t1 = ic.Scalar(2,deriv,t[1],1)
            t0 /= t1
            return t0.results()
        for counter in xrange(num):
            while True:
                a = numpy.random.normal(0,big,2)
                if a[-1] > small: break
            d = numpy.random.normal(0,eps,2) # random displacements
            b = a+d
            comp = numpy.random.randint(3)
            qa, ga, ha = my_div(a,2)
            qb, gb, hb = my_div(b,2)
            delta1 = qb - qa
            delta2 = numpy.dot(0.5*(ga+gb), d)
            error = abs(delta1-delta2)
            oom = abs(delta1)
            self.assert_(error*1e5 < oom)
            delta1 = (gb - ga)
            delta2 = numpy.dot(0.5*(ha+hb), d)
            error = ((delta1-delta2)**2).mean()
            oom = (delta1**2).mean()
            self.assert_(error*1e5 < oom)

    def test_diff_bend(self):
        for a1, a2, a3 in yield_bends():
            d = numpy.random.normal(0,eps,9)
            b1 = a1+d[:3]
            b2 = a2+d[3:6]
            b3 = a3+d[6:]
            for fn in [ic.bend_cos, ic.bend_angle]:
                qa, ga, ha = fn(a1, a2, a3, 2)
                qb, gb, hb = fn(b1, b2, b3, 2)
                delta1 = qb - qa
                delta2 = numpy.dot(0.5*(ga+gb).ravel(), d)
                error = abs(delta1-delta2)
                oom = abs(delta1)
                self.assert_(error*1e5 < oom)
                delta1 = (gb - ga).ravel()
                delta2 = numpy.dot(0.5*(ha+hb).reshape((9,9)), d)
                error = ((delta1-delta2)**2).mean()
                oom = (delta1**2).mean()
                self.assert_(error*1e5 < oom)

    def test_diff_imul(self):
        def my_mul(t,deriv):
            t0 = ic.Scalar(2,deriv,t[0],0)
            t1 = ic.Scalar(2,deriv,t[1],1)
            t0 *= t1
            return t0.results()
        for counter in xrange(num):
            while True:
                a = numpy.random.normal(0,big,2)
                if a[-1] > small: break
            d = numpy.random.normal(0,eps,2) # random displacements
            b = a+d
            comp = numpy.random.randint(3)
            qa, ga, ha = my_mul(a,2)
            qb, gb, hb = my_mul(b,2)
            delta1 = qb - qa
            delta2 = numpy.dot(0.5*(ga+gb), d)
            error = abs(delta1-delta2)
            oom = abs(delta1)
            self.assert_(error*1e5 < oom)
            delta1 = (gb - ga)
            delta2 = numpy.dot(0.5*(ha+hb), d)
            error = ((delta1-delta2)**2).mean()
            oom = (delta1**2).mean()
            self.assert_(error*1e5 < oom)

    def test_diff_dihed(self):
        for r1, r2, r3, r4 in yield_diheds():
            d = numpy.random.normal(0,eps,12)
            s1 = r1+d[:3]
            s2 = r2+d[3:6]
            s3 = r3+d[6:9]
            s4 = r4+d[9:]
            for fn in [ic.dihed_cos, ic.dihed_angle]:
                qr, gr, hr = fn(r1, r2, r3, r4, 2)
                qs, gs, hs = fn(s1, s2, s3, s4, 2)
                delta1 = qs - qr
                delta2 = numpy.dot(0.5*(gr+gs).ravel(), d)
                error = abs(delta1-delta2)
                oom = abs(delta1)
                self.assert_(error*1e5 < oom)
                delta1 = (gs - gr).ravel()
                delta2 = numpy.dot(0.5*(hs+hr).reshape((12,12)), d)
                error = ((delta1-delta2)**2).mean()
                oom = (delta1**2).mean()
                self.assert_(error*1e5 < oom)


    def test_dihedral_ethene(self):
        mol = XYZFile("input/ethene.xyz").get_molecule()
        c = mol.coordinates.copy()
        self.assertAlmostEqual(ic.dihed_cos(c[2], c[0], c[3], c[5])[0], 1.0)
        self.assertAlmostEqual(ic.dihed_angle(c[2], c[0], c[3], c[5])[0], 0.0)
        for i in xrange(1000):
            angle = numpy.random.uniform(-numpy.pi, numpy.pi)
            radius = numpy.random.uniform(0, 5*angstrom)
            offset = numpy.random.uniform(0, 5*angstrom)
            c[5] = [
                offset,
                -radius*numpy.cos(angle),
                -radius*numpy.sin(angle),
            ]
            self.assertAlmostEqual(ic.dihed_cos(c[2], c[0], c[3], c[5])[0], numpy.cos(angle))
            self.assertAlmostEqual(ic.dihed_angle(c[2], c[0], c[3], c[5])[0], angle)


