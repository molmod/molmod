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
import pkg_resources

import molmod.ic as ic
from molmod import *


num = 10
big = 5
small = 1
eps = 1e-5


def iter_bonds():
    counter = 0
    while counter < num:
        p1 = np.random.normal(0,big,3)
        p2 = np.random.normal(0,big,3)
        if np.linalg.norm(p1-p2) < small:
            continue
        yield p1, p2
        counter += 1


def iter_bends():
    counter = 0
    while counter < num:
        p1 = np.random.normal(0,big,3)
        p2 = np.random.normal(0,big,3)
        if np.linalg.norm(p1-p2) < small:
            continue
        p3 = np.random.normal(0,big,3)
        if np.linalg.norm(p1-p3) < small:
            continue
        if np.linalg.norm(p2-p3) < small:
            continue
        yield p1, p2, p3
        counter += 1


def iter_diheds():
    counter = 0
    while counter < num:
        p1 = np.random.normal(0,big,3)
        p2 = np.random.normal(0,big,3)
        if np.linalg.norm(p1-p2) < small:
            continue
        d = (p1-p2)/np.linalg.norm(p1-p2)
        p3 = np.random.normal(0,big,3)
        p3_ortho = p3 - d*np.dot(p3-p1,d)
        if np.linalg.norm(p3_ortho) < small:
            continue
        p4 = np.random.normal(0,big,3)
        p4_ortho = p4 - d*np.dot(p4-p2,d)
        if np.linalg.norm(p4_ortho) < small:
            continue
        yield p1, p2, p3, p4
        counter += 1


def iter_pairs():
    for counter in range(num):
        while True:
            x = np.random.normal(0,big,2)
            if x[-1] > small: break
        yield x


def check_diff_ic(icfn, iterp, shape=(-1,3), period=None):
    def fnv(x0, do_gradient=False):
        q, g = icfn(x0.reshape(shape),1)
        if do_gradient:
            return q, g.ravel()
        else:
            return q

    def fng(x0, do_gradient=False):
        q, g, h = icfn(x0.reshape(shape),2)
        if do_gradient:
            return g.ravel(), h.reshape(g.size,g.size)
        else:
            return g.ravel()

    for ps in iterp():
        x0 = np.array(ps).ravel()
        dxs = np.random.normal(0, eps, (100, len(x0)))
        check_delta(fnv, x0, dxs, period)
        check_delta(fng, x0, dxs, period)


def test_diff_bond():
    check_diff_ic(ic.bond_length, iter_bonds)


def test_diff_dot():
    def my_dot(rs,deriv):
        x = ic.Vector3(6,deriv,rs[0],(0,1,2))
        y = ic.Vector3(6,deriv,rs[1],(3,4,5))
        return ic.dot(x,y).results()
    check_diff_ic(my_dot, iter_bonds)


class MyCross(object):
    def __init__(self, index):
        self.index = index

    def __call__(self, rs, deriv):
        x = ic.Vector3(6,deriv,rs[0],(0,1,2))
        y = ic.Vector3(6,deriv,rs[1],(3,4,5))
        if self.index == 0:
            return ic.cross(x,y).x.results()
        elif self.index == 1:
            return ic.cross(x,y).y.results()
        elif self.index == 2:
            return ic.cross(x,y).z.results()
        else:
            raise NotImplementedError


def test_diff_cross():
    check_diff_ic(MyCross(0), iter_bonds)
    check_diff_ic(MyCross(1), iter_bonds)
    check_diff_ic(MyCross(2), iter_bonds)


def test_diff_idiv():
    def my_div(t,deriv):
        t0 = ic.Scalar(2,deriv,t[0],0)
        t1 = ic.Scalar(2,deriv,t[1],1)
        t0 /= t1
        return t0.results()
    check_diff_ic(my_div, iter_pairs, (2,))


def test_diff_bend_angle():
    check_diff_ic(ic.bend_angle, iter_bends)


def test_diff_bend_cos():
    check_diff_ic(ic.bend_cos, iter_bends)


def test_diff_imul():
    def my_mul(t,deriv):
        t0 = ic.Scalar(2,deriv,t[0],0)
        t1 = ic.Scalar(2,deriv,t[1],1)
        t0 *= t1
        return t0.results()
    check_diff_ic(my_mul, iter_pairs, (2,))


def test_diff_dihed_cos():
    check_diff_ic(ic.dihed_cos, iter_diheds)


def test_diff_dihed_angle():
    check_diff_ic(ic.dihed_angle, iter_diheds, period=2*np.pi)


def iter_diheds_special():
    yield np.array([[+1,0,1], [0,0,0], [0,0,1], [1,0,1]])
    yield np.array([[0,+1,1], [0,0,0], [0,0,1], [1,0,1]])
    yield np.array([[-1,0,1], [0,0,0], [0,0,1], [1,0,1]])
    yield np.array([[0,-1,1], [0,0,0], [0,0,1], [1,0,1]])


def test_diff_dihed_angle_special():
    check_diff_ic(ic.dihed_angle, iter_diheds_special, period=2*np.pi)


def test_diff_opbend_cos():
    check_diff_ic(ic.opbend_cos, iter_diheds)


def test_diff_opbend_angle():
    check_diff_ic(ic.opbend_angle, iter_diheds)

def test_diff_opbend_mcos():
    check_diff_ic(ic.opbend_mcos, iter_diheds)

def test_diff_opbend_mangle():
    check_diff_ic(ic.opbend_mangle, iter_diheds)

def test_diff_opbend_dist():
    check_diff_ic(ic.opbend_dist, iter_diheds)


def test_dihedral_ethene():
    mol = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/ethene.xyz"))
    c = mol.coordinates.copy()
    assert abs(ic.dihed_cos([c[2], c[0], c[3], c[5]])[0] - 1.0) < 1e-5
    assert abs(ic.dihed_angle([c[2], c[0], c[3], c[5]])[0]) < 1e-5
    for i in range(1000):
        angle = np.random.uniform(-np.pi, np.pi)
        radius = np.random.uniform(0, 5*angstrom)
        offset = np.random.uniform(0, 5*angstrom)
        c[5] = [
            offset,
            -radius*np.cos(angle),
            -radius*np.sin(angle),
        ]
        assert abs(ic.dihed_cos([c[2], c[0], c[3], c[5]])[0] - np.cos(angle)) < 1e-5
        assert abs(ic.dihed_angle([c[2], c[0], c[3], c[5]])[0] - angle) < 1e-5


def test_opbend_ethene():
    mol = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/ethene.xyz"))
    c = mol.coordinates.copy()
    assert abs(ic.opbend_cos([c[0], c[5], c[4], c[3]])[0] - 1.0) < 1e-5
    assert abs(ic.opbend_angle([c[0], c[5], c[4], c[3]])[0]) < 1e-5
    assert abs(ic.opbend_mcos([c[0], c[5], c[4], c[3]])[0] - 1.0) < 1e-5
    assert abs(ic.opbend_mangle([c[0], c[5], c[4], c[3]])[0]) < 1e-5
    assert abs(ic.opbend_dist([c[0], c[5], c[4], c[3]])[0]) < 1e-5
    for i in range(1000):
        angle = np.random.uniform(-np.pi/2, np.pi/2)
        radius = np.random.uniform(0, 5*angstrom)
        #offset = np.random.uniform(0, 5*angstrom)
        c[3] = [
            radius*np.cos(angle),
            0.0,
            radius*np.sin(angle),
        ]
        assert abs(ic.opbend_cos([c[0], c[5], c[4], c[3]])[0] - np.cos(angle)) < 1e-5
        assert abs(ic.opbend_angle([c[0], c[5], c[4], c[3]])[0] - angle) < 1e-5
