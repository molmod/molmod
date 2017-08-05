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
"""Cartesian vector manipulation routines"""


from __future__ import division

import numpy as np


__all__ = [
    "cosine", "angle", "random_unit", "random_orthonormal", "triangle_normal",
]



def cosine(a, b):
    """Compute the cosine between two vectors

       The result is clipped within the range [-1, 1]
    """
    result = np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)
    return np.clip(result, -1, 1)


def angle(a, b):
    """Compute the angle between two vectors

       The result is clipped within the range [-1, 1]
    """
    return np.arccos(cosine(a, b))


def random_unit(size=3):
    """Return a random unit vector of the given dimension

       Optional argument:
         size  --  the number of dimensions of the unit vector [default=3]
    """
    while True:
        result = np.random.normal(0, 1, size)
        length = np.linalg.norm(result)
        if length > 1e-3:
            return result/length


normal_fns = [
    lambda a: np.array([0.0, -a[2], a[1]]),
    lambda a: np.array([a[2], 0.0, -a[0]]),
    lambda a: np.array([-a[1], a[0], 0.0])
]

def random_orthonormal(normal):
    """Return a random normalized vector orthogonal to the given vector"""
    u = normal_fns[np.argmin(np.fabs(normal))](normal)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    v /= np.linalg.norm(v)
    alpha = np.random.uniform(0.0, np.pi*2)
    return np.cos(alpha)*u + np.sin(alpha)*v

def triangle_normal(a, b, c):
    """Return a vector orthogonal to the given triangle

       Arguments:
         a, b, c  --  three 3D numpy vectors
    """
    normal = np.cross(a - c, b - c)
    norm = np.linalg.norm(normal)
    return normal/norm
