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
"""Cartesian vector manipulation routines"""

import numpy


__all__= [
    "cosine", "angle", "random_unit", "random_normal", "random_orthonormal",
    "triangle_normal",
]



def cosine(a, b):
    """Compute the cosine between two vectors

       The result is clipped within the range [-1, 1]
    """
    result = numpy.dot(a, b) / numpy.linalg.norm(a) / numpy.linalg.norm(b)
    return numpy.clip(result, -1, 1)


def angle(a, b):
    """Compute the angle between two vectors

       The result is clipped within the range [-1, 1]
    """
    return numpy.arccos(cosine(a, b))


def random_unit(size):
    """Return a random unit vector of the given dimension"""
    while True:
        result = numpy.random.normal(0, 1, size)
        length = numpy.linalg.norm(result)
        if length > 1e-3:
            return result/length

def random_normal():
    """Return a random 3D unit vector"""
    return random_unit(3)


normal_fns = [
    lambda a: numpy.array([0.0, -a[2], a[1]]),
    lambda a: numpy.array([a[2], 0.0, -a[0]]),
    lambda a: numpy.array([-a[1], a[0], 0.0])
]

def random_orthonormal(normal):
    """Return a random normalized vector orthogonal to the given vector"""
    u = normal_fns[numpy.argmin(numpy.fabs(normal))](normal)
    u /= numpy.linalg.norm(u)
    v = numpy.cross(normal, u)
    v /= numpy.linalg.norm(v)
    alpha = numpy.random.uniform(0.0, numpy.pi*2)
    return numpy.cos(alpha)*u + numpy.sin(alpha)*v

def triangle_normal(a, b, c):
    """Return a vector orthogonal to the given triangle

       Arguments:
         a, b, c  --  three 3D numpy vectors
    """
    normal = numpy.cross(a - c, b - c)
    norm = numpy.linalg.norm(normal)
    return normal/norm



