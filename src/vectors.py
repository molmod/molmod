# Zeobuilder is an extensible GUI-toolkit for molecular model construction.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of Zeobuilder.
#
# Zeobuilder is free software; you can redistribute it and/or
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


import numpy

import math, random


__all__= ["cosine", "angle", "random_orthonormal", "trivial_orthonormal", "triangle_normal"]


def cosine(a, b):
    result = numpy.dot(a, b) / math.sqrt(numpy.dot(a, a) * numpy.dot(b, b))
    if result <= -1: return -1
    elif result >= 1: return 1
    else: return result


def angle(a, b):
    return math.acos(cosine(a, b))


normal_fns = [
    lambda a: numpy.array([0.0, -a[2], a[1]]),
    lambda a: numpy.array([a[2], 0.0, -a[0]]),
    lambda a: numpy.array([-a[1], a[0], 0.0])
]

def random_orthonormal(normal):
    u = normal_fns[numpy.argmin(numpy.fabs(normal))](normal)
    u /= math.sqrt(numpy.dot(u, u))
    v = numpy.cross(normal, u)
    angle = random.uniform(0.0, math.pi*2)
    return math.cos(angle)*u + math.sin(angle)*v

def trivial_orthonormal(normal):
    u = normal_fns[numpy.argmin(numpy.fabs(normal))](normal)
    u /= math.sqrt(numpy.dot(u, u))
    return u


def triangle_normal(a, b, c):
    normal = numpy.cross(a - c, b - c)
    norm = math.sqrt(numpy.dot(normal, normal))
    if norm <= 1e-8:
        return numpy.zeros(3, float)
    else:
        normal /= norm
        return normal

