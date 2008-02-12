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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


import numpy

import math, random


__all__= [
    "cosine", "angle", "random_orthonormal", "trivial_orthonormal",
    "triangle_normal", "random_normal",
]



def cosine(a, b):
    result = numpy.dot(a, b) / math.sqrt(numpy.dot(a, a) * numpy.dot(b, b))
    if result <= -1: return -1
    elif result >= 1: return 1
    else: return result


def angle(a, b):
    return math.acos(cosine(a, b))


def random_unit(size):
    while True:
        result = numpy.random.uniform(-1, 1, size)
        length = numpy.linalg.norm(result)
        if length <= 1.0 and length > 1e-3:
            return result/length


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


def random_normal():
    while True:
        result = numpy.random.uniform(-1,1,3)
        norm = numpy.linalg.norm(result)
        if norm < 1 and norm > 1e-4:
            return result / norm


