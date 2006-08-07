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

import math


__all__ = ["quaternion_product", "conjugate", "conjugated", "apply_to_vector",
           "quaternion_from_rotation_matrix", "quaternion_to_rotation_matrix"]


def quaternion_product(quat1, quat2):
    return numpy.array(
        [quat1[0]*quat2[0] - numpy.dot(quat1[1:], quat2[1:]),
         quat1[0]*quat2[1] + quat2[0]*quat1[1] + quat1[2]*quat2[3] - quat1[3]*quat2[2],
         quat1[0]*quat2[2] + quat2[0]*quat1[2] + quat1[3]*quat2[1] - quat1[1]*quat2[3],
         quat1[0]*quat2[3] + quat2[0]*quat1[3] + quat1[1]*quat2[2] - quat1[2]*quat2[1]],
        float)


def conjugate(quat):
    quat[1:] *= -1


def conjugated(quat):
    result = quat.copy()
    conjugate(result)
    return result


def apply_to_vector(quat, vector):
    dp = numpy.dot(quat[1:], vector)
    cos = (2*quat[0]*quat[0] - 1)
    return numpy.array(
        [2 * (quat[0] * (quat[2] * vector[2] - quat[3] * vector[1]) + quat[1] * dp) + cos * vector[0],
         2 * (quat[0] * (quat[3] * vector[0] - quat[1] * vector[2]) + quat[2] * dp) + cos * vector[1],
         2 * (quat[0] * (quat[1] * vector[1] - quat[2] * vector[0]) + quat[3] * dp) + cos * vector[2]],
        float)


off_diagonals = [[2,1], [0,2], [1,0]]

def quaternion_from_rotation_matrix(rotation_matrix):
    invert = (numpy.linalg.det(rotation_matrix) < 0)
    if invert:
        factor = -1
    else:
        factor = 1
    c2 = 0.25*(factor*numpy.trace(rotation_matrix) + 1)
    if c2 < 0:
        print c2
        c2 = 0.0
    c = math.sqrt(c2)
    r2 = 0.5*(1 + factor*numpy.diagonal(rotation_matrix)) - c2
    r = numpy.zeros(3, float)
    for index, r_comp in enumerate(r2):
        if r_comp < 0:
            continue
        else:
            row, col = off_diagonals[index]
            if (rotation_matrix[row, col] + rotation_matrix[col, row] < 0):
                r[index] = -r_comp
            else:
                r[index] = +r_comp
    return factor, numpy.array([c, r[0], r[1], r[2]], float)


def quaternion_to_rotation_matrix(quaternion):
    c, x, y, z = quaternion
    return 2*numpy.array(
        [[0.5 - y*y - z*z, x*y - c*z,       x*z + c*y      ],
         [x*y + c*z      , 0.5 - x*x - z*z, y*z - c*x      ],
         [x*z - c*y      , y*z + c*x      , 0.5 - x*x - y*y]],
        float)

