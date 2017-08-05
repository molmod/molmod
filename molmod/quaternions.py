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

"""Auxiliary routines to work with quaternions

Quaternions are represented by numpy vectors with four real floating point
values. The routines below focus on the link between quaternions and three-
dimensional rotations.
"""


import numpy as np

__all__ = [
  "quaternion_product", "conjugate", "quaternion_rotation",
  "rotation_matrix_to_quaternion", "quaternion_to_rotation_matrix"
]


def quaternion_product(quat1, quat2):
    """Return the quaternion product of the two arguments"""
    return np.array([
        quat1[0]*quat2[0] - np.dot(quat1[1:], quat2[1:]),
        quat1[0]*quat2[1] + quat2[0]*quat1[1] + quat1[2]*quat2[3] - quat1[3]*quat2[2],
        quat1[0]*quat2[2] + quat2[0]*quat1[2] + quat1[3]*quat2[1] - quat1[1]*quat2[3],
        quat1[0]*quat2[3] + quat2[0]*quat1[3] + quat1[1]*quat2[2] - quat1[2]*quat2[1]
    ], float)


def conjugate(quat):
    """Return the conjugate quaternion"""
    result = quat.copy()
    result[1:] *= -1
    conjugate(result)
    return result


def quaternion_rotation(quat, vector):
    """Apply the rotation represented by the quaternion to the vector

       Warning: This only works correctly for normalized quaternions.
    """
    dp = np.dot(quat[1:], vector)
    cos = (2*quat[0]*quat[0] - 1)
    return np.array([
        2 * (quat[0] * (quat[2] * vector[2] - quat[3] * vector[1]) + quat[1] * dp) + cos * vector[0],
        2 * (quat[0] * (quat[3] * vector[0] - quat[1] * vector[2]) + quat[2] * dp) + cos * vector[1],
        2 * (quat[0] * (quat[1] * vector[1] - quat[2] * vector[0]) + quat[3] * dp) + cos * vector[2]
    ], float)


off_diagonals = [[2, 1], [0, 2], [1, 0]]

def rotation_matrix_to_quaternion(rotation_matrix):
    """Compute the quaternion representing the rotation given by the matrix"""
    invert = (np.linalg.det(rotation_matrix) < 0)
    if invert:
        factor = -1
    else:
        factor = 1
    c2 = 0.25*(factor*np.trace(rotation_matrix) + 1)
    if c2 < 0:
        #print c2
        c2 = 0.0
    c = np.sqrt(c2)
    r2 = 0.5*(1 + factor*np.diagonal(rotation_matrix)) - c2
    #print "check", r2.sum()+c2
    r = np.zeros(3, float)
    for index, r2_comp in enumerate(r2):
        if r2_comp < 0:
            continue
        else:
            row, col = off_diagonals[index]
            if (rotation_matrix[row, col] - rotation_matrix[col, row] < 0):
                r[index] = -np.sqrt(r2_comp)
            else:
                r[index] = +np.sqrt(r2_comp)
    return factor, np.array([c, r[0], r[1], r[2]], float)


def quaternion_to_rotation_matrix(quaternion):
    """Compute the rotation matrix representated by the quaternion"""
    c, x, y, z = quaternion
    return np.array([
        [c*c + x*x - y*y - z*z, 2*x*y - 2*c*z,         2*x*z + 2*c*y        ],
        [2*x*y + 2*c*z,         c*c - x*x + y*y - z*z, 2*y*z - 2*c*x        ],
        [2*x*z - 2*c*y,         2*y*z + 2*c*x,         c*c - x*x - y*y + z*z]
    ], float)
