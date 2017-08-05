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
"""Data structures to handle 3D rotations and translations

In addition to Translation, Rotation and Complete classes, two utility
functions are provided: rotation_around_center and superpose. The latter is an
implementation of the Kabsch algorithm.
"""


from __future__ import division

from builtins import range
import numpy as np

from molmod.utils import cached, ReadOnly, ReadOnlyAttribute, compute_rmsd
from molmod.vectors import random_unit
from molmod.unit_cells import UnitCell


__all__ = [
    "Translation", "Rotation", "Complete", "superpose", "fit_rmsd"
]


eps = 1.0e-6

def check_matrix(m):
    """Check the sanity of the given 4x4 transformation matrix"""
    if m.shape != (4, 4):
        raise ValueError("The argument must be a 4x4 array.")
    if max(abs(m[3, 0:3])) > eps:
        raise ValueError("The given matrix does not have correct translational part")
    if abs(m[3, 3] - 1.0) > eps:
        raise ValueError("The lower right element of the given matrix must be 1.0.")


class Translation(ReadOnly):
    """Represents a translation in 3D

       The attribute t contains the actual translation vector, which is a numpy
       array with three elements.
    """
    t = ReadOnlyAttribute(np.ndarray, none=False, npdim=1, npshape=(3,),
        npdtype=float, doc="the translation vector")

    def __init__(self, t):
        """
           Argument:
            | ``t``  --  translation vector, a list-like object with three
                         numbers
        """
        self.t = t

    @classmethod
    def from_matrix(cls, m):
        """Initialize a translation from a 4x4 matrix"""
        check_matrix(m)
        return cls(m[0:3, 3])

    @classmethod
    def identity(cls):
        """Return the identity transformation"""
        return cls(np.zeros(3, float))

    @cached
    def matrix(self):
        """The 4x4 matrix representation of this translation"""
        result = np.identity(4, float)
        result[0:3, 3] = self.t
        return result

    @cached
    def inv(self):
        """The inverse translation"""
        result = Translation(-self.t)
        result._cache_inv = self
        return result

    def apply_to(self, x, columns=False):
        """Apply this translation to the given object

           The argument can be several sorts of objects:

           * ``np.array`` with shape (3, )
           * ``np.array`` with shape (N, 3)
           * ``np.array`` with shape (3, N), use ``columns=True``
           * ``Translation``
           * ``Rotation``
           * ``Complete``
           * ``UnitCell``

           In case of arrays, the 3D vectors are translated. In case of trans-
           formations, a new transformation is returned that consists of this
           translation applied AFTER the given translation. In case of a unit
           cell, the original object is returned.

           This method is equivalent to ``self*x``.
        """
        if isinstance(x, np.ndarray) and len(x.shape) == 2 and x.shape[0] == 3 and columns:
            return x + self.t.reshape((3,1))
        if isinstance(x, np.ndarray) and (x.shape == (3, ) or (len(x.shape) == 2 and x.shape[1] == 3)) and not columns:
            return x + self.t
        elif isinstance(x, Complete):
            return Complete(x.r, x.t + self.t)
        elif isinstance(x, Translation):
            return Translation(x.t + self.t)
        elif isinstance(x, Rotation):
            return Complete(x.r, self.t)
        elif isinstance(x, UnitCell):
            return x
        else:
            raise ValueError("Can not apply this translation to %s" % x)

    __mul__ = apply_to

    def compare(self, other, t_threshold=1e-3):
        """Compare two translations

           The RMSD of the translation vectors is computed. The return value
           is True when the RMSD is below the threshold, i.e. when the two
           translations are almost identical.
        """
        return compute_rmsd(self.t, other.t) < t_threshold


class Rotation(ReadOnly):
    """Represents a rotation in 3D about the origin

       The attribute r contains the actual rotation matrix, which is a numpy
       array with shape (3, 3).
    """
    def _check_r(self, r):
        """the columns must orthogonal"""
        if abs(np.dot(r[:, 0], r[:, 0]) - 1) > eps or \
            abs(np.dot(r[:, 0], r[:, 0]) - 1) > eps or \
            abs(np.dot(r[:, 0], r[:, 0]) - 1) > eps or \
            np.dot(r[:, 0], r[:, 1]) > eps or \
            np.dot(r[:, 1], r[:, 2]) > eps or \
            np.dot(r[:, 2], r[:, 0]) > eps:
            raise ValueError("The rotation matrix is significantly non-orthonormal.")


    r = ReadOnlyAttribute(np.ndarray, none=False, check=_check_r, npdim=2,
        npshape=(3,3), npdtype=float, doc="the rotation matrix")

    def __init__(self, r):
        """
           Argument:
            | ``r``  --  rotation matrix, a 3 by 3 orthonormal array-like object
        """
        self.r = r

    @classmethod
    def from_matrix(cls, m):
        """Initialize a rotation from a 4x4 matrix"""
        check_matrix(m)
        return cls(m[0:3, 0:3])

    @classmethod
    def identity(cls):
        """Return the identity transformation"""
        return cls(np.identity(3, float))

    @classmethod
    def random(cls):
        """Return a random rotation"""
        axis = random_unit()
        angle = np.random.uniform(0,2*np.pi)
        invert = bool(np.random.randint(0,2))
        return Rotation.from_properties(angle, axis, invert)

    @classmethod
    def from_properties(cls, angle, axis, invert):
        """Initialize a rotation based on the properties"""
        norm = np.linalg.norm(axis)
        if norm > 0:
            x = axis[0] / norm
            y = axis[1] / norm
            z = axis[2] / norm
            c = np.cos(angle)
            s = np.sin(angle)
            r = (1-2*invert) * np.array([
                [x*x*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s],
                [x*y*(1-c)+z*s, y*y*(1-c)+c  , y*z*(1-c)-x*s],
                [x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c  ]
            ])
        else:
            r = np.identity(3) * (1-2*invert)
        return cls(r)

    @cached
    def properties(self):
        """Rotation properties: angle, axis, invert"""
        # determine wether an inversion rotation has been applied
        invert = (np.linalg.det(self.r) < 0)
        factor = {True: -1, False: 1}[invert]
        # get the rotation data
        # trace(r) = 1+2*cos(angle)
        cos_angle = 0.5*(factor*np.trace(self.r) - 1)
        if cos_angle > 1: cos_angle = 1.0
        if cos_angle < -1: cos_angle = -1.0
        # the antisymmetric part of the non-diagonal vector tell us something
        # about sin(angle) and n.
        axis = 0.5*factor*np.array([-self.r[1, 2] + self.r[2, 1], self.r[0, 2] - self.r[2, 0], -self.r[0, 1] + self.r[1, 0]])
        sin_angle = np.linalg.norm(axis)
        # look for the best way to normalize the
        if (sin_angle == 0.0) and (cos_angle > 0):
            axis[2] = 1.0
        elif abs(sin_angle) < (1-cos_angle):
            for index in range(3):
                axis[index] = {True: -1, False: 1}[axis[index] < 0] * np.sqrt(abs((factor*self.r[index, index] - cos_angle) / (1 - cos_angle)))
        else:
            axis = axis / sin_angle

        # Finally calculate the angle:
        angle = np.arctan2(sin_angle, cos_angle)
        return angle, axis, invert

    @cached
    def matrix(self):
        """The 4x4 matrix representation of this rotation"""
        result = np.identity(4, float)
        result[0:3, 0:3] = self.r
        return result

    @cached
    def inv(self):
        """The inverse rotation"""
        result = Rotation(self.r.transpose())
        result._cache_inv = self
        return result

    def apply_to(self, x, columns=False):
        """Apply this rotation to the given object

           The argument can be several sorts of objects:

           * ``np.array`` with shape (3, )
           * ``np.array`` with shape (N, 3)
           * ``np.array`` with shape (3, N), use ``columns=True``
           * ``Translation``
           * ``Rotation``
           * ``Complete``
           * ``UnitCell``

           In case of arrays, the 3D vectors are rotated. In case of trans-
           formations, a transformation is returned that consists of this
           rotation applied AFTER the given translation. In case of a unit cell,
           a unit cell with rotated cell vectors is returned.

           This method is equivalent to ``self*x``.
        """
        if isinstance(x, np.ndarray) and len(x.shape) == 2 and x.shape[0] == 3 and columns:
            return np.dot(self.r, x)
        if isinstance(x, np.ndarray) and (x.shape == (3, ) or (len(x.shape) == 2 and x.shape[1] == 3)) and not columns:
            return np.dot(x, self.r.transpose())
        elif isinstance(x, Complete):
            return Complete(np.dot(self.r, x.r), np.dot(self.r, x.t))
        elif isinstance(x, Translation):
            return Complete(self.r, np.dot(self.r, x.t))
        elif isinstance(x, Rotation):
            return Rotation(np.dot(self.r, x.r))
        elif isinstance(x, UnitCell):
            return UnitCell(np.dot(self.r, x.matrix), x.active)
        else:
            raise ValueError("Can not apply this rotation to %s" % x)

    __mul__ = apply_to

    def compare(self, other, r_threshold=1e-3):
        """Compare two rotations

           The RMSD of the rotation matrices is computed. The return value
           is True when the RMSD is below the threshold, i.e. when the two
           rotations are almost identical.
        """
        return compute_rmsd(self.r, other.r) < r_threshold


class Complete(Translation, Rotation):
    """Represents a rotation and translation in 3D

       The attribute t contains the actual translation vector, which is a numpy
       array with three elements. The attribute r contains the actual rotation
       matrix, which is a numpy array with shape (3, 3).

       Internally the translation part is always applied after the rotation
       part.
    """
    def __init__(self, r, t):
        """
           Arguments:
            | ``r``  --  rotation matrix, a 3 by 3 orthonormal array-like object
            | ``t``  --  translation vector, a list-like object with three
                         numbers
        """
        Translation.__init__(self, t)
        Rotation.__init__(self, r)

    @classmethod
    def from_matrix(cls, m):
        """Initialize a complete transformation from a 4x4 matrix"""
        check_matrix(m)
        return cls(m[0:3, 0:3], m[0:3, 3])

    @classmethod
    def identity(cls):
        """Return the identity transformation"""
        return cls(np.identity(3, float), np.zeros(3, float))

    @classmethod
    def from_properties(cls, angle, axis, invert, translation):
        """Initialize a transformation based on the properties"""
        rot = Rotation.from_properties(angle, axis, invert)
        return Complete(rot.r, translation)

    @classmethod
    def cast(cls, c):
        """Convert the first argument into a Complete object"""
        if isinstance(c, Complete):
            return c
        elif isinstance(c, Translation):
            return Complete(np.identity(3, float), c.t)
        elif isinstance(c, Rotation):
            return Complete(c.r, np.zeros(3, float))

    @classmethod
    def about_axis(cls, center, angle, axis, invert=False):
        """Create transformation that represents a rotation about an axis

           Arguments:
            | ``center``  --  Point on the axis
            | ``angle``  --  Rotation angle
            | ``axis``  --  Rotation axis
            | ``invert``  --  When True, an inversion rotation is constructed
                              [default=False]
        """
        return Translation(center) * \
               Rotation.from_properties(angle, axis, invert) * \
               Translation(-center)

    @cached
    def matrix(self):
        """The 4x4 matrix representation of this transformation"""
        result = np.identity(4, float)
        result[0:3, 3] = self.t
        result[0:3, 0:3] = self.r
        return result

    @cached
    def properties(self):
        """Transformation properties: angle, axis, invert, translation"""
        rot = Rotation(self.r)
        angle, axis, invert = rot.properties
        return angle, axis, invert, self.t

    @cached
    def inv(self):
        """The inverse transformation"""
        result = Complete(self.r.transpose(), np.dot(self.r.transpose(), -self.t))
        result._cache_inv = self
        return result

    def apply_to(self, x, columns=False):
        """Apply this transformation to the given object

           The argument can be several sorts of objects:

           * ``np.array`` with shape (3, )
           * ``np.array`` with shape (N, 3)
           * ``np.array`` with shape (3, N), use ``columns=True``
           * ``Translation``
           * ``Rotation``
           * ``Complete``
           * ``UnitCell``

           In case of arrays, the 3D vectors are transformed. In case of trans-
           formations, a transformation is returned that consists of this
           transformation applied AFTER the given translation. In case of a unit
           cell, a unit cell with rotated cell vectors is returned. (The
           translational part does not affect the unit cell.)

           This method is equivalent to self*x.
        """
        if isinstance(x, np.ndarray) and len(x.shape) == 2 and x.shape[0] == 3 and columns:
            return np.dot(self.r, x) + self.t.reshape((3,1))
        if isinstance(x, np.ndarray) and (x.shape == (3, ) or (len(x.shape) == 2 and x.shape[1] == 3)) and not columns:
            return np.dot(x, self.r.transpose()) + self.t
        elif isinstance(x, Complete):
            return Complete(np.dot(self.r, x.r), np.dot(self.r, x.t) + self.t)
        elif isinstance(x, Translation):
            return Complete(self.r, np.dot(self.r, x.t) + self.t)
        elif isinstance(x, Rotation):
            return Complete(np.dot(self.r, x.r), self.t)
        elif isinstance(x, UnitCell):
            return UnitCell(np.dot(self.r, x.matrix), x.active)
        else:
            raise ValueError("Can not apply this rotation to %s" % x)

    __mul__ = apply_to

    def compare(self, other, t_threshold=1e-3, r_threshold=1e-3):
        """Compare two transformations

           The RMSD values of the rotation matrices and the translation vectors
           are computed. The return value is True when the RMSD values are below
           the thresholds, i.e. when the two transformations are almost
           identical.
        """
        return compute_rmsd(self.t, other.t) < t_threshold and compute_rmsd(self.r, other.r) < r_threshold


def superpose(ras, rbs, weights=None):
    """Compute the transformation that minimizes the RMSD between the points ras and rbs

       Arguments:
        | ``ras``  --  a ``np.array`` with 3D coordinates of geometry A,
                       shape=(N,3)
        | ``rbs``  --  a ``np.array`` with 3D coordinates of geometry B,
                       shape=(N,3)

       Optional arguments:
        | ``weights``  --  a numpy array with fitting weights for each
                           coordinate, shape=(N,)

       Return value:
        | ``transformation``  --  the transformation that brings geometry A into
                                  overlap with geometry B

       Each row in ras and rbs represents a 3D coordinate. Corresponding rows
       contain the points that are brought into overlap by the fitting
       procedure. The implementation is based on the Kabsch Algorithm:

       http://dx.doi.org/10.1107%2FS0567739476001873
    """
    if weights is None:
        ma = ras.mean(axis=0)
        mb = rbs.mean(axis=0)
    else:
        total_weight = weights.sum()
        ma = np.dot(weights, ras)/total_weight
        mb = np.dot(weights, rbs)/total_weight


    # Kabsch
    if weights is None:
        A = np.dot((rbs-mb).transpose(), ras-ma)
    else:
        weights = weights.reshape((-1, 1))
        A = np.dot(((rbs-mb)*weights).transpose(), (ras-ma)*weights)
    v, s, wt = np.linalg.svd(A)
    s[:] = 1
    if np.linalg.det(np.dot(v, wt)) < 0:
        s[2] = -1
    r = np.dot(wt.T*s, v.T)
    return Complete(r, np.dot(r, -mb) + ma)


def fit_rmsd(ras, rbs, weights=None):
    """Fit geometry rbs onto ras, returns more info than superpose

       Arguments:
        | ``ras``  --  a numpy array with 3D coordinates of geometry A,
                       shape=(N,3)
        | ``rbs``  --  a numpy array with 3D coordinates of geometry B,
                       shape=(N,3)

       Optional arguments:
        | ``weights``  --  a numpy array with fitting weights for each
                           coordinate, shape=(N,)

       Return values:
        | ``transformation``  --  the transformation that brings geometry A into
                                  overlap with geometry B
        | ``rbs_trans``  --  the transformed coordinates of geometry B
        | ``rmsd``  --  the rmsd of the distances between corresponding atoms in
                        geometry A and B

       This is a utility routine based on the function superpose. It just
       computes rbs_trans and rmsd after calling superpose with the same
       arguments
    """
    transformation = superpose(ras, rbs, weights)
    rbs_trans = transformation * rbs
    rmsd = compute_rmsd(ras, rbs_trans)
    return transformation, rbs_trans, rmsd
