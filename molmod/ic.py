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
"""Evaluation of internal coordinates with first and second order derivatives

This implementation is pure python and sacrifices computational efficiency on
the altar of programming flexibility. It is really easy to implement new types
of internal coordinates since one only has to enter the formula that evaluates
the internal coordinate. First and second order derivatives towards Cartesian
coordinates require only a minimum of extra work.

Two auxiliary classes Scalar and Vector3 support most of the mathematical
operations required to compute the internal coordinates. Additionally they also
know the chain rule for each operation and can therefore evaluate the
derivatives simultaneously.
"""


from __future__ import division

import numpy as np


__all__ = [
    "Scalar", "Vector3", "dot", "cross",
    "bond_length", "pair_distance",
    "bend_cos", "bend_angle",
    "dihed_cos", "dihed_angle",
    "opbend_cos", "opbend_angle", "opbend_dist",
    "opbend_mcos", "opbend_mangle",
]


#
# Auxiliary classes
#

class Scalar(object):
    """A scalar object with optional first and second order derivates

       Each input value to which the derivative is computed has its own index.
       The numerical value of the derivatives are stored in arrays self.d and
       self.dd. The value of the scalar itself if self.v
    """

    def __init__(self, size, deriv=0, value=0, index=None):
        """
           Arguments:
            | ``size`` -- The number of inputs on which this ic depends. e.g. a
                          distance depends on 6 Cartesian coordinates.
            | ``deriv`` -- Consider up to deriv order derivatives. (max=2)
            | ``value`` -- The initial value.
            | ``index`` -- If this scalar is one of the input variables, this is
                           its index.

           The scalar object supports several in place modifications.
        """
        self.deriv = deriv
        self.size = size
        self.v = value
        if deriv > 0:
            self.d = np.zeros(size, float)
            if index is not None:
                self.d[index] = 1
        if deriv > 1:
            self.dd = np.zeros((size, size), float)
        if deriv > 2:
            raise ValueError("This implementation (only) supports up to second order derivatives.")

    def copy(self):
        """Return a deep copy"""
        result = Scalar(self.size, self.deriv)
        result.v = self.v
        if self.deriv > 0: result.d[:] = self.d[:]
        if self.deriv > 1: result.dd[:] = self.dd[:]
        return result

    def results(self):
        """Return the value and optionally derivative and second order derivative"""
        if self.deriv == 0:
            return self.v,
        if self.deriv == 1:
            return self.v, self.d
        if self.deriv == 2:
            return self.v, self.d, self.dd

    def __iadd__(self, other):
        if self.deriv > 1: self.dd += other.dd
        if self.deriv > 0: self.d += other.d
        self.v += other.v
        return self

    def __add__(self, other):
        result = self.copy()
        result += other
        return result

    def __isub__(self, other):
        if self.deriv > 1: self.dd -= other.dd
        if self.deriv > 0: self.d -= other.d
        self.v -= other.v
        return self

    def __sub__(self, other):
        result = self.copy()
        result -= other
        return result

    def __imul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            self.v *= other
            if self.deriv > 0:
                self.d *= other
            if self.deriv > 1:
                self.dd *= other
        elif isinstance(other, Scalar):
            # trying to avoid temporaries as much as possible
            if self.deriv > 1:
                self.dd *= other.v
                self.dd += self.v*other.dd
                tmp = np.outer(self.d, other.d)
                self.dd += tmp
                self.dd += tmp.transpose()
            if self.deriv > 0:
                self.d *= other.v
                self.d += self.v*other.d
            self.v *= other.v
        else:
            raise TypeError("Second argument must be float, int or Scalar")
        return self

    def __mul__(self, other):
        result = self.copy()
        result *= other
        return result

    def __itruediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            self.v /= other
            if self.deriv > 0:
                self.d /= other
            if self.deriv > 1:
                self.dd /= other
        elif isinstance(other, Scalar):
            # trying to avoid temporaries as much as possible
            self.v /= other.v
            if self.deriv > 0:
                self.d -= self.v*other.d
                self.d /= other.v
            if self.deriv > 1:
                self.dd -= self.v*other.dd
                tmp = np.outer(self.d, other.d)
                self.dd -= tmp
                self.dd -= tmp.transpose()
                self.dd /= other.v
        else:
            raise TypeError("Second argument must be float, int or Scalar")
        return self

    def inv(self):
        """In place invert"""
        self.v = 1/self.v
        tmp = self.v**2
        if self.deriv > 1:
            self.dd[:] = tmp*(2*self.v*np.outer(self.d, self.d) - self.dd)
        if self.deriv > 0:
            self.d[:] = -tmp*self.d[:]


class Vector3(object):
    """A Three dimensional vector with optional first and second order derivatives.

       This object is nothing more than a tier for three Scalar objects.
    """

    def __init__(self, size, deriv=0, values=(0, 0, 0), indexes=(None, None, None)):
        """
           Arguments:
            | ``size`` -- The number of inputs on which this ic depends. e.g. a
                          distance depends on 6 Cartesian coordinates.
            | ``deriv`` -- Consider up to deriv order derivatives. (max=2)
            | ``values`` -- The initial values.
            | ``indexes`` -- If this vector is one of the input variables, these
                             are the indexes of the components.
        """
        self.deriv = deriv
        self.size = size
        self.x = Scalar(size, deriv, values[0], indexes[0])
        self.y = Scalar(size, deriv, values[1], indexes[1])
        self.z = Scalar(size, deriv, values[2], indexes[2])

    def copy(self):
        """Return a deep copy"""
        result = Vector3(self.size, self.deriv)
        result.x.v = self.x.v
        result.y.v = self.y.v
        result.z.v = self.z.v
        if self.deriv > 0:
            result.x.d[:] = self.x.d
            result.y.d[:] = self.y.d
            result.z.d[:] = self.z.d
        if self.deriv > 1:
            result.x.dd[:] = self.x.dd
            result.y.dd[:] = self.y.dd
            result.z.dd[:] = self.z.dd
        return result

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def __itruediv__(self, other):
        self.x /= other
        self.y /= other
        self.z /= other
        return self

    def norm(self):
        """Return a Scalar object with the norm of this vector"""
        result = Scalar(self.size, self.deriv)
        result.v = np.sqrt(self.x.v**2 + self.y.v**2 + self.z.v**2)
        if self.deriv > 0:
            result.d += self.x.v*self.x.d
            result.d += self.y.v*self.y.d
            result.d += self.z.v*self.z.d
            result.d /= result.v
        if self.deriv > 1:
            result.dd += self.x.v*self.x.dd
            result.dd += self.y.v*self.y.dd
            result.dd += self.z.v*self.z.dd
            denom = result.v**2
            result.dd += (1 - self.x.v**2/denom)*np.outer(self.x.d, self.x.d)
            result.dd += (1 - self.y.v**2/denom)*np.outer(self.y.d, self.y.d)
            result.dd += (1 - self.z.v**2/denom)*np.outer(self.z.d, self.z.d)
            tmp = -self.x.v*self.y.v/denom*np.outer(self.x.d, self.y.d)
            result.dd += tmp+tmp.transpose()
            tmp = -self.y.v*self.z.v/denom*np.outer(self.y.d, self.z.d)
            result.dd += tmp+tmp.transpose()
            tmp = -self.z.v*self.x.v/denom*np.outer(self.z.d, self.x.d)
            result.dd += tmp+tmp.transpose()
            result.dd /= result.v
        return result


#
# Auxiliary functions
#


def dot(r1, r2):
    """Compute the dot product

       Arguments:
        | ``r1``, ``r2``  -- two :class:`Vector3` objects

       (Returns a Scalar)
    """
    if r1.size != r2.size:
        raise ValueError("Both arguments must have the same input size.")
    if r1.deriv != r2.deriv:
        raise ValueError("Both arguments must have the same deriv.")
    return r1.x*r2.x + r1.y*r2.y + r1.z*r2.z


def cross(r1, r2):
    """Compute the cross product

       Arguments:
        | ``r1``, ``r2``  -- two :class:`Vector3` objects

       (Returns a Vector3)
    """
    if r1.size != r2.size:
        raise ValueError("Both arguments must have the same input size.")
    if r1.deriv != r2.deriv:
        raise ValueError("Both arguments must have the same deriv.")
    result = Vector3(r1.size, r1.deriv)
    result.x = r1.y*r2.z - r1.z*r2.y
    result.y = r1.z*r2.x - r1.x*r2.z
    result.z = r1.x*r2.y - r1.y*r2.x
    return result


#
# Internal coordinate functions
#

def bond_length(rs, deriv=0):
    """Compute the distance between the two points rs[0] and rs[1]

       Arguments:
        | ``rs``  --  two numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]

       When derivatives are computed a tuple with a single result is returned
    """
    return _bond_transform(rs, _bond_length_low, deriv)

pair_distance = bond_length


def bend_cos(rs, deriv=0):
    """Compute the cosine of the angle between the vectors rs[0]-rs[1] and rs[2]-rs[1]

       Arguments:
        | ``rs``  --  three numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]

       When derivatives are computed a tuple with a single result is returned
    """
    return _bend_transform(rs, _bend_cos_low, deriv)


def bend_angle(rs, deriv=0):
    """Compute the angle between the vectors rs[0]-rs[1] and rs[2]-rs[1]

       Arguments:
        | ``rs``  --  three numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]

       When derivatives are computed a tuple with a single result is returned
    """
    return _bend_transform(rs, _bend_angle_low, deriv)


def dihed_cos(rs, deriv=0):
    """Compute the cosine of the angle between the planes rs[0], rs[1], rs[2] and rs[1], rs[2], rs[3]

       Arguments:
        | ``rs``  --  four numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]
    """
    return _dihed_transform(rs, _dihed_cos_low, deriv)


def dihed_angle(rs, deriv=0):
    """Compute the angle between the planes rs[0], rs[1], rs[2] and rs[1], rs[2], rs[3]

       The sign convention corresponds to the IUPAC definition of the torsion
       angle: http://dx.doi.org/10.1351/goldbook.T06406

       Arguments:
        | ``rs``  --  four numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]

       When derivatives are computed a tuple with a single result is returned
    """
    return _dihed_transform(rs, _dihed_angle_low, deriv)


def opbend_dist(rs, deriv=0):
    """Compute the out-of-plane distance, i.e. the distance between atom rs[3] and plane rs[0],rs[1],rs[2]

       Arguments:
        | ``rs``  --  four numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]
    """
    return _opbend_transform(rs, _opdist_low, deriv)


def opbend_cos(rs, deriv=0):
    """Compute the cosine of the angle between the vector (rs[0],rs[3]) and plane rs[0],rs[1],rs[2]

       Arguments:
        | ``rs``  --  four numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]
    """
    return _opbend_transform(rs, _opbend_cos_low, deriv)


def opbend_angle(rs, deriv=0):
    """Compute the angle between the vector rs[0], rs[3] and the plane rs[0], rs[1], rs[2]

       The sign convention is as follows: positive if rs[3] lies in the space
       above plane rs[0], rs[1], rs[2] and negative if rs[3] lies below. Above
       is defined by right hand rule from rs[0]-rs[1] to rs[0]-rs[2].

       Arguments:
        | ``rs``  --  four numpy array with three elements
        | ``deriv``  --  the derivatives to be computed: 0, 1 or 2 [default=0]

       When no derivatives are computed a tuple with a single result is returned.
    """
    return _opbend_transform(rs, _opbend_angle_low, deriv)


def opbend_mangle(rs, deriv=0):
    """Compute the mean value of the 3 opbend_angles
    """
    return _opbend_transform_mean(rs, _opbend_angle_low, deriv)


def opbend_mcos(rs, deriv=0):
    """Compute the mean cos of the 3 opbend_angles
    """
    return _opbend_transform_mean(rs, _opbend_cos_low, deriv)


#
# Transformers
#


def _bond_transform(rs, fn_low, deriv):
    r = rs[0] - rs[1]
    result = fn_low(r, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = np.zeros((2, 3), float)
    d[0] = result[1]
    d[1] = -result[1]
    if deriv == 1:
        return v, d
    dd = np.zeros((2, 3, 2, 3), float)
    dd[0, :, 0, :] = result[2]
    dd[1, :, 1, :] = result[2]
    dd[0, :, 1, :] = -result[2]
    dd[1, :, 0, :] = -result[2]
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")


def _bend_transform(rs, fn_low, deriv):
    a = rs[0] - rs[1]
    b = rs[2] - rs[1]
    result = fn_low(a, b, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = np.zeros((3, 3), float)
    d[0] = result[1][:3]
    d[1] = -result[1][:3]-result[1][3:]
    d[2] = result[1][3:]
    if deriv == 1:
        return v, d
    dd = np.zeros((3, 3, 3, 3), float)
    aa = result[2][:3, :3]
    ab = result[2][:3, 3:]
    ba = result[2][3:, :3]
    bb = result[2][3:, 3:]
    dd[0, :, 0, :] =   aa
    dd[0, :, 1, :] = - aa - ab
    dd[0, :, 2, :] =   ab
    dd[1, :, 0, :] = - aa - ba
    dd[1, :, 1, :] =   aa + ba + ab + bb
    dd[1, :, 2, :] = - ab - bb
    dd[2, :, 0, :] =   ba
    dd[2, :, 1, :] = - ba - bb
    dd[2, :, 2, :] =   bb
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")


def _dihed_transform(rs, fn_low, deriv):
    a = rs[0] - rs[1]
    b = rs[2] - rs[1]
    c = rs[3] - rs[2]
    result = fn_low(a, b, c, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = np.zeros((4, 3), float)
    d[0] = result[1][:3]
    d[1] = -result[1][:3]-result[1][3:6]
    d[2] = result[1][3:6]-result[1][6:]
    d[3] = result[1][6:]
    if deriv == 1:
        return v, d
    dd = np.zeros((4, 3, 4, 3), float)
    aa = result[2][:3, :3]
    ab = result[2][:3, 3:6]
    ac = result[2][:3, 6:]
    ba = result[2][3:6, :3]
    bb = result[2][3:6, 3:6]
    bc = result[2][3:6, 6:]
    ca = result[2][6:, :3]
    cb = result[2][6:, 3:6]
    cc = result[2][6:, 6:]

    dd[0, :, 0, :] =   aa
    dd[0, :, 1, :] = - aa - ab
    dd[0, :, 2, :] =   ab - ac
    dd[0, :, 3, :] =   ac

    dd[1, :, 0, :] = - aa - ba
    dd[1, :, 1, :] =   aa + ba + ab + bb
    dd[1, :, 2, :] = - ab - bb + ac + bc
    dd[1, :, 3, :] = - ac - bc

    dd[2, :, 0, :] =   ba - ca
    dd[2, :, 1, :] = - ba + ca - bb + cb
    dd[2, :, 2, :] =   bb - cb - bc + cc
    dd[2, :, 3, :] =   bc - cc

    dd[3, :, 0, :] =   ca
    dd[3, :, 1, :] = - ca - cb
    dd[3, :, 2, :] =   cb - cc
    dd[3, :, 3, :] =   cc
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")


def _opbend_transform(rs, fn_low, deriv):
    a = rs[1] - rs[0]
    b = rs[2] - rs[0]
    c = rs[3] - rs[0]
    result = fn_low(a, b, c, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = np.zeros((4, 3), float)
    d[0] = -result[1][:3]-result[1][3:6]-result[1][6:]
    d[1] = result[1][:3]
    d[2] = result[1][3:6]
    d[3] = result[1][6:]
    if deriv == 1:
        return v, d
    dd = np.zeros((4, 3, 4, 3), float)
    aa = result[2][:3, :3]
    ab = result[2][:3, 3:6]
    ac = result[2][:3, 6:]
    ba = result[2][3:6, :3]
    bb = result[2][3:6, 3:6]
    bc = result[2][3:6, 6:]
    ca = result[2][6:, :3]
    cb = result[2][6:, 3:6]
    cc = result[2][6:, 6:]

    dd[0, :, 0, :] =   aa + ab + ac + ba + bb + bc + ca + cb + cc
    dd[0, :, 1, :] = - aa - ba - ca
    dd[0, :, 2, :] = - ab - bb - cb
    dd[0, :, 3, :] = - ac - bc - cc

    dd[1, :, 0, :] = - aa - ab - ac
    dd[1, :, 1, :] =   aa
    dd[1, :, 2, :] =   ab
    dd[1, :, 3, :] =   ac

    dd[2, :, 0, :] = - ba - bb - bc
    dd[2, :, 1, :] =   ba
    dd[2, :, 2, :] =   bb
    dd[2, :, 3, :] =   bc

    dd[3, :, 0, :] = - ca - cb - cc
    dd[3, :, 1, :] =   ca
    dd[3, :, 2, :] =   cb
    dd[3, :, 3, :] =   cc
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")


def _opbend_transform_mean(rs, fn_low, deriv=0):
    """Compute the mean of the 3 opbends
    """
    v = 0.0
    d = np.zeros((4,3), float)
    dd = np.zeros((4,3,4,3), float)
    #loop over the 3 cyclic permutations
    for p in np.array([[0,1,2], [2,0,1], [1,2,0]]):
        opbend = _opbend_transform([rs[p[0]], rs[p[1]], rs[p[2]], rs[3]], fn_low, deriv)
        v += opbend[0]/3
        index0 = np.where(p==0)[0][0] #index0 is the index of the 0th atom (rs[0])
        index1 = np.where(p==1)[0][0]
        index2 = np.where(p==2)[0][0]
        index3 = 3
        if deriv>0:
            d[0] += opbend[1][index0]/3
            d[1] += opbend[1][index1]/3
            d[2] += opbend[1][index2]/3
            d[3] += opbend[1][index3]/3
        if deriv>1:
            dd[0, :, 0, :] += opbend[2][index0, :, index0, :]/3
            dd[0, :, 1, :] += opbend[2][index0, :, index1, :]/3
            dd[0, :, 2, :] += opbend[2][index0, :, index2, :]/3
            dd[0, :, 3, :] += opbend[2][index0, :, index3, :]/3

            dd[1, :, 0, :] += opbend[2][index1, :, index0, :]/3
            dd[1, :, 1, :] += opbend[2][index1, :, index1, :]/3
            dd[1, :, 2, :] += opbend[2][index1, :, index2, :]/3
            dd[1, :, 3, :] += opbend[2][index1, :, index3, :]/3

            dd[2, :, 0, :] += opbend[2][index2, :, index0, :]/3
            dd[2, :, 1, :] += opbend[2][index2, :, index1, :]/3
            dd[2, :, 2, :] += opbend[2][index2, :, index2, :]/3
            dd[2, :, 3, :] += opbend[2][index2, :, index3, :]/3

            dd[3, :, 0, :] += opbend[2][index3, :, index0, :]/3
            dd[3, :, 1, :] += opbend[2][index3, :, index1, :]/3
            dd[3, :, 2, :] += opbend[2][index3, :, index2, :]/3
            dd[3, :, 3, :] += opbend[2][index3, :, index3, :]/3
    if deriv==0:
        return v,
    elif deriv==1:
        return v, d
    elif deriv==2:
        return v, d, dd
    else:
        raise ValueError("deriv must be 0, 1 or 2.")


#
# Low level Internal coordinate functions
#


def _bond_length_low(r, deriv):
    """Similar to bond_length, but with a relative vector"""
    r = Vector3(3, deriv, r, (0, 1, 2))
    d = r.norm()
    return d.results()


def _bend_cos_low(a, b, deriv):
    """Similar to bend_cos, but with relative vectors"""
    a = Vector3(6, deriv, a, (0, 1, 2))
    b = Vector3(6, deriv, b, (3, 4, 5))
    a /= a.norm()
    b /= b.norm()
    return dot(a, b).results()


def _bend_angle_low(a, b, deriv):
    """Similar to bend_angle, but with relative vectors"""
    result = _bend_cos_low(a, b, deriv)
    return _cos_to_angle(result, deriv)


def _dihed_cos_low(a, b, c, deriv):
    """Similar to dihed_cos, but with relative vectors"""
    a = Vector3(9, deriv, a, (0, 1, 2))
    b = Vector3(9, deriv, b, (3, 4, 5))
    c = Vector3(9, deriv, c, (6, 7, 8))
    b /= b.norm()
    tmp = b.copy()
    tmp *= dot(a, b)
    a -= tmp
    tmp = b.copy()
    tmp *= dot(c, b)
    c -= tmp
    a /= a.norm()
    c /= c.norm()
    return dot(a, c).results()


def _dihed_angle_low(av, bv, cv, deriv):
    """Similar to dihed_cos, but with relative vectors"""
    a = Vector3(9, deriv, av, (0, 1, 2))
    b = Vector3(9, deriv, bv, (3, 4, 5))
    c = Vector3(9, deriv, cv, (6, 7, 8))
    b /= b.norm()
    tmp = b.copy()
    tmp *= dot(a, b)
    a -= tmp
    tmp = b.copy()
    tmp *= dot(c, b)
    c -= tmp
    a /= a.norm()
    c /= c.norm()
    result = dot(a, c).results()
    # avoid trobles with the gradients by either using arccos or arcsin
    if abs(result[0]) < 0.5:
        # if the cosine is far away for -1 or +1, it is safe to take the arccos
        # and fix the sign of the angle.
        sign = 1-(np.linalg.det([av, bv, cv]) > 0)*2
        return _cos_to_angle(result, deriv, sign)
    else:
        # if the cosine is close to -1 or +1, it is better to compute the sine,
        # take the arcsin and fix the sign of the angle
        d = cross(b, a)
        side = (result[0] > 0)*2-1 # +1 means angle in range [-pi/2,pi/2]
        result = dot(d, c).results()
        return _sin_to_angle(result, deriv, side)


def _opdist_low(av, bv, cv, deriv):
    """Similar to opdist, but with relative vectors"""
    a = Vector3(9, deriv, av, (0, 1, 2))
    b = Vector3(9, deriv, bv, (3, 4, 5))
    c = Vector3(9, deriv, cv, (6, 7, 8))
    n  = cross(a, b)
    n /= n.norm()
    dist = dot(c, n)
    return dist.results()


def _opbend_cos_low(a, b, c, deriv):
    """Similar to opbend_cos, but with relative vectors"""
    a = Vector3(9, deriv, a, (0, 1, 2))
    b = Vector3(9, deriv, b, (3, 4, 5))
    c = Vector3(9, deriv, c, (6, 7, 8))
    n  = cross(a,b)
    n /= n.norm()
    c /= c.norm()
    temp = dot(n,c)
    result = temp.copy()
    result.v = np.sqrt(1.0-temp.v**2)
    if result.deriv > 0:
        result.d *= -temp.v
        result.d /= result.v
    if result.deriv > 1:
        result.dd *= -temp.v
        result.dd /= result.v
        temp2 = np.array([temp.d]).transpose()*temp.d
        temp2 /= result.v**3
        result.dd -= temp2
    return result.results()


def _opbend_angle_low(a, b, c, deriv=0):
    """Similar to opbend_angle, but with relative vectors"""
    result = _opbend_cos_low(a, b, c, deriv)
    sign = np.sign(np.linalg.det([a, b, c]))
    return _cos_to_angle(result, deriv, sign)


#
# Cosine and sine to angle conversion
#


def _cos_to_angle(result, deriv, sign=1):
    """Convert a cosine and its derivatives to an angle and its derivatives"""
    v = np.arccos(np.clip(result[0], -1, 1))
    if deriv == 0:
        return v*sign,
    if abs(result[0]) >= 1:
        factor1 = 0
    else:
        factor1 = -1.0/np.sqrt(1-result[0]**2)
    d = factor1*result[1]
    if deriv == 1:
        return v*sign, d*sign
    factor2 = result[0]*factor1**3
    dd = factor2*np.outer(result[1], result[1]) + factor1*result[2]
    if deriv == 2:
        return v*sign, d*sign, dd*sign
    raise ValueError("deriv must be 0, 1 or 2.")


def _sin_to_angle(result, deriv, side=1):
    """Convert a sine and its derivatives to an angle and its derivatives"""
    v = np.arcsin(np.clip(result[0], -1, 1))
    sign = side
    if sign == -1:
        if v < 0:
            offset = -np.pi
        else:
            offset = np.pi
    else:
        offset = 0.0
    if deriv == 0:
        return v*sign + offset,
    if abs(result[0]) >= 1:
        factor1 = 0
    else:
        factor1 = 1.0/np.sqrt(1-result[0]**2)
    d = factor1*result[1]
    if deriv == 1:
        return v*sign + offset, d*sign
    factor2 = result[0]*factor1**3
    dd = factor2*np.outer(result[1], result[1]) + factor1*result[2]
    if deriv == 2:
        return v*sign + offset, d*sign, dd*sign
    raise ValueError("deriv must be 0, 1 or 2.")
