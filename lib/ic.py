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

import numpy


__all__ = [
    "Scalar", "Vector3", "dot",
    "bond_length", "pair_distance",
    "bend_cos", "bend_angle",
    "dihed_cos", "dihed_angle",
]


class Scalar(object):
    """A scalar object with optional first and second order derivates

       Each input value to which the derivative is computed has its own index.
       The numerical value of the derivatives are stored in arrays self.d and
       self.dd. The value of the scalar itself if self.v
    """

    def __init__(self, size, deriv=0, value=0, index=None):
        """Initialize a scalar internal coordinate.

           Arguments:
             size -- The number of inputs on which this ic depends. e.g. a distance
                     depends on 6 Cartesian coordinates.
             deriv -- Consider up to deriv order derivatives. (max=2)
             value -- The initial value.
             index -- If this scalar is one of the input variables, this is its index.

           The scalar object supports several in place modifications.
        """
        self.deriv = deriv
        self.size = size
        self.v = value
        if deriv > 0:
            self.d = numpy.zeros(size, float)
            if index is not None:
                self.d[index] = 1
        if deriv > 1:
            self.dd = numpy.zeros((size,size), float)
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

    def __isub__(self, other):
        if self.deriv > 1: self.dd -= other.dd
        if self.deriv > 0: self.d -= other.d
        self.v -= other.v
        return self

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
                tmp = numpy.outer(self.d, other.d)
                self.dd += tmp
                self.dd += tmp.transpose()
            if self.deriv > 0:
                self.d *= other.v
                self.d += self.v*other.d
            self.v *= other.v
        else:
            raise TypeError("Second argument must be float, int or Scalar")
        return self

    def __idiv__(self, other):
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
                tmp = numpy.outer(self.d, other.d)
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
            self.dd[:] = tmp*(2*self.v*numpy.outer(self.d, self.d) - self.dd)
        if self.deriv > 0:
            self.d[:] = -tmp*self.d[:]


class Vector3(object):
    """A Three dimensional vector with optional first and second order derivatives.

       This object is nothing more than a tier for three Scalar objects.
    """

    def __init__(self, size, deriv=0, values=(0,0,0), indexes=(None,None,None)):
        """Initialize a Vector3 object

           Arguments:
             size -- The number of inputs on which this ic depends. e.g. a distance
                     depends on 6 Cartesian coordinates.
             deriv -- Consider up to deriv order derivatives. (max=2)
             values -- The initial values.
             indexes -- If this vector is one of the input variables, these are
                        the indexes of the components.
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

    def __idiv__(self, other):
        self.x /= other
        self.y /= other
        self.z /= other
        return self

    def norm(self):
        """Return a Scalar object with the norm of this vector"""
        result = Scalar(self.size, self.deriv)
        result.v = numpy.sqrt(self.x.v**2 + self.y.v**2 + self.z.v**2)
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
            result.dd += (1 - self.x.v**2/denom)*numpy.outer(self.x.d, self.x.d)
            result.dd += (1 - self.y.v**2/denom)*numpy.outer(self.y.d, self.y.d)
            result.dd += (1 - self.z.v**2/denom)*numpy.outer(self.z.d, self.z.d)
            tmp = -self.x.v*self.y.v/denom*numpy.outer(self.x.d, self.y.d)
            result.dd += tmp+tmp.transpose()
            tmp = -self.y.v*self.z.v/denom*numpy.outer(self.y.d, self.z.d)
            result.dd += tmp+tmp.transpose()
            tmp = -self.z.v*self.x.v/denom*numpy.outer(self.z.d, self.x.d)
            result.dd += tmp+tmp.transpose()
            result.dd /= result.v
        return result

    #def norm2(self):
    #    """Return a Scalar object with the square of the norm of this vector"""
    #    result = Scalar(self.size, self.deriv)
    #    result.v = self.x.v**2 + self.y.v**2 + self.z.v**2
    #    if self.deriv > 0:
    #        result.d += self.x.v*self.x.d
    #        result.d += self.y.v*self.y.d
    #        result.d += self.z.v*self.z.d
    #        result.d *= 2
    #    if self.deriv > 1:
    #        result.dd += self.x.v*self.x.dd
    #        result.dd += self.y.v*self.y.dd
    #        result.dd += self.z.v*self.z.dd
    #        result.dd += numpy.outer(self.x.d, self.x.d)
    #        result.dd += numpy.outer(self.y.d, self.y.d)
    #        result.dd += numpy.outer(self.z.d, self.z.d)
    #        tmp = denom*numpy.outer(self.x.d, self.y.d)
    #        result.dd += tmp+tmp.transpose()
    #        tmp = denom*numpy.outer(self.y.d, self.z.d)
    #        result.dd += tmp+tmp.transpose()
    #        tmp = denom*numpy.outer(self.z.d, self.x.d)
    #        result.dd += tmp+tmp.transpose()
    #        result.dd *= 2
    #    return result


def dot(r1, r2):
    """Compute the dot product between two Vector3 Objects (Returns a Scalar)"""
    if r1.size != r2.size:
        raise ValueError("Both arguments must have the same input size.")
    if r1.deriv != r2.deriv:
        raise ValueError("Both arguments must have the same deriv.")
    result = Scalar(r1.size, r1.deriv)
    result.v = r1.x.v*r2.x.v + r1.y.v*r2.y.v + r1.z.v*r2.z.v
    if result.deriv > 0:
        result.d[:] += r1.x.v*r2.x.d
        result.d[:] += r2.x.v*r1.x.d
        result.d[:] += r1.y.v*r2.y.d
        result.d[:] += r2.y.v*r1.y.d
        result.d[:] += r1.z.v*r2.z.d
        result.d[:] += r2.z.v*r1.z.d
    if result.deriv > 1:
        result.dd[:] += numpy.outer(r1.x.d, r2.x.d)
        result.dd[:] += numpy.outer(r2.x.d, r1.x.d)
        result.dd[:] += numpy.outer(r1.y.d, r2.y.d)
        result.dd[:] += numpy.outer(r2.y.d, r1.y.d)
        result.dd[:] += numpy.outer(r1.z.d, r2.z.d)
        result.dd[:] += numpy.outer(r2.z.d, r1.z.d)
        result.dd[:] += r1.x.v*r2.x.dd
        result.dd[:] += r2.x.v*r1.x.dd
        result.dd[:] += r1.y.v*r2.y.dd
        result.dd[:] += r2.y.v*r1.y.dd
        result.dd[:] += r1.z.v*r2.z.dd
        result.dd[:] += r2.z.v*r1.z.dd
    return result


def bond_length(r0, r1, deriv=0):
    """Compute the distance between the two points r0 and r1

       Arguments:
         r0  --  a numpy array with three elements
         r1  --  a numpy array with three elements
         derive=0  --  the derivatives to be computed (0, 1 or 2)
    """
    r = r0 - r1
    result = _bond_length_low(r, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = numpy.zeros((2,3), float)
    d[0] = result[1]
    d[1] = -result[1]
    if deriv == 1:
        return v, d
    dd = numpy.zeros((2,3,2,3), float)
    dd[0,:,0,:] = result[2]
    dd[1,:,1,:] = result[2]
    dd[0,:,1,:] = -result[2]
    dd[1,:,0,:] = -result[2]
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")

pair_distance = bond_length

def _bond_length_low(r, deriv):
    """Similar to bond_length, but with a relative vector"""
    r = Vector3(3, deriv, r, (0,1,2))
    d = r.norm()
    return d.results()


def bend_cos(r0, r1, r2, deriv=0):
    """Compute the cosine of the angle between the vectors r0-r1 and r2-r1

       Arguments:
         r0  --  a numpy array with three elements
         r1  --  a numpy array with three elements
         r2  --  a numpy array with three elements
         derive=0  --  the derivatives to be computed (0, 1 or 2)
    """
    a = r0 - r1
    b = r2 - r1
    result = _bend_cos_low(a, b, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = numpy.zeros((3,3), float)
    d[0] = result[1][:3]
    d[1] = -result[1][:3]-result[1][3:]
    d[2] = result[1][3:]
    if deriv == 1:
        return v, d
    dd = numpy.zeros((3,3,3,3), float)
    aa = result[2][:3,:3]
    ab = result[2][:3,3:]
    ba = result[2][3:,:3]
    bb = result[2][3:,3:]
    dd[0,:,0,:] =   aa
    dd[0,:,1,:] = - aa - ab
    dd[0,:,2,:] =   ab
    dd[1,:,0,:] = - aa - ba
    dd[1,:,1,:] =   aa + ba + ab + bb
    dd[1,:,2,:] = - ab - bb
    dd[2,:,0,:] =   ba
    dd[2,:,1,:] = - ba - bb
    dd[2,:,2,:] =   bb
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")

def _bend_cos_low(a, b, deriv):
    """Similar to bend_cos, but with relative vectors"""
    a = Vector3(6, deriv, a, (0,1,2))
    b = Vector3(6, deriv, b, (3,4,5))
    a /= a.norm()
    b /= b.norm()
    return dot(a,b).results()

def _cos_to_angle(result, deriv, sign=1):
    """Convert a cosine and its derivatives to an angle"""
    v = numpy.arccos(numpy.clip(result[0],-1,1))
    if deriv == 0:
        return v*sign,
    if abs(result[0]) >= 1:
        factor1 = 0
    else:
        factor1 = -1.0/numpy.sqrt(1-result[0]**2)
    d = factor1*result[1]
    if deriv == 1:
        return v*sign, d*sign
    factor2 = result[0]*factor1**3
    size = result[2].shape[0]
    dd = (
        factor2*numpy.outer(result[1].ravel(), result[1].ravel()).reshape((size,3,size,3)) +
        factor1*result[2]
    )
    if deriv == 2:
        return v*sign, d*sign, dd*sign
    raise ValueError("deriv must be 0, 1 or 2.")

def bend_angle(r0, r1, r2, deriv=0):
    """Compute the angle between the vectors r0-r1 and r2-r1

       Arguments:
         r0  --  a numpy array with three elements
         r1  --  a numpy array with three elements
         r2  --  a numpy array with three elements
         derive=0  --  the derivatives to be computed (0, 1 or 2)
    """
    result = bend_cos(r0, r1, r2, deriv)
    return _cos_to_angle(result, deriv)


def dihed_cos(r0, r1, r2, r3, deriv=0):
    """Compute the cosine of the angle between the planes r0,r1,r2 and r1,r2,r3

       Arguments:
         r0  --  a numpy array with three elements
         r1  --  a numpy array with three elements
         r2  --  a numpy array with three elements
         r3  --  a numpy array with three elements
         derive=0  --  the derivatives to be computed (0, 1 or 2)
    """
    a = r0 - r1
    b = r2 - r1
    c = r3 - r2
    result = _dihed_cos_low(a, b, c, deriv)
    v = result[0]
    if deriv == 0:
        return v,
    d = numpy.zeros((4,3), float)
    d[0] = result[1][:3]
    d[1] = -result[1][:3]-result[1][3:6]
    d[2] = result[1][3:6]-result[1][6:]
    d[3] = result[1][6:]
    if deriv == 1:
        return v, d
    dd = numpy.zeros((4,3,4,3), float)
    aa = result[2][:3,:3]
    ab = result[2][:3,3:6]
    ac = result[2][:3,6:]
    ba = result[2][3:6,:3]
    bb = result[2][3:6,3:6]
    bc = result[2][3:6,6:]
    ca = result[2][6:,:3]
    cb = result[2][6:,3:6]
    cc = result[2][6:,6:]

    dd[0,:,0,:] =   aa
    dd[0,:,1,:] = - aa - ab
    dd[0,:,2,:] =   ab - ac
    dd[0,:,3,:] =   ac

    dd[1,:,0,:] = - aa - ba
    dd[1,:,1,:] =   aa + ba + ab + bb
    dd[1,:,2,:] = - ab - bb + ac + bc
    dd[1,:,3,:] = - ac - bc

    dd[2,:,0,:] =   ba - ca
    dd[2,:,1,:] = - ba + ca - bb + cb
    dd[2,:,2,:] =   bb - cb - bc + cc
    dd[2,:,3,:] =   bc - cc

    dd[3,:,0,:] =   ca
    dd[3,:,1,:] = - ca - cb
    dd[3,:,2,:] =   cb - cc
    dd[3,:,3,:] =   cc
    if deriv == 2:
        return v, d, dd
    raise ValueError("deriv must be 0, 1 or 2.")

def _dihed_cos_low(a, b, c, deriv):
    """Similar to dihed_cos, but with relative vectors"""
    a = Vector3(9, deriv, a, (0,1,2))
    b = Vector3(9, deriv, b, (3,4,5))
    c = Vector3(9, deriv, c, (6,7,8))
    b /= b.norm()
    tmp = b.copy()
    tmp *= dot(a,b)
    a -= tmp
    tmp = b.copy()
    tmp *= dot(c,b)
    c -= tmp
    a /= a.norm()
    c /= c.norm()
    return dot(a,c).results()

def dihed_angle(r0, r1, r2, r3, deriv=0):
    """Compute the angle between the planes r0,r1,r2 and r1,r2,r3

       The sign convention corresponds to the IUPAC definition of the torsion
       angle: http://dx.doi.org/10.1351/goldbook.T06406

       Arguments:
         r0  --  a numpy array with three elements
         r1  --  a numpy array with three elements
         r2  --  a numpy array with three elements
         r3  --  a numpy array with three elements
         derive=0  --  the derivatives to be computed (0, 1 or 2)
    """
    result = dihed_cos(r0, r1, r2, r3, deriv)
    a = r0 - r1
    b = r2 - r1
    c = r3 - r2
    sign = 1-(numpy.linalg.det([a,b,c]) > 0)*2
    return _cos_to_angle(result, deriv, sign)


