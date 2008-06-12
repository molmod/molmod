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



import numpy, math


__all__ = [
    "Base", "Translation", "Rotation", "Complete", "rotation_around_center",
    "CoincideError", "coincide",
]


class Base(object):
    def clear(self):
        raise NotImplementedError

    def from_matrix(self, matrix):
        raise NotImplementedError

    def get_matrix(self):
        raise NotImplementedError

    def get_inverse_matrix(self):
        raise NotImplementedError

    def invert(self):
        raise NotImplementedError

    def vector_apply(self, v):
        raise NotImplementedError

    def vector_apply_inverse(self, v):
        raise NotImplementedError

    def vector_apply_translation(self, v):
        raise NotImplementedError

    def matrix_apply_before(self, m):
        raise NotImplementedError

    def matrix_apply_inverse_before(self, m):
        raise NotImplementedError

    def matrix_apply_after(self, m):
        raise NotImplementedError

    def matrix_apply_inverse_after(self, m):
        raise NotImplementedError

    def apply_after(self, parent): # self -> parent AFTER self
        raise NotImplementedError

    def apply_inverse_after(self, parent): # self -> !parent AFTER self
        raise NotImplementedError

    def apply_before(self, child): # self -> self AFTER child
        raise NotImplementedError

    def apply_inverse_before(self, child): # self -> self AFTER !child
        raise NotImplementedError

    def compare(self, other, translation_threshold=1e-3, rotation_threshold=1e-3):
        raise NotImplementedError

    def assign_shallow(self, other):
        raise NotImplementedError


class Translation(Base):
    def __init__(self):
        self.t = numpy.zeros(3, float)

    def __str__(self):
        result = "TRANSLATION\n"
        for i in range(3):
            result += "% 10.7f\n" % self.t[i]
        return result[:-1]

    def clear(self):
        self.t[:] = 0

    def from_matrix(self, m):
        # check wether the t part is ok
        z = m[3, 0:3]
        numpy.power(z,2,z)
        assert max(z) < 1.0e-6, "The given matrix doesn't have correct translational part"
        assert m[3,3] == 1.0, "The lower right element of the given matrix must be 1.0."
        # get the translational part
        self.t = m[0:3, 3]

    def get_matrix(self):
        temp = numpy.identity(4, float)
        temp[0:3, 3] = self.t
        return temp

    def get_inverse_matrix(self):
        temp = numpy.identity(4, float)
        temp[0:3, 3] = -self.t
        return temp

    def invert(self):
        self.t *= -1

    def vector_apply(self, v):
        return v + self.t

    def vector_apply_inverse(self, v):
        return v - self.t

    def vector_apply_translation(self, v):
        return v + self.t

    def matrix_apply_before(self, m):
        return m

    def matrix_apply_inverse_before(self, m):
        return m

    def matrix_apply_after(self, m):
        return m

    def matrix_apply_inverse_after(self, m):
        return m

    def apply_after(self, parent): # self -> parent AFTER self
        self.t = parent.vector_apply(self.t)

    def apply_inverse_after(self, parent): # self -> !parent AFTER self
        self.t = parent.vector_apply_inverse(self.t)

    def apply_before(self, child): # self -> self AFTER child
        self.t = self.vector_apply(child.vector_apply(numpy.zeros(3, float)))

    def apply_inverse_before(self, child): # self -> self AFTER !child
        self.t = self.vector_apply(child.vector_apply_inverse(numpy.zeros(3, float)))

    def compare(self, other, translation_threshold=1e-3):
        return sum((self.t - other.t)**2) < translation_threshold

    def assign_shallow(self, other):
        if isinstance(other, Translation):
            self.t = other.t


class Rotation(Base):
    def __init__(self):
        self.r = numpy.identity(3, float)

    def __str__(self):
        result = "ROTATION\n"
        for i in range(3):
            result += "[ % 10.7f \t % 10.7f \t % 10.7f ]\n" % tuple(self.r[i])
        result += "det: %3.2f" % numpy.linalg.det(self.r)
        return result

    def clear(self):
        self.r[:] = 0
        self.r.ravel()[::4] = 1

    def from_matrix(self, m):
        self.r = m[0:3, 0:3]

    def get_matrix(self):
        temp = numpy.identity(4, float)
        temp[0:3, 0:3] = self.r
        return temp

    def get_inverse_matrix(self):
        temp = numpy.identity(4, float)
        temp[0:3, 0:3] = self.r.transpose()
        return temp

    def invert(self):
        self.r = self.r.transpose()

    def inversion_rotation(self):
        self.r *= -1

    def vector_apply(self, v):
        return numpy.dot(self.r, v)

    def vector_apply_inverse(self, v):
        return numpy.dot(self.r.transpose(), v)

    def vector_apply_translation(self, v):
        return v

    def matrix_apply_before(self, m):
        return numpy.dot(self.r, m)

    def matrix_apply_inverse_before(self, m):
        return numpy.dot(self.r.transpose(), m)

    def matrix_apply_after(self, m):
        return numpy.dot(m, self.r)

    def matrix_apply_inverse_after(self, m):
        return numpy.dot(m, self.r.transpose())

    def apply_after(self, parent): # self -> parent AFTER self
        self.r = parent.matrix_apply_before(self.r)

    def apply_inverse_after(self, parent): # self -> !parent AFTER self
        self.r = parent.matrix_apply_inverse_before(self.r)

    def apply_before(self, child): # self -> self AFTER child
        self.r = child.matrix_apply_after(self.r)

    def apply_inverse_before(self, child): # self -> self AFTER !child
        self.r = child.matrix_apply_inverse_after(self.r)

    def compare(self, other, rotation_threshold=1e-3):
        return sum((self.r - other.r).ravel()**2) < rotation_threshold

    def assign_shallow(self, other):
        if isinstance(other, Rotation):
            self.r = other.r

    def get_rotation_properties(self):
        # determine wether an inversion rotation has been applied
        invert = (numpy.linalg.det(self.r) < 0)
        factor = {True: -1, False: 1}[invert]
        # get the rotation data
        # trace(r) = 1+2*cos(angle)
        cos_angle = 0.5*(factor*numpy.trace(self.r) - 1)
        if cos_angle > 1: cos_angle = 1.0
        if cos_angle < -1: cos_angle = -1.0
        # the antisymmetric part of the non-diagonal vector tell us something
        # about sin(angle) and n.
        axis = 0.5*factor*numpy.array([-self.r[1, 2] + self.r[2, 1], self.r[0, 2] - self.r[2, 0], -self.r[0, 1] + self.r[1, 0]])
        sin_angle = math.sqrt(numpy.dot(axis, axis))
        # look for the best way to normalize the
        if (sin_angle == 0.0) and (cos_angle > 0):
            axis[2] = 1.0
        elif abs(sin_angle) < (1-cos_angle):
            for index in range(3):
                axis[index] = {True: -1, False: 1}[axis[index] < 0] * math.sqrt(abs((factor*self.r[index, index] - cos_angle) / (1 - cos_angle)))
        else:
            axis = axis / sin_angle

        # Finally calculate the angle:
        angle = math.atan2(sin_angle, cos_angle)
        return angle, axis, invert

    def set_rotation_properties(self, angle, axis, invert):
        norm = math.sqrt(numpy.dot(axis, axis))
        if norm > 0:
            x = axis[0] / norm
            y = axis[1] / norm
            z = axis[2] / norm
            c = math.cos(angle)
            s = math.sin(angle)
            self.r = (1-2*invert) * numpy.array([
                [x*x*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s],
                [x*y*(1-c)+z*s, y*y*(1-c)+c  , y*z*(1-c)-x*s],
                [x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c  ]
            ])
        else:
            self.r = numpy.identity(3) * (1-2*invert)


class Complete(Translation, Rotation):
    def __init__(self):
        self.t = numpy.zeros(3, float)
        self.r = numpy.identity(3, float)

    def __str__(self):
        result = "COMPLETE\n"
        for i in range(3):
            result += "[ % 10.7f \t % 10.7f \t % 10.7f ] \t % 10.7f \n" % (tuple(self.r[i]) + (self.t[i],))
        result += "det: %3.2f" % numpy.linalg.det(self.r)
        return result

    def clear(self):
        Translation.clear(self)
        Rotation.clear(self)

    def from_matrix(self, m):
        Rotation.from_matrix(self, m)
        Translation.from_matrix(self, m)

    def get_matrix(self):
        temp = Translation.get_matrix(self)
        temp[0:3, 0:3] = self.r
        return temp

    def get_inverse_matrix(self):
        invtrans = numpy.dot(-self.t, self.r)
        temp = Rotation.get_inverse_matrix(self)
        temp[0:3, 3] = invtrans
        return temp

    def invert(self):
        self.r = self.r.transpose()
        self.t = numpy.dot(self.r, -self.t)

    def vector_apply(self, v):
        return numpy.dot(self.r, v) + self.t

    def vector_apply_inverse(self, v):
        return numpy.dot(self.r.transpose(), v - self.t)

    def vector_apply_translation(self, v):
        return v + self.t

    def matrix_apply_before(self, m):
        return numpy.dot(self.r, m)

    def matrix_apply_inverse_before(self, m):
        return numpy.dot(self.r.transpose(), m)

    def matrix_apply_after(self, m):
        return numpy.dot(m, self.r)

    def matrix_apply_inverse_after(self, m):
        return numpy.dot(m, self.r.transpose())

    def apply_after(self, parent): # self -> parent AFTER self
        Translation.apply_after(self, parent)
        Rotation.apply_after(self, parent)

    def apply_inverse_after(self, parent): # self -> !parent AFTER self
        Translation.apply_inverse_after(self, parent)
        Rotation.apply_inverse_after(self, parent)

    def apply_before(self, child): # self -> self AFTER child
        Translation.apply_before(self, child)
        Rotation.apply_before(self, child)

    def apply_inverse_before(self, child): # self -> self AFTER !child
        Translation.apply_inverse_before(self, child)
        Rotation.apply_inverse_before(self, child)

    def compare(self, other, translation_threshold=1e-3, rotation_threshold=1e-3):
        return (
            sum((self.t - other.t)**2) < translation_threshold and
            sum((self.r - other.r).ravel()**2) < rotation_threshold
        )

    def assign_shallow(self, other):
        Translation.assign_shallow(self, other)
        Rotation.assign_shallow(self, other)


def rotation_around_center(center, angle, axis, invert=False):
    result = Complete()
    result.t = -center

    rotation = Rotation()
    rotation.set_rotation_properties(angle, axis, invert)
    result.apply_after(rotation)

    translation = Translation()
    translation.t = center
    result.apply_after(translation)

    return result


def random_rotation():
    from molmod.vectors import random_unit, trivial_orthonormal
    result = Rotation()
    # first generate a random unit vector, the new x-axis
    result.r[:,0] = random_unit(3)
    x = result.r[:,0]
    # generate a not so random y-axis and z-axis
    y = trivial_orthonormal(x)
    z = numpy.cross(x, y)
    # rotate y,z with about the x-axis by a random angle
    angle = numpy.random.uniform(0, 2*numpy.pi)
    result.r[:,1] = numpy.cos(angle)*y - numpy.sin(angle)*z
    result.r[:,2] = numpy.sin(angle)*y + numpy.cos(angle)*z
    return result


class CoincideError(Exception):
    pass

def coincide(ras, rbs, maxiter=100):
    """Compute the transformation that minimizes the RMSD between the points ras and rbs

    Both ras and rbs are Nx3 numpy arrays. Each row corresponds to a 3D
    coordinate and corresponding rows contain the points that are brought
    into overlap.
    """
    # This is a constrained least squares problem. The constraints are non-
    # linear, while the overlap equations are linear. The rbs are the atom
    # coordinates to be transformed, while the ras remain fixed in space.
    dm = [] # design_matrix
    ev = [] # expected_values
    # parameters are ordered as follows: R11, R12, R13, R21, R22, R23, R31, R32, R33, T1, T2, T3
    for ra, rb in zip(ras, rbs):
        # per atom pair, there are three equations, equaling x, y and z.
        dm.append(numpy.array([rb[0], rb[1], rb[2], 0, 0, 0, 0, 0, 0, 1, 0, 0]))
        dm.append(numpy.array([0, 0, 0, rb[0], rb[1], rb[2], 0, 0, 0, 0, 1, 0]))
        dm.append(numpy.array([0, 0, 0, 0, 0, 0, rb[0], rb[1], rb[2], 0, 0, 1]))
        ev.extend(ra)
    dm = numpy.array(dm)
    ev = numpy.array(ev)
    # We take the (ugly) approach with the normal equations to simplify the
    # treatment of the lagrange multipliers. An SVD based approach should
    # also be possible, but is definitely more complicated and total
    # overkill for the purpose of this function.
    Ap = numpy.dot(dm.transpose(), dm)
    Bp = numpy.dot(dm.transpose(), ev)
    Cp = numpy.dot(ev, ev)
    # contruct the terms for the lagrange multipliers
    # The matrices M and constants m are constructed in such a way that the
    # constraints reduce to the following simple form:
    # 0.5 x^T M x - m = 0
    # where x containts the first nign elements of the parameter vector
    # above. The derivative of this constraint is trivial.
    Ms = []
    ms = []
    for i1 in xrange(3):
        for i2 in xrange(i1,3):
            M = numpy.zeros((9,9), float)
            Ms.append(M)
            for k1 in xrange(3):
                M[k1*3+i1,k1*3+i2] += 1
                M[k1*3+i2,k1*3+i1] += 1
            ms.append(int(i1==i2))
    Ms = numpy.array(Ms)
    ms = numpy.array(ms)
    # Find the lagrange multiplieres (lambdas), starting from zero with the
    # Newton-Raphson scheme.
    from molmod.linalg import safe_solve, extended_solve
    #safe_solve = numpy.linalg.solve
    lambdas = numpy.zeros(6, float) # lagrange multipliers
    g = numpy.zeros(6, float) # constraint functions
    G = numpy.zeros((6,6), float) # derivatives of constraint functions
    # towards lagrange multipliers. Each row corresponds to one constraint
    # function while each column corresponds to a lagrange multiplier.
    K = numpy.zeros((9,6), float) # derivatives of the first nign elements
    # of the solution towards the lambdas. The rows correspond to components
    # of the solution, while the columns correspond to the lagrange
    # multipliers.
    tmp = numpy.zeros(12, float)
    A = numpy.zeros(Ap.shape, float)
    done = False

    for counter in xrange(maxiter):
        # Make a linear approximation of the lambdas that solve the constraints.
        # (this is the Newton-Raphson algo)
        A[:] = Ap
        for i in xrange(6):
            # loop over the lagrange multipliers
            A[:9,:9] += lambdas[i]*Ms[i]
        # solve locally with SVD, in case someone is aligning two linear
        # molecules.
        x, null_space = extended_solve(A, Bp, 1e-10)
        if len(null_space) > 0:
            # this is a nasty situation. we must use the vectors from the null
            # space to satisfy the constraints as good as possible
            alpha = opt_constraints(Ms, ms, x[:9].copy(), null_space[:9])
            x += numpy.dot(null_space, alpha)
        #x += numpy.dot(null_space, numpy.random.normal(0, 1, null_space.shape[1]))
        # xp is a shortcut
        xp = x[:9]
        # compute the derivative of x towards the lambdas
        for i in xrange(6):
            # loop over the lagrange multipliers
            tmp[:9] = numpy.dot(Ms[i], xp)
            tmp[9:] = 0
            tmp[:] = -safe_solve(A, tmp, 1e-10)
            K[:,i] = tmp[:9]
        # compute constraint functions and compute the derivative of these
        # constraint functions with respect to the lambdas
        for i in xrange(6):
            # loop over constraint functions
            Mxp = numpy.dot(Ms[i], xp)
            g[i] = 0.5*numpy.dot(xp, Mxp) - ms[i]
            G[i] = numpy.dot(Mxp, K)

        #print "norm(g)", numpy.linalg.norm(g)
        if numpy.linalg.norm(g) < 1e-9:
            #ssd = numpy.dot(x, numpy.dot(Ap, x) -2*Bp) + Cp
            #print "ssd", ssd
            done = True
            break

        step = -safe_solve(G, g, 1e-10)
        lambdas += step
        #low = 0.1
        #lambdas[lambdas < low] = low

    if not done:
        raise CoincideError("The constrained optimization could not be solved!")
    # wrap the solution in a fancy data structure
    trans = Complete()
    trans.r[:] = x[:9].reshape((3,3))
    trans.t[:] = x[9:]
    return trans


def opt_constraints(Ms, ms, x0, N):
    gradient_x = numpy.zeros(9, float)
    hessian_x = numpy.zeros((9, 9), float)
    while True:
        # Try a few initial guesses until an orthogonal system is found.
        alpha = numpy.random.normal(0, 1, N.shape[1])
        while True:
            x = x0 + numpy.dot(N, alpha)
            error = 0
            gradient_x[:] = 0
            hessian_x[:] = 0
            #print x
            for i in xrange(6):
                Mx = numpy.dot(Ms[i], x)
                xMxm = 0.5*numpy.dot(x, Mx) - ms[i]
                #print "xMxm", xMxm
                error += xMxm**2
                gradient_x[:] += 2*xMxm*Mx
                hessian_x[:] += 2*(numpy.outer(Mx, Mx) + xMxm*Ms[i])
            #print "error", error
            gradient_alpha = numpy.dot(gradient_x, N)
            hessian_alpha = numpy.dot(numpy.dot(N.transpose(), hessian_x), N)
            #print numpy.linalg.norm(gradient_alpha)
            if numpy.linalg.norm(gradient_alpha) < 1e-10:
                break
            #print "evals", numpy.linalg.eigvals(hessian_alpha)
            step = -0.5*numpy.linalg.solve(hessian_alpha, gradient_alpha)
            if numpy.dot(step, gradient_alpha) > 0:
                #print "SD"
                step = -0.1*gradient_alpha/numpy.linalg.norm(gradient_alpha)
            #print "step", step
            alpha += step
        #print "error", error
        if abs(error) < 1e-10:
            break
    return alpha

