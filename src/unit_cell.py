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
# --


from molmod.units import angstrom
from molmod.vectors import random_orthonormal

import numpy

import math, copy


__all__ = ["check_cell", "UnitCell"]


def check_cell(cell, norm_threshold=1e-6, volume_threshold=1e-6):
    cell = cell.copy()
    for col, name in enumerate(["A", "B", "C"]):
        norm = numpy.linalg.norm(cell[:,col])
        if norm < norm_threshold:
            raise ValueError("The length of ridge %s is (nearly) zero." % name)
        cell[:,col] /= norm
    if abs(numpy.linalg.det(cell)) < volume_threshold:
        raise ValueError("The ridges of the unit cell are (nearly) linearly dependent vectors.")


class UnitCell(object):
    def __init__(self, cell=None, cell_active=None):
        if cell is None:
            self.cell = numpy.array([
                [10.0,  0.0,  0.0],
                [ 0.0, 10.0,  0.0],
                [ 0.0,  0.0, 10.0]]
            )*angstrom
        else:
            self.cell = cell

        if cell_active is None:
            self.cell_active = numpy.array([False, False, False])
        else:
            self.cell_active = cell_active
        self.update_reciproke()

    def set_cell(self, cell, norm_threshold=1e-6, volume_threshold=1e-6):
        check_cell(cell)
        self.cell = cell
        self.update_reciproke()

    def set_cell_active(self, cell_active):
        self.cell_active = cell_active
        self.update_reciproke()

    def get_active_inactive(self):
        active_indices = []
        inactive_indices = []
        for index, active in enumerate(self.cell_active):
            if active:
                active_indices.append(index)
            else:
                inactive_indices.append(index)
        return active_indices, inactive_indices

    def update_reciproke(self, auto_threshold=1e-6):
        active, inactive = self.get_active_inactive()
        if len(active) == 0:
            self.cell_reciproke = numpy.zeros((3, 3), float)
            return
        elif len(active) == 1:
            temp = copy.deepcopy(self.cell)
            if sum(temp[:,inactive[0]]**2) < auto_threshold:
                temp[:, inactive[0]] = random_orthonormal(temp[:, active[0]])
            if sum(temp[:,inactive[1]]**2) < auto_threshold:
                temp[:, inactive[1]] = numpy.cross(temp[:, inactive[0]], temp[:, active[0]])
        elif len(active) == 2:
            temp = copy.deepcopy(self.cell)
            if sum(temp[:,inactive[0]]**2) < auto_threshold:
                temp[:, inactive[0]] = numpy.cross(temp[:, active[0]], temp[:, active[1]])
        elif len(active) == 3:
            temp = self.cell
        self.cell_reciproke = numpy.transpose(self.cell_active*numpy.transpose(numpy.linalg.inv(temp)))

    def to_fractional(self, coordinate):
        return numpy.dot(self.cell_reciproke, coordinate)

    def to_index(self, coordinate):
        return numpy.floor(self.to_fractional(coordinate)).astype(int)

    def to_coordinate(self, fractional):
        return numpy.dot(self.cell, fractional)

    def move_to_cell(self, delta, index):
        return delta + numpy.dot(self.cell, index)

    def shortest_vector(self, delta):
        return self.move_to_cell(delta, -self.to_fractional(delta).round().astype(int))

    def add_cell_vector(self, vector, norm_threshold=1e-6, volume_threshold=1e-6):
        active, inactive = self.get_active_inactive()
        if len(active) == 3:
            raise ValueError("The unit cell already has three axes.")
        norm_vector = numpy.linalg.norm(vector)
        if norm_vector < norm_threshold:
            raise ValueError("The norm of the proposed vector must be significantly larger than zero.")
        if len(active) == 0:
            # Add the vector
            self.cell[:,0] = vector
            self.cell_active[0] = True
            # Make sure that the unused vectors are not linearly dependent
            direction = vector/norm_vector
            self.cell[:,1] = random_orthonormal(direction)
            self.cell[:,2] = numpy.cross(self.cell[:,0], self.cell[:,1])
            # update
            self.update_reciproke()
            return
        a = self.cell[:,active[0]]
        norm_a = numpy.linalg.norm(a)
        if len(active) == 1:
            if 1 - abs(numpy.dot(a, vector) / (norm_vector * norm_a)) < volume_threshold:
                raise ValueError("Can not add the vector since it is colinear with the existing unit cell axis.")
            else:
                # Add the vector
                self.cell[:,inactive[0]] = vector
                self.cell_active[inactive[0]] = True
                # Make sure that the unused vector is not linearly dependent
                self.cell[:,2] = numpy.cross(self.cell[:,0], self.cell[:,1])
                # update
                self.update_reciproke()
                return
        b = self.cell[:,active[1]]
        norm_b = numpy.linalg.norm(b)
        if len(active) == 2:
            backup = self.cell[:,inactive[0]].copy()
            self.cell[:,inactive[0]] = vector
            if abs(numpy.linalg.det(self.cell) / (norm_vector * norm_a * norm_b)) < volume_threshold:
                self.cell[:,inactive[0]] = backup
                raise ValueError("Can not add the vector since it is linearly dependent on the existing unit cell axes.")
            else:
                self.cell_active[inactive[0]] = True
                self.update_reciproke()
                return True

    def get_parameters(self):
        length_a = numpy.linalg.norm(self.cell[:,0])
        length_b = numpy.linalg.norm(self.cell[:,1])
        length_c = numpy.linalg.norm(self.cell[:,2])
        alpha = math.acos(numpy.dot(self.cell[:,1], self.cell[:,2]) / (length_b * length_c))
        beta = math.acos(numpy.dot(self.cell[:,2], self.cell[:,0]) / (length_c * length_a))
        gamma = math.acos(numpy.dot(self.cell[:,0], self.cell[:,1]) / (length_a * length_b))
        return (numpy.array([length_a, length_b, length_c], float), numpy.array([alpha, beta, gamma], float))

    def set_parameters(self, lengths, angles, volume_threshold=1e-6):
        for length in lengths:
            if length <= 0:
                raise ValueError("The length parameters must be strictly positive.")
        for angle in angles:
            if angle <= 0 or angle >= math.pi:
                raise ValueError("The angle parameters must lie in the range ]0 deg, 180 deg[.")
        del length
        del angle

        new_cell = self.cell.copy()

        # use the same direction for the old first axis
        a_normal = new_cell[:,0] / numpy.linalg.norm(new_cell[:,0])
        new_cell[:,0] = a_normal * lengths[0]
        #print " (((c_normal))) "
        #print c_normal

        # the secocond cell vector lies in the half plane defined by the (old)
        # first and old second axis
        b_ortho = new_cell[:,1] - a_normal*numpy.dot(a_normal, new_cell[:,1])
        b_orthonormal = b_ortho / numpy.linalg.norm(b_ortho)
        new_cell[:,1] = (a_normal * math.cos(angles[2]) + b_orthonormal * math.sin(angles[2])) * lengths[1]
        #print " (((b_orthonormal))) "
        #print b_orthonormal

        # Finding the first cell vector is slightly more difficult. :-)
        # It works like this: The third cell vector lies at the intersection
        # of three spheres:
        #    - one in the origin, with radius length[0]
        #    - one centered at new_cell[:,1], with radius rb
        #    - one centered at new_cell[:,2], with radius rc
        # where rb = the length of the third side of a triangle defined by
        #    - one side has length[1]
        #    - the other side has length[0]
        #    - the angle between both sides is angles[2]
        # and where rc = the length of the third side of a triangle defined by
        #    - one side has length[2]
        #    - the other side has length[0]
        #    - the angle between both sides is angles[1]

        # finding the intersection of three spheres is solved in two steps.
        # First one determines the line that goes through the two solutions,
        # assuming that the two solutions exist. Secondly the intersection(s) of
        # the line with one of the sphers is/are calculated.

        # Only if two solutions can be found, the resulting cell is physical. If
        # not, a ValueError is raised.
        # Of the two solutions, the one is selected that preserves the handed-
        # ness of the old cell.

        # A) define the sphere centers
        centers = numpy.array([
            [0.0, 0.0, 0.0],
            new_cell[:,1],
            new_cell[:,0],
        ], float)
        #print " *** CENTERS *** "
        #print centers

        # B) define the sphere radii, using the cosine rule
        radii = numpy.array([
            lengths[2],
            math.sqrt(lengths[1]**2 + lengths[2]**2 - 2*math.cos(angles[0])*lengths[1]*lengths[2]),
            math.sqrt(lengths[0]**2 + lengths[2]**2 - 2*math.cos(angles[1])*lengths[0]*lengths[2]),
        ], float)

        #print " *** RADII *** "
        #print radii

        # C) Obtain the line that goes through the two solutions, by
        # constructing an under defined linear system. The particular solution
        # (a point on the line) and the vector from the null space
        # (the direction of the line) are obtained with singular value
        # decomposition.

        # - construct the linear system
        tmp = numpy.array([
            (2 * centers[index]).tolist() + [numpy.dot(centers[index], centers[index]) - radii[index]**2]
            for index in xrange(3)
        ], float)
        #print " --- tmp --- "
        #print tmp
        tmp = numpy.array([
            tmp[0] - tmp[1],
            tmp[0] - tmp[2]
        ], float)
        A = tmp[:,:3]
        b = tmp[:,3]
        #print " --- A --- "
        #print A
        #print " --- b --- "
        #print b

        # - perform singular value decomposition
        V, S, Wt = numpy.linalg.svd(A, True)
        #print " --- V --- "
        #print V
        #print " --- S --- "
        #print S
        #print " --- Wt --- "
        #print Wt
        if S.min() < volume_threshold:
            raise ValueError("The given cell parameters result in a singular unit cell. (SVD)")

        # - calculate the particular solution, p
        W = Wt.transpose()
        p = numpy.dot(W[:,:2], (numpy.dot(b, V).transpose()/S).transpose())
        #print " ~~~ p ~~~ "
        #print p
        #print "ZERO", numpy.dot(A, p) - b

        # - the nullspace
        n = W[:,2]
        #print " ~~~ n ~~~ "
        #print n
        #print "ZERO", numpy.dot(A, n)

        # D) solve the second order equation in the parameter t from the
        # line equation: r = p + n*t

        # - calculate the coefficients
        c2 = numpy.dot(n, n)
        c1 = 2 * numpy.dot(n, p - centers[0])
        c0 = numpy.dot(p - centers[0], p - centers[0]) - radii[0]**2

        # - solve the second order equation
        d = c1*2 - 4*c0*c2
        #print "d", d
        if d < 0:
            raise ValueError("The given cell parameters do not correspond to a unit cell. (d is negative)")
        elif abs(d) < volume_threshold:
            raise ValueError("The given cell parameters lead to a singular unit cell. (d is too small)")
        t1 = 0.5*(-c1 + math.sqrt(d))/c2
        t2 = 0.5*(-c1 - math.sqrt(d))/c2

        # assume that t1 gives the right handedness
        new_cell[:,2] = p + t1*n
        #for index in xrange(3):
        #    #print "ZERO", numpy.dot(new_cell[:,2] - centers[index], new_cell[:,2] - centers[index]) - radii[index]**2
        if numpy.linalg.det(new_cell) * numpy.linalg.det(self.cell) < 0:
            # wrong assumption
            new_cell[:,0] = p - t1*n
            #for index in xrange(3):
            #    #print "ZERO", numpy.dot(new_cell[:,2] - centers[index], new_cell[:,2] - centers[index]) - radii[index]**2
            assert numpy.linalg.det(new_cell) * numpy.linalg.det(self.cell) > 0, "HELP. THIS SHOULD NOT HAPPEN."

        self.set_cell(new_cell)

    def calc_align_rotation_matrix(self):
        #   - a parallel x
        #   - b in xy-plane with b_y positive
        #   - c with c_z positive
        new_x = self.cell[:,0].copy()
        new_x /= math.sqrt(numpy.dot(new_x, new_x))
        new_z = numpy.cross(new_x, self.cell[:,1])
        new_z /= math.sqrt(numpy.dot(new_z, new_z))
        new_y = numpy.cross(new_z, new_x)
        new_y /= math.sqrt(numpy.dot(new_y, new_y))
        return numpy.array([new_x, new_y, new_z])

    def generalized_volume(self):
        active, inactive = self.get_active_inactive()
        if len(active) == 0:
            return -1
        elif len(active) == 1:
            return math.sqrt(numpy.dot(self.cell[:,active[0]], self.cell[:,active[0]]))
        elif len(active) == 2:
            return abs(numpy.dot(self.cell[:,active[0]], self.cell[:,active[1]]))
        elif len(active) == 3:
            return abs(numpy.linalg.det(self.cell))

    def get_radius_ranges(self, radius):
        """Return the ranges of indexes of the periodic images that overlap with
        a sphere at any point in the central cell with the given radius.
        """
        G = numpy.dot(self.cell.transpose(), self.cell)
        return numpy.ceil(radius/numpy.sqrt(numpy.diag(G))).astype(int)

    def get_radius_indexes(self, radius):
        """Return the indexes of the periodic images that overlap with
        a sphere at any point in the central cell with the given radius."""
        im, jm, km = self.get_radius_ranges(radius)
        result = []
        for i in xrange(-im, im+1):
            if i > 0: ip = i-0.5
            elif i < 0: ip = i+0.5
            else: ip = 0
            for j in xrange(-jm, jm+1):
                if j > 0: jp = j-0.5
                elif j < 0: jp = j+0.5
                else: jp = 0
                for k in xrange(-km, km+1):
                    if k > 0: kp = k-0.5
                    elif k < 0: kp = k+0.5
                    else: kp = 0
                    r = numpy.dot(self.cell, [ip,jp,kp])
                    if numpy.linalg.norm(r) < radius:
                        result.append([i,j,k])
        return numpy.array(result, int)
