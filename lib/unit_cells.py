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
"""Data structure & tools to work with periodic systems"""


from molmod.units import angstrom
from molmod.vectors import random_orthonormal
from molmod.utils import cached, ReadOnly

import numpy


__all__ = ["UnitCell"]


class UnitCell(ReadOnly):
    """Extensible representation of a unit cell.

       Most attributes of the UnitCell object are treated as constants. If you
       want to modify the unit cell, just create a modified UnitCell
       object. This facilitates the caching of derived quantities such as the
       distance matrices, while it imposes a cleaner coding style without
       a signifacant computational overhead.
    """
    eps = 1e-6 # small positive number, below this value is approximately zero

    def __init__(self, matrix=None, active=None):
        if matrix is None:
            matrix = numpy.array([
                [10.0,  0.0,  0.0],
                [ 0.0, 10.0,  0.0],
                [ 0.0,  0.0, 10.0]]
            )*angstrom
        if active is None:
            active = numpy.array([True, True, True])
        ReadOnly.__init__(self)
        self._init_attributes({"matrix": matrix, "active": active}, {})
        # sanity checks for the unit cell
        for col, name in enumerate(["a", "b", "c"]):
            if self.active[col]:
                norm = numpy.linalg.norm(matrix[:, col])
                if norm < self.eps:
                    raise ValueError("The length of ridge %s is (nearly) zero." % name)
        if abs(self.generalized_volume) < self.eps:
            raise ValueError("The ridges of the unit cell are (nearly) linearly dependent vectors.")

    def __div__(self, i):
        if not isinstance(i, int) or i <= 0:
            raise ValueError("Can only divide a unit cell by a strictly positive integer.")
        return UnitCell(self.matrix/i, self.active)

    @classmethod
    def from_parameters3(cls, lengths, angles):
        """Construct a 3D unit cell with the given parameters

           The a vector is always parallel with the x-axis and they point in the
           same direction. The b vector is always in the xy plane and points
           towards the positive y-direction. The c vector points towards the
           positive z-direction.
        """
        for length in lengths:
            if length <= 0:
                raise ValueError("The length parameters must be strictly positive.")
        for angle in angles:
            if angle <= 0 or angle >= numpy.pi:
                raise ValueError("The angle parameters must lie in the range ]0 deg, 180 deg[.")
        del length
        del angle

        matrix = numpy.zeros((3, 3), float)

        # first cell vector
        matrix[0, 0] = lengths[0]

        # second cell vector
        matrix[0, 1] = numpy.cos(angles[2])*lengths[1]
        matrix[1, 1] = numpy.sin(angles[2])*lengths[1]

        # Finding the third cell vector is slightly more difficult. :-)
        # It works like this:
        # The dot products of a with c, b with c and c with c are known. the
        # vector a has only an x component, b has no z component. This results
        # in the following equations:
        u_a = lengths[0]*lengths[2]*numpy.cos(angles[1])
        u_b = lengths[1]*lengths[2]*numpy.cos(angles[0])
        matrix[0, 2] = u_a/matrix[0, 0]
        matrix[1, 2] = (u_b - matrix[0, 1]*matrix[0, 2])/matrix[1, 1]
        u_c = lengths[2]**2 - matrix[0, 2]**2 - matrix[1, 2]**2
        if u_c < 0:
            raise ValueError("The given cell parameters do not correspond to a unit cell.")
        matrix[2, 2] = numpy.sqrt(u_c)

        active = numpy.ones(3, bool)
        return cls(matrix, active)

    @cached
    def generalized_volume(self):
        """The volume of the unit cell

           The actual definition of the volume depends on the number of active
           directions:
             num_active == 0  --  always -1
             num_active == 1  --  length of the cell vector
             num_active == 2  --  surface of the parallellogram
             num_active == 3  --  volume of the parallelepiped
        """
        active = self.active_inactive[0]
        if len(active) == 0:
            return -1
        elif len(active) == 1:
            return numpy.linalg.norm(self.matrix[:, active[0]])
        elif len(active) == 2:
            return numpy.linalg.norm(numpy.cross(self.matrix[:, active[0]], self.matrix[:, active[1]]))
        elif len(active) == 3:
            return abs(numpy.linalg.det(self.matrix))

    @cached
    def active_inactive(self):
        """The indexes of the active and the inactive cell vectors"""
        active_indices = []
        inactive_indices = []
        for index, active in enumerate(self.active):
            if active:
                active_indices.append(index)
            else:
                inactive_indices.append(index)
        return active_indices, inactive_indices

    @cached
    def reciprocal(self):
        """The reciprocal of the unit cell

           In case of a three-dimensional periodic system, this is trivially the
           transpose of the inverse of the cell matrix. This means that each
           column of the matrix corresponds to a reciprocal cell vector. In case
           of lower-dimensional periodicity, the inactive columns are zero.
        """
        active, inactive = self.active_inactive
        if len(active) == 0:
            return numpy.zeros((3, 3), float)
        elif len(active) == 1:
            temp = self.matrix.copy()
            if numpy.linalg.norm(temp[:, inactive[0]]) < self.eps:
                temp[:, inactive[0]] = random_orthonormal(temp[:, active[0]])
            if numpy.linalg.norm(temp[:, inactive[1]]) < self.eps:
                temp[:, inactive[1]] = numpy.cross(temp[:, inactive[0]], temp[:, active[0]])
        elif len(active) == 2:
            temp = self.matrix.copy()
            if numpy.linalg.norm(temp[:, inactive[0]]) < self.eps:
                temp[:, inactive[0]] = numpy.cross(temp[:, active[0]], temp[:, active[1]])
        elif len(active) == 3:
            temp = self.matrix
        return self.active*numpy.transpose(numpy.linalg.inv(temp))

    @cached
    def parameters(self):
        """The cell parameters (lengths and angles)"""
        length_a = numpy.linalg.norm(self.matrix[:, 0])
        length_b = numpy.linalg.norm(self.matrix[:, 1])
        length_c = numpy.linalg.norm(self.matrix[:, 2])
        alpha = numpy.arccos(numpy.dot(self.matrix[:, 1], self.matrix[:, 2]) / (length_b * length_c))
        beta = numpy.arccos(numpy.dot(self.matrix[:, 2], self.matrix[:, 0]) / (length_c * length_a))
        gamma = numpy.arccos(numpy.dot(self.matrix[:, 0], self.matrix[:, 1]) / (length_a * length_b))
        return (numpy.array([length_a, length_b, length_c], float), numpy.array([alpha, beta, gamma], float))

    @cached
    def alignment_a(self):
        """Computes the rotation matrix that aligns the unit cell with the
           Cartesian axes, starting with cell vector a

              a parallel to x
              b in xy-plane with b_y positive
              c with c_z positive
        """
        from molmod.transformations import Rotation
        new_x = self.matrix[:, 0].copy()
        new_x /= numpy.linalg.norm(new_x)
        new_z = numpy.cross(new_x, self.matrix[:, 1])
        new_z /= numpy.linalg.norm(new_z)
        new_y = numpy.cross(new_z, new_x)
        new_y /= numpy.linalg.norm(new_y)
        return Rotation(numpy.array([new_x, new_y, new_z]))

    @cached
    def alignment_c(self):
        """Computes the rotation matrix that aligns the unit cell with the
           Cartesian axes, starting with cell vector c

              c parallel to z
              b in zy-plane with b_y positive
              a with a_x positive
        """
        from molmod.transformations import Rotation
        new_z = self.matrix[:, 2].copy()
        new_z /= numpy.linalg.norm(new_z)
        new_x = numpy.cross(self.matrix[:, 1], new_z)
        new_x /= numpy.linalg.norm(new_x)
        new_y = numpy.cross(new_z, new_x)
        new_y /= numpy.linalg.norm(new_y)
        return Rotation(numpy.array([new_x, new_y, new_z]))

    @cached
    def spacings(self):
        """Computes the distances between neighboring crystal planes"""
        return (self.reciprocal**2).sum(axis=0)**(-0.5)

    def to_fractional(self, cartesian):
        """Convert Cartesian to fractional coordinates

           Argument:
             cartesian  --  Can be a numpy array with shape (3, ) or with shape
                            (N, 3).

           The return value has the same shape as the argument. This function is
           the inverse of to_cartesian.
        """
        return numpy.dot(cartesian, self.reciprocal)

    def to_cartesian(self, fractional):
        """Converts fractional to Cartesian coordinates

           Argument:
             fractional  --  Can be a numpy array with shape (3, ) or with shape
                             (N, 3).

           The return value has the same shape as the argument. This function is
           the inverse of to_fractional.
        """
        return numpy.dot(fractional, self.matrix.transpose())

    def shortest_vector(self, delta):
        """Compute the shortest vector between two points under periodic
           boundary conditions.

           Argument:
             delta  --  the relative vector between two points
           Returns:
             The shortest relative vector in this unit cell.
        """
        fractional = self.to_fractional(delta)
        fractional -= fractional.round()
        return self.to_cartesian(fractional)

    def add_cell_vector(self, vector):
        """Returns a new unit cell with an additional cell vector"""
        act = self.active_inactive[0]
        if len(act) == 3:
            raise ValueError("The unit cell already has three active cell vectors.")
        matrix = numpy.zeros((3, 3), float)
        active = numpy.zeros(3, bool)
        if len(act) == 0:
            # Add the new vector
            matrix[:, 0] = vector
            active[0] = True
            return UnitCell(matrix, active)

        a = self.matrix[:, act[0]]
        matrix[:, 0] = a
        active[0] = True
        if len(act) == 1:
            # Add the new vector
            matrix[:, 1] = vector
            active[1] = True
            return UnitCell(matrix, active)

        b = self.matrix[:, act[1]]
        matrix[:, 1] = b
        active[1] = True
        if len(act) == 2:
            # Add the new vector
            matrix[:, 2] = vector
            active[2] = True
            return UnitCell(matrix, active)

    def get_radius_ranges(self, radius):
        """Return the ranges of indexes of the interacting neighboring unit cells

           Interacting neigboring unit cells have at least one point in their
           box volume that has a distance smaller or equal than radius to at
           least one point in the central cell. This concept is of importance
           when computing pair wise long-range interactions in periodic systems.
        """
        result = numpy.ceil(radius/self.spacings).astype(int)
        result[True^self.active] = 0
        return result

    def get_radius_indexes(self, radius):
        """Return the indexes of the interacting neighboring unit cells

           Interacting neigboring unit cells have at least one point in their
           box volume that has a distance smaller or equal than radius to at
           least one point in the central cell. This concept is of importance
           when computing pair wise long-range interactions in periodic systems.
        """

        def iter_indexes():
            """Iterate over all unit cell indexes in the range"""
            ranges = self.get_radius_ranges(radius)
            for i0 in xrange(-ranges[0], ranges[0]+1):
                for i1 in xrange(-ranges[1], ranges[1]+1):
                    for i2 in xrange(0, ranges[2]+1):
                        yield numpy.array([i0, i1, i2])

        flag_combinations = numpy.array([
            [ 1,  1,  1], [ 1,  1,  0], [ 1,  1, -1],
            [ 1,  0,  1], [ 1,  0,  0], [ 1,  0, -1],
            [ 1, -1,  1], [ 1, -1,  0], [ 1, -1, -1],
            [ 0,  1,  1], [ 0,  1,  0], [ 0,  1, -1],
            [ 0,  0,  1], [ 0,  0,  0], [ 0,  0, -1],
            [ 0, -1,  1], [ 0, -1,  0], [ 0, -1, -1],
            [-1,  1,  1], [-1,  1,  0], [-1,  1, -1],
            [-1,  0,  1], [-1,  0,  0], [-1,  0, -1],
            [-1, -1,  1], [-1, -1,  0], [-1, -1, -1],
        ])
        sign_combinations = numpy.array([
            [ 1,  1,  1], [ 1,  1, -1], [ 1, -1,  1], [ 1, -1, -1],
            [-1,  1,  1], [-1,  1, -1], [-1, -1,  1], [-1, -1, -1],
        ])

        def close_enough(center):
            """Test if two cell volumes could have interactions within radius

               Arguments:
                 center  --  the relative vector between the two cell volumes
            """
            # The problem can be simplified to finding the shortest distance
            # from the origin to a point in a double unit cell box centered
            # around the relative vector between the two unit cell boxes under
            # scrunity. Therefore the relative vector is called 'center' in
            # this routine.
            #
            # This is a constrained optimization problem that would in principle
            # be solved with the active set algorithm. We prefer to keep the
            # implementation simple and will try all reasonable combinations of
            # constraints.
            #
            # There are 6 constraints, i.e. the six faces of the box, which
            # naively leads to 2**6 possible combinations of turning constraints
            # of and on. Constrainging the solution to two oposite faces at the
            # same time does not make sense. Consequently, there are only 27
            # reasonable combinations of constraints.
            #
            # For each combination, only those that have a solution in, or on
            # the edge of the box are considered. From the latter, we test if
            # it is within the radius.

            for flags in flag_combinations:
                corner = center + numpy.dot(self.matrix,flags)
                if sum(abs(flags)) == 3:
                    # the solution is already known
                    pos = corner
                else:
                    # The basis vectors for the null space of the active
                    # constraints are the cell vectors corresponding to the
                    # inactive constraints.
                    N = self.matrix[:,flags==0]
                    Nt = N.transpose()
                    # The particular solution is the corner.
                    # The solution with the smallest norm:
                    H = numpy.dot(numpy.dot(N, numpy.linalg.inv(numpy.dot(Nt, N))), Nt)
                    pos = corner - numpy.dot(H, corner)
                # test if pos is inside the box AND the distance is below radius
                if abs(self.to_fractional(pos - center)).max() < 1.000001 and \
                   numpy.linalg.norm(pos) <= radius:
                    return True
            # found nothing
            return False

        result = []
        for index in iter_indexes():
            # Run over all potential indexes and check if they are in the
            # interaction sphere.
            if close_enough(self.to_cartesian(index)):
                result.append(index)
                if index[2] != 0:
                    result.append(-index)
        return numpy.array(result)


