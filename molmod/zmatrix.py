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
"""Conversion between Cartesian coordinates and ZMatrices"""


from builtins import range
import numpy as np

import molmod.ic as ic
from molmod.vectors import random_orthonormal


__all__ = ["ZMatrixGenerator", "zmat_to_cart"]


class ZMatrixGenerator(object):
    """ZMatrixGenerator converts cartesian coordinates to zmatrix notation.

    One detail is very different from conventional zmatrices: the references to
    other atoms are relative. This makes it easy to concatenate zmatrices.
    """

    dtype = [
        ("number", int),
        ("distance", float), ("rel1", int),
        ("angle", float), ("rel2", int),
        ("dihed", float), ("rel3", int),
        # refs are strictly positive integers. they refer to an atom that is
        # located 'ref' rows earlier in the zmatrix
    ]
    def __init__(self, graph):
        """Initialize a ZMatrix generator

           Argument:
             graph  --  an inseparable molecular graph
        """
        # Check that the graph is inseparable.
        if len(graph.independent_vertices) > 1:
            raise ValueError("The graph must be inseparable.")
        self.graph = graph
        # First step is to reorder the atoms so that the zmatrix is somewhat
        # workable. The constraint is that each subsequent atom must be bonded
        # to one of the previous atoms. It is the user's responsability to make
        # sure that heavier atoms come first.
        new_order = [0]
        # We will try to take the original order as long as it satisfies the
        # constraint.
        for i in range(1, graph.num_vertices):
            if (graph.distances[i][new_order] == 1).any():
                new_order.append(i)
            else:
                break
        # If not all atoms are listed in new_order, we continue adding the
        # remaining atoms in an order that satisfies the constraint. This will
        # break the original order.
        remaining = list(range(len(new_order), graph.num_vertices))
        while len(remaining) > 0:
            pivot = remaining.pop()
            if (graph.distances[pivot][new_order] == 1).any():
                new_order.append(pivot)
            else:
                remaining.insert(0, pivot)
        # store the orders as indexes
        self.old_index = np.array(new_order)
        self.new_index = np.zeros(self.old_index.shape, int)
        for i, j in enumerate(self.old_index):
            self.new_index[j] = i

    def _get_new_ref(self, existing_refs):
        """Get a new reference atom for a row in the ZMatrix

           The reference atoms should obey the following conditions:
             - They must be different
             - They must be neighbours in the bond graph
             - They must have an index lower than the current atom

           If multiple candidate refs can be found, take the heaviest atom
        """
        # ref0 is the atom whose position is defined by the current row in the
        # zmatrix.
        ref0 = existing_refs[0]
        for ref in existing_refs:
            # try to find a neighbor of the ref that can serve as the new ref
            result = None
            for n in sorted(self.graph.neighbors[ref]):
                if self.new_index[n] > self.new_index[ref0]:
                    # index is too high, zmatrix rows can't refer to future
                    # atoms
                    continue
                if n in existing_refs:
                    # ref is already in use
                    continue
                if result is None or self.graph.numbers[n] <= self.graph.numbers[result]:
                    # acceptable ref, prefer heaviest atom
                    result = n
            if result is not None:
                return result
        raise RuntimeError("Could not find new reference.")

    def cart_to_zmat(self, coordinates):
        """Convert cartesian coordinates to ZMatrix format

           Argument:
             coordinates  --  Cartesian coordinates (numpy array Nx3)

           The coordinates must match with the graph that was used to initialize
           the ZMatrixGenerator object.
        """
        N = len(self.graph.numbers)
        if coordinates.shape != (N, 3):
            raise ValueError("The shape of the coordinates must be (%i, 3)" % N)
        result = np.zeros(N, dtype=self.dtype)
        for i in range(N):
            ref0 = self.old_index[i]
            rel1 = -1
            rel2 = -1
            rel3 = -1
            distance = 0
            angle = 0
            dihed = 0
            if i > 0:
                ref1 = self._get_new_ref([ref0])
                distance = np.linalg.norm(coordinates[ref0]-coordinates[ref1])
                rel1 = i - self.new_index[ref1]
            if i > 1:
                ref2 = self._get_new_ref([ref0, ref1])
                angle, = ic.bend_angle(coordinates[[ref0, ref1, ref2]])
                rel2 = i - self.new_index[ref2]
            if i > 2:
                ref3 = self._get_new_ref([ref0, ref1, ref2])
                dihed, = ic.dihed_angle(coordinates[[ref0, ref1, ref2, ref3]])
                rel3 = i - self.new_index[ref3]
            result[i] = (self.graph.numbers[i], distance, rel1, angle, rel2, dihed, rel3)
        return result


def zmat_to_cart(zmat):
    """Converts a ZMatrix back to cartesian coordinates."""

    numbers = zmat["number"]
    N = len(numbers)
    coordinates = np.zeros((N, 3), float)

    # special cases for the first coordinates
    coordinates[1, 2] = zmat["distance"][1]
    if zmat["rel1"][2] == 1:
        sign = -1
    else:
        sign = 1
    coordinates[2, 2] = zmat["distance"][2]*sign*np.cos(zmat["angle"][2])
    coordinates[2, 1] = zmat["distance"][2]*sign*np.sin(zmat["angle"][2])
    coordinates[2] += coordinates[2-zmat["rel1"][2]]

    ref0 = 3
    for (number, distance, rel1, angle, rel2, dihed, rel3) in zmat[3:]:
        ref1 = ref0 - rel1
        ref2 = ref0 - rel2
        ref3 = ref0 - rel3
        if ref1 < 0: ref1 = 0
        if ref2 < 0: ref2 = 0
        if ref3 < 0: ref3 = 0
        # define frame axes
        origin = coordinates[ref1]
        new_z = coordinates[ref2] - origin
        norm_z = np.linalg.norm(new_z)
        if norm_z < 1e-15:
            new_z = np.array([0, 0, 1], float)
        else:
            new_z /= np.linalg.norm(new_z)
        new_x = coordinates[ref3] - origin
        new_x -= np.dot(new_x, new_z)*new_z
        norm_x = np.linalg.norm(new_x)
        if norm_x < 1e-15:
            new_x = random_orthonormal(new_z)
        else:
            new_x /= np.linalg.norm(new_x)
        # we must make our axes frame left handed due to the poor IUPAC
        # definition of the sign of a dihedral angle.
        new_y = -np.cross(new_z, new_x)

        # coordinates of new atom:
        x = distance*np.cos(dihed)*np.sin(angle)
        y = distance*np.sin(dihed)*np.sin(angle)
        z = distance*np.cos(angle)
        coordinates[ref0] = origin + x*new_x + y*new_y + z*new_z
        # loop
        ref0 += 1

    return numbers, coordinates
