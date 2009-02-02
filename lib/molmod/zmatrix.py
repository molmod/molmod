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


import molmod.ic as ic

import numpy


__all__ = ["ZMatrixError", "ZMatrixGenerator", "zmat_to_cart"]


class ZMatrixError(Exception):
    pass


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
        # Check that the graph is inseparable.
        if len(graph.independent_nodes) > 1:
            raise ZMatrixError("The graph must be inseparable.")
        self.graph = graph
        # First step is to reorder the atoms so that the zmatrix is somewhat
        # workable. The constraint is that each subsequent atom must be bonded
        # to one of the previous atoms. It is the user's responsability to make
        # sure that heavier atoms come first.
        new_order = [0]
        # We will try to take the original order as long as it satisfies the
        # constraint.
        for i in xrange(1, graph.num_nodes):
            if (graph.distances[i][new_order] == 1).any():
                new_order.append(i)
            else:
                break
        # If not all atoms are listed in new_order, we continue adding the
        # remaining atoms in an order that satisfies the constraint. This will
        # break the original order.
        remaining = range(len(new_order), graph.num_nodes)
        while len(remaining) > 0:
            pivot = remaining.pop()
            if (graph.distances[pivot][new_order] == 1).any():
                new_order.append(pivot)
            else:
                remaining.insert(0, pivot)
        # store the orders as indexes
        self.old_index = numpy.array(new_order)
        self.new_index = numpy.zeros(self.old_index.shape, int)
        for i,j in enumerate(self.old_index):
            self.new_index[j] = i

    def get_new_ref(self, existing_refs):
        ref0 = existing_refs[0]
        for ref in existing_refs:
            result = None
            for n in sorted(self.graph.neighbors[ref]):
                if self.new_index[n] > self.new_index[ref0]: continue
                if n in existing_refs: continue
                if result is None or self.graph.numbers[n] <= self.graph.numbers[result]:
                    result = n
            if result is not None:
                return result
        raise ZMatrixException("Could not find new reference.")

    def cart_to_zmat(self, numbers, coordinates):
        N = len(numbers)
        result = numpy.zeros(N, dtype=self.dtype)
        for i in xrange(N):
            ref0 = self.old_index[i]
            rel1 = -1
            rel2 = -1
            rel3 = -1
            distance = 0
            angle = 0
            dihed = 0
            if i > 0:
                ref1 = self.get_new_ref([ref0])
                distance = numpy.linalg.norm(coordinates[ref0]-coordinates[ref1])
                rel1 = i - self.new_index[ref1]
            if i > 1:
                ref2 = self.get_new_ref([ref0, ref1])
                angle, = ic.bend_angle(coordinates[ref0], coordinates[ref1], coordinates[ref2])
                rel2 = i - self.new_index[ref2]
            if i > 2:
                ref3 = self.get_new_ref([ref0, ref1, ref2])
                dihed, = ic.dihed_angle(coordinates[ref0], coordinates[ref1], coordinates[ref2], coordinates[ref3])
                rel3 = i - self.new_index[ref3]
            result[i] = (numbers[i], distance, rel1, angle, rel2, dihed, rel3)
        return result


def zmat_to_cart(zmat):
    """Converts a ZMatrix back to cartesian coordinates."""

    numbers = zmat["number"]
    N = len(numbers)
    coordinates = numpy.zeros((N,3), float)

    # special cases for the first coordinates
    coordinates[1,2] = zmat["distance"][1]
    if zmat["rel1"][2] == 1:
        sign = -1
    else:
        sign = 1
    coordinates[2,2] = zmat["distance"][2]*sign*numpy.cos(zmat["angle"][2])
    coordinates[2,1] = zmat["distance"][2]*sign*numpy.sin(zmat["angle"][2])
    coordinates[2] += coordinates[2-zmat["rel1"][2]]

    ref0 = 3
    for (number, distance, rel1, angle, rel2, dihed, rel3) in zmat[3:]:
        ref1 = ref0 - rel1
        ref2 = ref0 - rel2
        ref3 = ref0 - rel3
        # define frame axes
        origin = coordinates[ref1]
        new_z = coordinates[ref2] - origin
        new_z /= numpy.linalg.norm(new_z)
        new_x = coordinates[ref3] - origin
        new_x -= numpy.dot(new_x, new_z)*new_z
        new_x /= numpy.linalg.norm(new_x)
        new_y = numpy.cross(new_z, new_x)

        # coordinates of new atom:
        x = distance*numpy.cos(dihed)*numpy.sin(angle)
        y = distance*numpy.sin(dihed)*numpy.sin(angle)
        z = distance*numpy.cos(angle)
        coordinates[ref0] = origin + x*new_x + y*new_y + z*new_z
        # loop
        ref0 += 1

    return numbers, coordinates


