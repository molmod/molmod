#!/usr/bin/env python

from __future__ import print_function

import numpy

from molmod import *
from molmod.io import FCHKFile


class InternalCoordinate(object):
    """Abstract base class for all internal coordinates."""
    def __init__(self, indexes, icfn, conversion=1.0):
        """
           Arguments:
            | ``indexes`` -- The indexes of the atoms in the internal
                             coordinate. The order must be the same as the order
                             of the mandatory arguments of icfn.
            | ``icfn`` -- a function from molmod.ic that can compute the
                          internal coordinate and its derivatives.
            | ``conversion`` -- In case the internal coordinate does not have a
                                unit of length, then this conversion factor is
                                used to convert it to a length unit. This way,
                                the Jacobian becomes a dimensionless constant.

           All the Jacobian-logic is implemented in this abstract class.
        """
        self.indexes = indexes
        self.icfn = icfn
        self.conversion = conversion

    def fill_jacobian_column(self, jaccol, coordinates):
        """Fill in a column of the Jacobian.

           Arguments:
            | ``jaccol`` -- The column of Jacobian to which the result must be
                            added.
            | ``coordinates`` -- A numpy array with Cartesian coordinates,
                                 shape=(N,3)
        """
        q, g = self.icfn(coordinates[list(self.indexes)], 1)
        for i, j in enumerate(self.indexes):
            jaccol[3*j:3*j+3] += g[i]
        return jaccol


class BondLength(InternalCoordinate):
    def __init__(self, i0, i1):
        InternalCoordinate.__init__(self, (i0, i1), bond_length)


class BendingAngle(InternalCoordinate):
    def __init__(self, i0, i1, i2):
        InternalCoordinate.__init__(self, (i0, i1, i2), bend_angle, angstrom/(5*deg))


class DihedralAngle(InternalCoordinate):
    def __init__(self, i0, i1, i2, i3):
        InternalCoordinate.__init__(self, (i0, i1, i2, i3), dihed_angle, angstrom/(5*deg))


def setup_ics(graph):
    """Make a list of internal coordinates based on the graph

       Argument:
        | ``graph`` -- A Graph instance.

       The list of internal coordinates will include all bond lengths, all
       bending angles, and all dihedral angles.
    """
    ics = []
    # A) Collect all bonds.
    for i0, i1 in graph.edges:
        ics.append(BondLength(i0, i1))
    # B) Collect all bends. (see b_bending_angles.py for the explanation)
    for i1 in range(graph.num_vertices):
        n = list(graph.neighbors[i1])
        for index, i0 in enumerate(n):
            for i2 in n[:index]:
                ics.append(BendingAngle(i0, i1, i2))
    # C) Collect all dihedrals.
    for i1, i2 in graph.edges:
        for i0 in graph.neighbors[i1]:
            if i0==i2:
                # All four indexes must be different.
                continue
            for i3 in graph.neighbors[i2]:
                if i3==i1 or i3==i0:
                    # All four indexes must be different.
                    continue
                ics.append(DihedralAngle(i0, i1, i2, i3))
    return ics


def compute_jacobian(ics, coordinates):
    """Construct a Jacobian for the given internal and Cartesian coordinates

       Arguments:
        | ``ics`` -- A list of internal coordinate objects.
        | ``coordinates`` -- A numpy array with Cartesian coordinates,
                             shape=(N,3)

       The return value will be a numpy array with the Jacobian matrix. There
       will be a column for each internal coordinate, and a row for each
       Cartesian coordinate (3*N rows).
    """
    N3 = coordinates.size
    jacobian = numpy.zeros((N3, len(ics)), float)
    for j, ic in enumerate(ics):
        # Let the ic object fill in each column of the Jacobian.
        ic.fill_jacobian_column(jacobian[:,j], coordinates)
    return jacobian


# This if block is only executed when this file is ran as a program, and not
# when it is loaded as a module.
if __name__ == "__main__":
    # Load the formatted checkpoint file with the frequency computation. This
    # file also contains the atomic numbers and the coordinates of the atoms,
    # and therefore one can access the dopamine molecule object through
    # fchk.molecule.
    fchk = FCHKFile("dopamine.fchk")
    # Set the default graph for the construction of the internal coordinates:
    fchk.molecule.set_default_graph()
    # Setup a list of internal coordinates
    ics = setup_ics(fchk.molecule.graph)
    # Compute the Jacobian.
    J = compute_jacobian(ics, fchk.molecule.coordinates)
    # Compute the pseudo-inverse, using a loose threshold for the singular
    # values to filter out equivalent internal coordinates.
    Jinv = numpy.linalg.pinv(J, 1e-5)
    # Get the Hessian in Cartesian coordinates.
    H = fchk.get_hessian()
    # Transform to internal coordinates.
    K = numpy.dot(Jinv, numpy.dot(H, Jinv.transpose()))
    # Make a nice printout of K.
    print("The Hessian in internal coordinates in kcal/mol/angstrom**2")
    unit = kcalmol/angstrom**2
    for row in K:
        print(" ".join("% 5.0f" % (v/unit) for v in row))
