#!/usr/bin/env python

from __future__ import print_function

import numpy

from molmod import *


class HarmonicEnergyTerm(object):
    """Base class for all energy terms."""
    def __init__(self, force_constant, indexes, icfn):
        """
           Arguments:
            | ``force_constant`` -- The force constant in atomic units.
            | ``indexes`` -- The indexes of the atoms in the internal
                             coordinate. The order must be the same as the order
                             of the mandatory arguments of icfn.
            | ``icfn`` -- a function from molmod.ic that can compute the
                          internal coordinate and its derivatives.
        """
        self.force_constant = force_constant
        self.indexes = indexes
        self.icfn = icfn

    def add_to_hessian(self, coordinates, hessian):
        """Add the contributions of this energy term to the Hessian

           Arguments:
            | ``coordinates`` -- A numpy array with 3N Cartesian coordinates.
            | ``hessian`` -- A matrix for the full Hessian to which this energy
                             term has to add its contribution.
        """
        # Compute the derivatives of the bond stretch towards the two cartesian
        # coordinates. The bond length is computed too, but not used.
        q, g = self.icfn(coordinates[list(self.indexes)], 1)
        # Add the contribution to the Hessian (an outer product)
        for ja, ia in enumerate(self.indexes):
            # ja is 0, 1, 2, ...
            # ia is i0, i1, i2, ...
            for jb, ib in enumerate(self.indexes):
                contrib = 2*self.force_constant*numpy.outer(g[ja], g[jb])
                hessian[3*ia:3*ia+3, 3*ib:3*ib+3] += contrib


class BondStretchTerm(HarmonicEnergyTerm):
    """A single bond-stretch energy term in the force-field term."""

    def __init__(self, force_constant, i0, i1):
        """
           Arguments:
            | ``force_constant`` -- The force constant in atomic units.
            | ``i0`` -- The atom index of the first atom in the bond.
            | ``i1`` -- The atom index of the second atom in the bond.
        """
        HarmonicEnergyTerm.__init__(self, force_constant, (i0, i1), bond_length)


class BendAngleTerm(HarmonicEnergyTerm):
    """A single bond-stretch energy term in the force-field term."""

    def __init__(self, force_constant, i0, i1, i2):
        """
           Arguments:
            | ``force_constant`` -- The force constant in atomic units.
            | ``i0`` -- The atom index of the first atom in the angle.
            | ``i1`` -- The atom index of the central atom in the angle.
            | ``i2`` -- The atom index of the third atom in the angle.
        """
        HarmonicEnergyTerm.__init__(self, force_constant, (i0, i1, i2), bend_angle)


class ForceField(object):
    """A container object for all force field terms."""
    def __init__(self, terms):
        """
           Argument:
            | ``terms`` -- a list of force-field terms
        """
        self.terms = terms
        # just print out the energy terms
        for term in self.terms:
            print("Energy term", term.icfn, term.force_constant, term.indexes)

    def hessian(self, coordinates):
        """Compute the force-field Hessian for the given coordinates.

           Argument:
            | ``coordinates`` -- A numpy array with the Cartesian atom
                                 coordinates, with shape (N,3).

           Returns:
            | ``hessian`` -- A numpy array with the Hessian, with shape (3*N,
                             3*N).
        """
        # N3 is 3 times the number of atoms.
        N3 = coordinates.size
        # Start with a zero hessian.
        hessian = numpy.zeros((N3,N3), float)
        # Add the contribution of each term.
        for term in self.terms:
            term.add_to_hessian(coordinates, hessian)
        return hessian


def setup_hydrocarbon_ff(graph):
    """Create a simple ForceField object for hydrocarbons based on the graph."""
    # A) Define parameters.
    # the bond parameters:
    bond_params = {
        (6, 1): 310*kcalmol/angstrom**2,
        (6, 6): 220*kcalmol/angstrom**2,
    }
    # for every (a, b), also add (b, a)
    for key, val in list(bond_params.items()):
        if key[0] != key[1]:
            bond_params[(key[1], key[0])] = val
    # the bend parameters
    bend_params = {
        (1, 6, 1): 35*kcalmol/rad**2,
        (1, 6, 6): 30*kcalmol/rad**2,
        (6, 6, 6): 60*kcalmol/rad**2,
    }
    # for every (a, b, c), also add (c, b, a)
    for key, val in list(bend_params.items()):
        if key[0] != key[2]:
            bend_params[(key[2], key[1], key[0])] = val

    # B) detect all internal coordinates and corresponding energy terms.
    terms = []
    # bonds
    for i0, i1 in graph.edges:
        K = bond_params[(graph.numbers[i0], graph.numbers[i1])]
        terms.append(BondStretchTerm(K, i0, i1))
    # bends (see b_bending_angles.py for the explanation)
    for i1 in range(graph.num_vertices):
        n = list(graph.neighbors[i1])
        for index, i0 in enumerate(n):
            for i2 in n[:index]:
                K = bend_params[(graph.numbers[i0], graph.numbers[i1], graph.numbers[i2])]
                terms.append(BendAngleTerm(K, i0, i1, i2))

    # C) Create and return the force field
    return ForceField(terms)


# This if block is only executed when this file is ran as a program, and not
# when it is loaded as a module.
if __name__ == "__main__":
    propane = Molecule.from_file("propane.xyz")
    propane.set_default_graph()
    ff = setup_hydrocarbon_ff(propane.graph)
    hessian = ff.hessian(propane.coordinates)
    print("The Hessian in kcal/mol/angstrom**2")
    unit = kcalmol/angstrom**2
    for row in hessian:
        print(" ".join("% 5.0f" % (v/unit) for v in row))
