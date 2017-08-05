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
"""Representation, analysis and manipulation of molecular systems."""


from __future__ import division

from builtins import range
import numpy as np

from molmod.periodic import periodic
from molmod.units import angstrom
from molmod.utils import cached, ReadOnly, ReadOnlyAttribute
from molmod.molecular_graphs import MolecularGraph
from molmod.unit_cells import UnitCell
from molmod.transformations import fit_rmsd
from molmod.symmetry import compute_rotsym


__all__ = ["Molecule"]


class Molecule(ReadOnly):
    """Extensible class for molecular systems.

       Most attributes of the molecule object are treated as constants. If you
       want to modify the molecular geometry, just create a modified molecule
       object with the method :meth:`molmod.utils.ReadOnly.copy_with`. This facilitates the caching of
       derived quantities such as the distance matrix, while it imposes a
       cleaner coding style without a signifacant computational overhead.
    """
    def _check_coordinates(self, coordinates):
        """the number of rows must be the same as the length of the array numbers"""
        if len(coordinates) != self.size:
            raise TypeError("The number of coordinates does not match the "
                "length of the atomic numbers array.")

    def _check_masses(self, masses):
        """the size must be the same as the length of the array numbers"""
        if len(masses) != self.size:
            raise TypeError("The number of masses does not match the length of "
                "the atomic numbers array.")

    def _check_graph(self, graph):
        """the atomic numbers must match"""
        if graph.num_vertices != self.size:
            raise TypeError("The number of vertices in the graph does not "
                "match the length of the atomic numbers array.")
        # In practice these are typically the same arrays using the same piece
        # of memory. Just checking to be sure.
        if (self.numbers != graph.numbers).any():
            raise TypeError("The atomic numbers in the graph do not match the "
                "atomic numbers in the molecule.")

    def _check_symbols(self, symbols):
        """the size must be the same as the length of the array numbers and all elements must be strings"""
        if len(symbols) != self.size:
            raise TypeError("The number of symbols in the graph does not "
                "match the length of the atomic numbers array.")
        for symbol in symbols:
            if not isinstance(symbol, str):
                raise TypeError("All symbols must be strings.")

    numbers = ReadOnlyAttribute(np.ndarray, none=False, npdim=1, npdtype=int,
        doc="the atomic numbers")
    coordinates = ReadOnlyAttribute(np.ndarray, npdim=2, npshape=(None,3),
        npdtype=float, check=_check_coordinates, doc="atomic Cartesian "
        "coordinates")
    title = ReadOnlyAttribute(str, doc="a short description of the system")
    masses = ReadOnlyAttribute(np.ndarray, npdim=1, npdtype=float,
        check=_check_masses, doc="the atomic masses")
    graph = ReadOnlyAttribute(MolecularGraph, check=_check_graph,
        doc="the molecular graph with the atom connectivity")
    symbols = ReadOnlyAttribute(tuple, _check_symbols, doc="symbols for the "
        "atoms, which can be element names for force-field atom types")
    unit_cell = ReadOnlyAttribute(UnitCell, doc="description of the periodic "
        "boundary conditions")

    def __init__(self, numbers, coordinates=None, title=None, masses=None, graph=None, symbols=None, unit_cell=None):
        """
           Mandatory arguments:
            | ``numbers``  --  numpy array (1D, N elements) with the atomic numbers

           Optional keyword arguments:
            | ``coordinates``  --  numpy array (2D, Nx3 elements) Cartesian coordinates
            | ``title``  --  a string with the name of the molecule
            | ``massess``  --  a numpy array with atomic masses in atomic units
            | ``graph``  --  a MolecularGraph instance
            | ``symbols``  --  atomic elements or force-field atom-types
            | ``unit_cell``  --  the unit cell in case the system is periodic
        """
        self.numbers = numbers
        self.coordinates = coordinates
        self.title = title
        self.masses = masses
        self.graph = graph
        self.symbols = symbols
        self.unit_cell = unit_cell

    @classmethod
    def from_file(cls, filename):
        """Construct a molecule object read from the given file.

           The file format is inferred from the extensions. Currently supported
           formats are: ``*.cml``, ``*.fchk``, ``*.pdb``, ``*.sdf``, ``*.xyz``

           If a file contains more than one molecule, only the first one is
           read.

           Argument:
            | ``filename``  --  the name of the file containing the molecule

           Example usage::

             >>> mol = Molecule.from_file("foo.xyz")
        """
        # TODO: many different API's to load files. brrr...
        if filename.endswith(".cml"):
            from molmod.io import load_cml
            return load_cml(filename)[0]
        elif filename.endswith(".fchk"):
            from molmod.io import FCHKFile
            fchk = FCHKFile(filename, field_labels=[])
            return fchk.molecule
        elif filename.endswith(".pdb"):
            from molmod.io import load_pdb
            return load_pdb(filename)
        elif filename.endswith(".sdf"):
            from molmod.io import SDFReader
            return next(SDFReader(filename))
        elif filename.endswith(".xyz"):
            from molmod.io import XYZReader
            xyz_reader = XYZReader(filename)
            title, coordinates = next(xyz_reader)
            return Molecule(xyz_reader.numbers, coordinates, title, symbols=xyz_reader.symbols)
        else:
            raise ValueError("Could not determine file format for %s." % filename)

    size = property(lambda self: self.numbers.shape[0],
        doc="*Read-only attribute:* the number of atoms.")

    @cached
    def distance_matrix(self):
        """the matrix with all atom pair distances"""
        from molmod.ext import molecules_distance_matrix
        return molecules_distance_matrix(self.coordinates)

    @cached
    def mass(self):
        """the total mass of the molecule"""
        return self.masses.sum()

    @cached
    def com(self):
        """the center of mass of the molecule"""
        return (self.coordinates*self.masses.reshape((-1,1))).sum(axis=0)/self.mass

    @cached
    def inertia_tensor(self):
        """the intertia tensor of the molecule"""
        result = np.zeros((3,3), float)
        for i in range(self.size):
            r = self.coordinates[i] - self.com
            # the diagonal term
            result.ravel()[::4] += self.masses[i]*(r**2).sum()
            # the outer product term
            result -= self.masses[i]*np.outer(r,r)
        return result

    @cached
    def chemical_formula(self):
        """the chemical formula of the molecule"""
        counts = {}
        for number in self.numbers:
            counts[number] = counts.get(number, 0)+1
        items = []
        for number, count in sorted(counts.items(), reverse=True):
            if count == 1:
                items.append(periodic[number].symbol)
            else:
                items.append("%s%i" % (periodic[number].symbol, count))
        return "".join(items)

    def set_default_masses(self):
        """Set self.masses based on self.numbers and periodic table."""
        self.masses = np.array([periodic[n].mass for n in self.numbers])

    def set_default_graph(self):
        """Set self.graph to the default graph.

           This method is equivalent to::

              mol.graph = MolecularGraph.from_geometry(mol)

           with the default options, and only works if the graph object is not
           present yet.
           See :meth:`molmod.molecular_graphs.MolecularGraph.from_geometry`
           for more fine-grained control over the assignment of bonds.
        """
        self.graph = MolecularGraph.from_geometry(self)

    def set_default_symbols(self):
        """Set self.symbols based on self.numbers and the periodic table."""
        self.symbols = tuple(periodic[n].symbol for n in self.numbers)

    def write_to_file(self, filename):
        """Write the molecular geometry to a file.

           The file format is inferred from the extensions. Currently supported
           formats are: ``*.xyz``, ``*.cml``

           Argument:
            | ``filename``  --  a filename
        """
        # TODO: give all file format writers the same API
        if filename.endswith('.cml'):
            from molmod.io import dump_cml
            dump_cml(filename, [self])
        elif filename.endswith('.xyz'):
            from molmod.io import XYZWriter
            symbols = []
            for n in self.numbers:
                atom = periodic[n]
                if atom is None:
                    symbols.append("X")
                else:
                    symbols.append(atom.symbol)
            xyz_writer = XYZWriter(filename, symbols)
            xyz_writer.dump(self.title, self.coordinates)
            del xyz_writer
        else:
            raise ValueError("Could not determine file format for %s." % filename)

    def rmsd(self, other):
        """Compute the RMSD between two molecules.

           Arguments:
            | ``other``  --  Another molecule with the same atom numbers

           Return values:
            | ``transformation``  --  the transformation that brings 'self' into
                                  overlap with 'other'
            | ``other_trans``  --  the transformed coordinates of geometry 'other'
            | ``rmsd``  --  the rmsd of the distances between corresponding atoms in
                            'self' and 'other'

           Make sure the atoms in `self` and `other` are in the same order.

           Usage::

             >>> print molecule1.rmsd(molecule2)[2]/angstrom
        """
        if self.numbers.shape != other.numbers.shape or \
           (self.numbers != other.numbers).all():
            raise ValueError("The other molecule does not have the same numbers as this molecule.")
        return fit_rmsd(self.coordinates, other.coordinates)

    def compute_rotsym(self, threshold=1e-3*angstrom):
        """Compute the rotational symmetry number.

           Optional argument:
            | ``threshold``  --  only when a rotation results in an rmsd below the given
                                 threshold, the rotation is considered to transform the
                                 molecule onto itself.
        """
        # Generate a graph with a more permissive threshold for bond lengths:
        # (is convenient in case of transition state geometries)
        graph = MolecularGraph.from_geometry(self, scaling=1.5)
        try:
            return compute_rotsym(self, graph, threshold)
        except ValueError:
            raise ValueError("The rotational symmetry number can only be computed when the graph is fully connected.")
