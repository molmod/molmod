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
"""Tool for comparing molecules and molecular geometries

In order to measure the similarity between two molecules or graphs, one must
first create a similarity descriptor of both molecules:

a = SimilarityDescriptor.from_molecule(mol_a)
b = SimilarityDescriptor.from_molecule(mol_b)
print compute_similarity(a, b)

There are several ways to construct a similarity_descriptor:

a1 = SimilarityDescriptor.from_molecule(molecule)
a1 = SimilarityDescriptor.from_molecular_graph(molecular_graph)
a1 = SimilarityDescriptor.from_coordinates(coordinates, labels)
a1 = SimilarityDescriptor(distance_matrix, labels)
"""


from __future__ import print_function

import numpy as np

from molmod.ext import similarity_table_labels, similarity_table_distances, \
    similarity_measure


__all__ = ["SimilarityDescriptor", "compute_similarity"]


class SimilarityDescriptor(object):
    """A descriptor used to compute the similarity between two molecules"""
    def __init__(self, distance_matrix, labels):
        """Initialize a similarity descriptor

           Arguments:
             distance_matrix  --  a matrix with interatomic distances, this can
                                  also be distances in a graph
             labels  --  a list with integer labels used to identify atoms of
                         the same type
        """
        self.table_distances = similarity_table_distances(distance_matrix.astype(float))
        self.table_labels = similarity_table_labels(labels.astype(int))
        print(len(labels), len(distance_matrix))
        order = np.lexsort([self.table_labels[:, 1], self.table_labels[:, 0]])
        self.table_labels = self.table_labels[order]
        self.table_distances = self.table_distances[order]

    @classmethod
    def from_molecule(cls, molecule, labels=None):
        """Initialize a similarity descriptor

           Arguments:
             molecule  --  a Molecules object
             labels  --  a list with integer labels used to identify atoms of
                         the same type. When not given, the atom numbers from
                         the molecule are used.
        """
        if labels is None:
            labels = molecule.numbers
        return cls(molecule.distance_matrix, labels)

    @classmethod
    def from_molecular_graph(cls, molecular_graph, labels=None):
        """Initialize a similarity descriptor

           Arguments:
             molecular_graphs  --  A MolecularGraphs object
             labels  --  a list with integer labels used to identify atoms of
                         the same type. When not given, the atom numbers from
                         the molecular graph are used.
        """
        if labels is None:
            labels = molecular_graph.numbers.astype(int)
        return cls(molecular_graph.distances, labels)

    @classmethod
    def from_coordinates(cls, coordinates, labels):
        """Initialize a similarity descriptor

           Arguments:
             coordinates  --  a Nx3 numpy array
             labels  --  a list with integer labels used to identify atoms of
                         the same type
        """
        from molmod.ext import molecules_distance_matrix
        distance_matrix = molecules_distance_matrix(coordinates)
        return cls(distance_matrix, labels)


def compute_similarity(a, b, margin=1.0, cutoff=10.0):
    """Compute the similarity between two molecules based on their descriptors

       Arguments:
         a  --  the similarity measure of the first molecule
         b  --  the similarity measure of the second molecule
         margin  --  the sensitivity when comparing distances (default = 1.0)
         cutoff  --  don't compare distances longer than the cutoff (default = 10.0 au)

       When comparing two distances (always between two atom pairs with
       identical labels), the folowing formula is used:

       dav = (distance1+distance2)/2
       delta = abs(distance1-distance2)

       When the delta is within the margin and dav is below the cutoff:

         (1-dav/cutoff)*(cos(delta/margin/np.pi)+1)/2

       and zero otherwise. The returned value is the sum of such terms over all
       distance pairs with matching atom types. When comparing similarities it
       might be useful to normalize them in some way, e.g.

       similarity(a, b)/(similarity(a, a)*similarity(b, b))**0.5
    """
    return similarity_measure(
        a.table_labels, a.table_distances,
        b.table_labels, b.table_distances,
        margin, cutoff
    )
