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


from molmod.ext import similarity_table_labels, similarity_table_distances, similarity_measure
from molmod.molecular_graphs import MolecularGraph
from molmod.molecules import Molecule

import numpy

__all__ = ["DistanceDescriptor", "distances_cor", "distances_dm", "compute_similarity"]


class DistanceDescriptor(object):
    def __init__(self, mol_or_graph, labels=None):
        if labels is None:
            self.labels = mol_or_graph.numbers.astype(numpy.int32)
        else:
            self.labels = labels.astype(numpy.int32)
        self.table_labels = similarity_table_labels(self.labels)
        if isinstance(mol_or_graph, Molecule):
            mol_or_graph.init_distance_matrix()
            self.table_distances = similarity_table_distances(mol_or_graph.distance_matrix)
        elif isinstance(mol_or_graph, MolecularGraph):
            self.table_distances = similarity_table_distances(mol_or_graph.distances)
        #order = self.table_labels.argsort(axis=0,kind='heapsort')
        order = numpy.lexsort([self.table_labels[:,1], self.table_labels[:,0]])
        self.table_labels = self.table_labels[order]
        self.table_distances = self.table_distances[order]

    def similarity(self, other, margin=1.0, cutoff=10.0):
        return similarity_measure(
            self.table_labels, self.table_distances,
            other.table_labels, other.table_distances,
            margin, cutoff
        )

    def norm(self, margin=1.0, cutoff=10.0):
        return numpy.sqrt(self.similarity(self, margin, cutoff))


def distances_cor(coordinates, labels):
    """Computes all interatomic distances, puts them into a table and sorts the table."""
    N = len(coordinates)
    if len(labels) != N:
        raise SimilarityError("The number of labels must match the size of the molecule.")
    all_distances = numpy.zeros((N*(N-1))/2, [("label1",int),("label2",int),("distance",float)])
    counter = 0
    for i1, l1 in enumerate(labels):
        for i2, l2 in enumerate(labels[:i1]):
            d = numpy.linalg.norm(coordinates[i1] - coordinates[i2])
            if l1 < l2:
                all_distances[counter] = (l1,l2,d)
            else:
                all_distances[counter] = (l2,l1,d)
            counter += 1
    all_distances.sort()
    return all_distances


def distances_dm(distance_matrix, labels):
    """Loads all interatomic distances, puts them into a table and sorts the table."""
    N = len(distance_matrix)
    if len(labels) != N:
        raise SimilarityError("The number of labels must match the size of the molecule.")
    all_distances = numpy.zeros((N*(N-1))/2, [("label1",int),("label2",int),("distance",float)])
    counter = 0
    for i1, l1 in enumerate(labels):
        for i2, l2 in enumerate(labels[:i1]):
            d = distance_matrix[i1,i2]
            if l1 < l2:
                all_distances[counter] = (l1,l2,d)
            else:
                all_distances[counter] = (l2,l1,d)
            counter += 1
    all_distances.sort()
    return all_distances


def compute_similarity(table_dist1, table_dist2, margin=1.0, cutoff=10.0):
    similarity = 0.0

    #print "table1"
    #print table_dist1
    #print "table2"
    #print table_dist2

    #print "="*20
    #print "="*20
    start2 = 0
    la2,lb2,d2 = table_dist2[start2]
    for la1,lb1,d1 in table_dist1:
        #print "TRY", la1,lb1, "  ", la2,lb2, "  ", cmp((la1,lb1),(la2,lb2))
        if (la1,lb1) < (la2,lb2):
            continue
        while (la1,lb1) > (la2,lb2):
            #print "HERE", start2
            start2 += 1
            if start2 == len(table_dist2):
                start2 = None
                break
            la2,lb2,d2 = table_dist2[start2]
        if (la1,lb1) < (la2,lb2):
            continue
        #print start2
        if start2 is None:
            break
        #print la1,lb1, "  ", la2,lb2, "  ", cmp((la1,lb1),(la2,lb2))
        current2 = start2
        while (la1,lb1) == (la2,lb2):
            dav = 0.5*(d1+d2)
            if dav < cutoff:
                delta = abs(d1-d2)
                if abs(delta) < margin:
                    scale = 1-dav/cutoff
                    similarity += scale*(numpy.cos(delta/margin/numpy.pi)+1)/2
            current2 += 1
            if current2 == len(table_dist2):
                break
            la2,lb2,d2 = table_dist2[current2]
        la2,lb2,d2 = table_dist2[start2]
    #print "="*20
    #print "="*20
    return similarity







