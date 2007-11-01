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


from molmod.helpers import fsim
from molmod.binning import PositionedObject, SparseBinnedObjects, IntraAnalyseNeighboringObjects
from molmod.molecular_graphs import MolecularGraph

import numpy

__all__ = ["all_distances", "cutoff_distances", "calculate_similarity"]


def all_distances(molecule):
    all_distances = {}
    for i1, n1 in enumerate(molecule.numbers):
        for i2, n2 in enumerate(molecule.numbers[:i1]):
            if n1 < n2:
                tag = (n1,n2)
            else:
                tag = (n2,n1)
            l = all_distances.setdefault(tag, [])
            l.append(numpy.linalg.norm(molecule.coordinates[i1] - molecule.coordinates[i2]))
    for tag in all_distances.keys():
        all_distances[tag] = numpy.array(all_distances[tag])
    return all_distances


def cutoff_distances(molecule, radius, inter=False, graph=None):
    sbo = SparseBinnedObjects([
        PositionedObject(index, coordinate)
        for index, coordinate
        in enumerate(molecule.coordinates)
    ], radius*1.01)

    if inter and graph is None:
        graph = MolecularGraph(molecule)

    def delta_distance(positioned1, positioned2):
        delta = positioned2.coordinate - positioned1.coordinate
        distance = numpy.linalg.norm(delta)
        if distance <= radius:
            if inter and graph.get_distance(positioned1.id, positioned2.id) > 0:
                return
            return distance

    cutoff_distances = {}
    for (positioned1, positioned2), distance in IntraAnalyseNeighboringObjects(sbo, delta_distance)():
        n1 = molecule.numbers[positioned1.id]
        n2 = molecule.numbers[positioned2.id]
        if n1 < n2:
            tag = (n1,n2)
        else:
            tag = (n2,n1)
        l = cutoff_distances.setdefault(tag, [])
        l.append(distance)
    for tag in cutoff_distances.keys():
        cutoff_distances[tag] = numpy.array(cutoff_distances[tag])
    return cutoff_distances


def calculate_similarity(distances1, distances2, margin, radius=None):
    if radius is None:
        reciproke = 0
    else:
        reciproke = 1/radius
    similarity = 0.0
    for pair, row1 in distances1.iteritems():
        row2 = distances2.get(pair)
        if row2 is not None:
            similarity += fsim.similarity_rows(row1, row2, margin, reciproke)
    return similarity


