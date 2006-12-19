# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from molmod.graphs2 import Graph, SubgraphMatchDefinition, Match
from molmod.binning import IntraAnalyseNeighboringObjects, PositionedObject, SparseBinnedObjects
from molmod.data import bonds

import math, numpy


class MolecularGraph(Graph):
    def __init__(self, molecule):
        self.molecule = molecule

        def yield_positioned_atoms():
            for index in xrange(len(self.molecule.numbers)):
                yield PositionedObject(index, self.molecule.coordinates[index])

        binned_atoms = SparseBinnedObjects(yield_positioned_atoms(), bonds.max_length*bonds.bond_tolerance)

        def compare_function(positioned1, positioned2):
            delta = positioned2.coordinate - positioned1.coordinate
            distance = math.sqrt(numpy.dot(delta, delta))
            if distance < binned_atoms.gridsize:
                bond_order = bonds.bonded(self.molecule.numbers[positioned1.id], self.molecule.numbers[positioned2.id], distance)
                if bond_order != None:
                    return bond_order, distance

        bond_data = list(
            (frozenset([positioned.id for positioned in key]), data)
            for key, data
            in IntraAnalyseNeighboringObjects(binned_atoms, compare_function)()
        )
        pairs = set(key for key, data in bond_data)
        self.bond_orders = dict([(key, data[0]) for key, data in bond_data])
        self.bond_lengths = dict([(key, data[1]) for key, data in bond_data])
        Graph.__init__(self, pairs, range(len(molecule.numbers)))


# molecular criteria


class Anything(object):
    def __call__(self, id):
        return True


class MolecularCriterion(object):
    def init_graph(self, graph):
        self.graph = graph


class MolecularOr(object):
    def __init__(self, *criteria):
        self.criteria = criteria

    def init_graph(self, graph):
        for c in self.criteria:
            c.init_graph(graph)

    def __call__(self, id):
        for c in self.criteria:
            if c(id):
                return True
        return False


class MolecularAnd(object):
    def __init__(self, *criteria):
        self.criteria = criteria

    def init_graph(self, graph):
        for c in self.criteria:
            c.init_graph(graph)

    def __call__(self, id):
        for c in self.criteria:
            if not c(id):
                return False
        return True


class HasAtomNumber(MolecularCriterion):
    def __init__(self, number):
        self.number = number

    def __call__(self, atom):
        return self.graph.molecule.numbers[atom] == self.number


class HasNumNeighbors(MolecularCriterion):
    def __init__(self, number):
        self.number = number

    def __call__(self, atom):
        return len(self.graph.neighbors[atom]) == self.number


class HasNeighborNumbers(MolecularCriterion):
    def __init__(self, numbers):
        self.numbers = list(numbers)
        self.numbers.sort()

    def __call__(self, atom):
        neighbors = self.graph.neighbors[atom]
        if not len(neighbors) == len(self.numbers):
            return
        neighbors = [self.graph.molecule.numbers[neighbor] for neighbor in neighbors]
        neighbors.sort()
        return neighbors == self.numbers


class BondLongerThan(MolecularCriterion):
    def __init__(self, length):
        self.length = length

    def __call__(self, pair):
        return self.graph.bond_lengths[pair] > self.length


def atom_criteria(*params):
    result = {}
    for index, param in enumerate(params):
        if param is None:
            continue
        elif isinstance(param, int):
            result[index] = HasAtomNumber(param)
        else:
            result[index] = param
    return result


# match definitions


class MolecularMatchDefinition(SubgraphMatchDefinition):
    def init_graph(self, graph):
        assert isinstance(graph, MolecularGraph)
        for criteria_set in self.criteria_sets:
            for c in criteria_set.thing_criteria.itervalues():
                c.init_graph(graph)
            for c in criteria_set.relation_criteria.itervalues():
                c.init_graph(graph)
            for c in criteria_set.global_criteria:
                c.init_graph(graph)
        SubgraphMatchDefinition.init_graph(self, graph)


class BondMatchDefinition(MolecularMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([frozenset([0, 1])], [0, 1])
        MolecularMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class BendingAngleMatchDefinition(MolecularMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
        ], [0, 1, 2])
        MolecularMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class DihedralAngleMatchDefinition(MolecularMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
            frozenset([2, 3]),
        ], [0, 1, 2, 3])
        MolecularMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class OutOfPlaneMatchDefinition(MolecularMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
        ], [0, 1, 2, 3])
        MolecularMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class TetraMatchDefinition(MolecularMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
            frozenset([0, 4]),
        ], [0, 1, 2, 3, 4])
        MolecularMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)
