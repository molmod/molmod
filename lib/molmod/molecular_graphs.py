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


from molmod.molecules import Molecule
from molmod.graphs import Graph, GraphError, SubgraphMatchDefinition, ExactMatchDefinition, Match, OneToOne, MatchGenerator, CriteriaSet
from molmod.binning import IntraAnalyseNeighboringObjects, PositionedObject, SparseBinnedObjects
from molmod.data.bonds import bonds

import numpy, copy


__all__ = [
    "MolecularGraph", "generate_molecular_graph",
    "Anything", "MolecularCriterion", "MolecularOr", "MolecularAnd",
    "HasAtomNumber", "HasNumNeighbors", "HasNeighborNumbers", "BondLongerThan",
    "atom_criteria",
    "MolecularMixinMatchDefinition", "MolecularSubgraphMatchDefinition",
    "MolecularExactMatchDefinition", "BondMatchDefinition",
    "BendingAngleMatchDefinition", "DihedralAngleMatchDefinition",
    "OutOfPlaneMatchDefinition", "TetraMatchDefinition",
    "FullMatchError", "full_match",
]


class MolecularGraph(Graph):
    def __init__(self, pairs, numbers, ordered_nodes=None):
        Graph.__init__(self, pairs, ordered_nodes)
        self.numbers = numbers

    def __mul__(self, other):
        result = Graph.__mul__(self, other)
        result.__class__ = MolecularGraph
        numbers = numpy.zeros((other, len(self.numbers)), int)
        numbers[:] = self.numbers
        result.numbers = numbers.ravel()
        return result

    __rmul__ = __mul__

    def subgraph(self, indexes=None):
        result = Graph.subgraph(self, indexes)
        result.__class__ = MolecularGraph
        result.numbers = self.numbers[indexes]
        return result


def generate_molecular_graph(molecule, labels=None, unit_cell=None):
    if labels is None:
        labels = range(len(molecule.numbers))

    def yield_positioned_atoms():
        for index in xrange(len(labels)):
            yield PositionedObject(index, molecule.coordinates[index])

    binned_atoms = SparseBinnedObjects(yield_positioned_atoms(), bonds.max_length*bonds.bond_tolerance)

    def compare_function(positioned1, positioned2):
        delta = positioned2.coordinate - positioned1.coordinate
        if unit_cell is not None:
            delta = unit_cell.shortest_vector(delta)
        distance = numpy.linalg.norm(delta)
        if distance < binned_atoms.gridsize:
            bond_order = bonds.bonded(molecule.numbers[positioned1.id], molecule.numbers[positioned2.id], distance)
            if bond_order != None:
                return bond_order, distance

    bond_data = list(
        (frozenset([labels[positioned.id] for positioned in key]), data)
        for key, data
        in IntraAnalyseNeighboringObjects(binned_atoms, compare_function)(unit_cell)
    )
    pairs = set(key for key, data in bond_data)

    result = MolecularGraph(pairs, molecule.numbers, labels)
    result.bond_data = bond_data
    result.bond_orders = dict([(key, data[0]) for key, data in bond_data])
    result.bond_lengths = dict([(key, data[1]) for key, data in bond_data])
    return result


# molecular criteria


class MolecularCriterion(object):
    def init_graph(self, graph):
        self.graph = graph
        graph.init_index()

    def __call__(self, id):
        raise NotImplementedError


class Anything(MolecularCriterion):
    def __call__(self, id):
        return True


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
        return self.graph.numbers[self.graph.index[atom]] == self.number


class HasNumNeighbors(MolecularCriterion):
    def __init__(self, number):
        self.number = number

    def init_graph(self, graph):
        graph.init_neighbors()
        MolecularCriterion.init_graph(self, graph)

    def __call__(self, atom):
        return len(self.graph.neighbors[self.graph.index[atom]]) == self.number


class HasNeighborNumbers(MolecularCriterion):
    def __init__(self, numbers):
        self.numbers = list(numbers)
        self.numbers.sort()

    def init_graph(self, graph):
        graph.init_neighbors()
        MolecularCriterion.init_graph(self, graph)

    def __call__(self, atom):
        neighbors = self.graph.neighbors[atom]
        if not len(neighbors) == len(self.numbers):
            return
        neighbors = [self.graph.numbers[self.graph.index[neighbor]] for neighbor in neighbors]
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


class MolecularMixinMatchDefinition(object):
    def init_graph(self, graph):
        assert isinstance(graph, MolecularGraph)
        if self.criteria_sets is None:
            return
        for criteria_set in self.criteria_sets:
            for c in criteria_set.thing_criteria.itervalues():
                c.init_graph(graph)
            for c in criteria_set.relation_criteria.itervalues():
                c.init_graph(graph)
            for c in criteria_set.global_criteria:
                c.init_graph(graph)


class MolecularSubgraphMatchDefinition(SubgraphMatchDefinition, MolecularMixinMatchDefinition):
    def init_graph(self, graph, one_match):
        MolecularMixinMatchDefinition.init_graph(self, graph)
        SubgraphMatchDefinition.init_graph(self, graph, one_match)


class MolecularExactMatchDefinition(ExactMatchDefinition, MolecularMixinMatchDefinition):
    def init_graph(self, graph, one_match):
        MolecularMixinMatchDefinition.init_graph(self, graph)
        ExactMatchDefinition.init_graph(self, graph, one_match)


class BondMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([frozenset([0, 1])], [0, 1])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class BendingAngleMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
        ], [0, 1, 2])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class DihedralAngleMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
            frozenset([2, 3]),
        ], [0, 1, 2, 3])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class OutOfPlaneMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
        ], [0, 1, 2, 3])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class TetraMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
            frozenset([0, 4]),
        ], [0, 1, 2, 3, 4])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)



class FullMatchError(Exception):
    pass


def full_match(graph1, graph2):
    # given the graphs of two geometries of the same set of molecules, return
    # the global match between the two numberings
    mgs1 = [graph1.subgraph(group) for group in graph1.get_indexes_per_independent_graph()]
    mgs2 = [graph2.subgraph(group) for group in graph2.get_indexes_per_independent_graph()]

    if len(mgs1) != len(mgs2):
        #print "not same number of molecules"
        return

    matches = []

    while len(mgs1) > 0:
        subgraph1 = mgs1.pop()
        atom_criteria = dict((index, HasAtomNumber(number)) for index, number in zip(subgraph1.nodes, subgraph1.numbers))
        md = MolecularExactMatchDefinition(subgraph1, [CriteriaSet(atom_criteria)])
        matched = False
        for subgraph2 in mgs2:
            if len(subgraph1.nodes) != len(subgraph2.nodes):
                #print "size differ", len(subgraph1.nodes), len(subgraph2.nodes)
                continue
            if len(subgraph1.pairs) != len(subgraph2.pairs):
                #print "bonds differ", len(subgraph1.pairs), len(subgraph2.pairs)
                continue
            try:
                match = MatchGenerator(md,debug=False)(subgraph2,one_match=True).next()
                matches.append(match)
                mgs2.remove(subgraph2)
                matched = True
                #print "match"
                break
            except StopIteration:
                pass
            #print "not match"
        if not matched:
            return


    result = OneToOne()
    for match in matches:
        result.add_relations(match.forward.items())
    return result





