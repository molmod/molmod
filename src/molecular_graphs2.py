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
        Graph.__init__(self, pairs)

    def yield_subgraphs(self, criteria_sets):
        subgraph = SymmetricGraph(criteria_sets.subpairs, criteria_sets.initiator)
        for tag, atom_criteria, bond_criteria, filter_tags in criteria_sets.yield_criteria():
            for atom_criterium in atom_criteria.itervalues():
                atom_criterium.set_molecular_graph(self)
            for bond_criterium in bond_criteria.itervalues():
                bond_criterium.set_molecular_graph(self)
            graph_filter = MatchFilterParameterized(
                subgraph,
                criteria_sets.calculation_tags,
                atom_criteria,
                bond_criteria,
                filter_tags
            )

            for match in subgraph.yield_matching_subgraphs(self):
                for parsed in graph_filter.parse(match):
                    yield tag, parsed


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
    def __init__(self, subgraph, atom_criteria, bond_criteria, node_tags={}):
        self.atom_criteria = atom_criteria
        self.bond_criteria = bond_criteria
        SubgraphMatchDefinition.__init__(self, subgraph, node_tags)

    def init_graph(self, graph):
        assert isinstance(graph, MolecularGraph)
        for c in self.atom_criteria.itervalues():
            c.init_graph(graph)
        for c in self.bond_criteria.itervalues():
            c.init_graph(graph)
        SubgraphMatchDefinition.init_graph(self, graph)

    def test_final_match(self, final_match):
        for node0, c in self.atom_criteria.iteritems():
            node1 = final_match.forward[node0]
            if not c(node1): return False
        for (node0a, node0b), c in self.bond_criteria.iteritems():
            node1a = final_match.forward[node0a]
            node1b = final_match.forward[node0b]
            if not c(frozenset([node1a, node1b])): return False
        return True


class BondMatchDefinition(MolecularMatchDefinition):
    def __init__(self, atom_criteria={}, bond_criteria={}, node_tags={}):
        subgraph = Graph([frozenset([0, 1])], [0, 1])
        MolecularMatchDefinition.__init__(self, subgraph, atom_criteria, bond_criteria, node_tags)


class BendingAngleMatchDefinition(MolecularMatchDefinition):
    def __init__(self, atom_criteria={}, bond_criteria={}, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
        ], [0, 1, 2])
        MolecularMatchDefinition.__init__(self, subgraph, atom_criteria, bond_criteria, node_tags)


class DihedralAngleMatchDefinition(MolecularMatchDefinition):
    def __init__(self, atom_criteria={}, bond_criteria={}, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
            frozenset([2, 3]),
        ], [0, 1, 2, 3])
        MolecularMatchDefinition.__init__(self, subgraph, atom_criteria, bond_criteria, node_tags)


class OutOfPlaneMatchDefinition(MolecularMatchDefinition):
    def __init__(self, atom_criteria={}, bond_criteria={}, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
        ], [0, 1, 2, 3])
        MolecularMatchDefinition.__init__(self, subgraph, atom_criteria, bond_criteria, node_tags)


class TetraMatchDefinition(MolecularMatchDefinition):
    def __init__(self, atom_criteria={}, bond_criteria={}, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
            frozenset([0, 4]),
        ], [0, 1, 2, 3, 4])
        MolecularMatchDefinition.__init__(self, subgraph, atom_criteria, bond_criteria, node_tags)
