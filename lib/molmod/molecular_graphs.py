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


from molmod.molecules import Molecule
from molmod.graphs import Graph, GraphError, SubgraphMatchDefinition, ExactMatchDefinition, Match, OneToOne, MatchGenerator, CriteriaSet
from molmod.binning import IntraAnalyseNeighboringObjects, PositionedObject, SparseBinnedObjects
from molmod.data import bonds, periodic
from molmod.transformations import rotation_around_center
from molmod.vectors import random_orthonormal

import numpy, copy


class MolecularGraph(Graph):
    def __init__(self, molecule, labels=None, unit_cell=None, pairs=None):
        self.molecule = molecule
        if labels is None:
            labels = range(len(molecule.numbers))

        if pairs is None:
            def yield_positioned_atoms():
                for index in xrange(len(labels)):
                    yield PositionedObject(index, self.molecule.coordinates[index])

            binned_atoms = SparseBinnedObjects(yield_positioned_atoms(), bonds.max_length*bonds.bond_tolerance)

            def compare_function(positioned1, positioned2):
                delta = positioned2.coordinate - positioned1.coordinate
                if unit_cell is not None:
                    delta = unit_cell.shortest_vector(delta)
                distance = numpy.linalg.norm(delta)
                if distance < binned_atoms.gridsize:
                    bond_order = bonds.bonded(self.molecule.numbers[positioned1.id], self.molecule.numbers[positioned2.id], distance)
                    if bond_order != None:
                        return bond_order, distance

            bond_data = list(
                (frozenset([labels[positioned.id] for positioned in key]), data)
                for key, data
                in IntraAnalyseNeighboringObjects(binned_atoms, compare_function)(unit_cell)
            )
            pairs = set(key for key, data in bond_data)
        else:
            pairs = set(frozenset(pair) for pair in pairs)
            bond_data = []
            for pair in pairs:
                a,b = pair
                delta = molecule.coordinates[a] - molecule.coordinates[b]
                if unit_cell is not None:
                    delta = unit_cell.shortest_vector(delta)
                distance = numpy.linalg.norm(delta)
                bond_order = bonds.bonded(self.molecule.numbers[a], self.molecule.numbers[b], distance)
                bond_data.append((pair, (bond_order, distance)))

        self.bond_orders = dict([(key, data[0]) for key, data in bond_data])
        self.bond_lengths = dict([(key, data[1]) for key, data in bond_data])
        Graph.__init__(self, pairs, labels)

    def subgraph(self, subnodes=None, subindices=None):
        subnodes, subindices = self._complete_args(subnodes, subindices)
        molecule = Molecule()
        molecule.numbers = self.molecule.numbers[subindices]
        molecule.coordinates = self.molecule.coordinates[subindices]
        return MolecularGraph(molecule, subnodes)

    def randomized_molecule(self, max_tries=1000, bond_fraction=0.2, dihedral_rotation=numpy.pi, bending_rotation=0.14, nonbond_threshold_factor=2.0):
        self.init_distances()
        radii = numpy.array([periodic[number].radius for number in self.molecule.numbers], float)

        def check(molecule):
            # check that no atoms overlap
            for index1, atom1 in enumerate(self.nodes):
                for index2, atom2 in enumerate(self.nodes[:index1]):
                    if self.get_distance(atom1, atom2) > 2:
                        distance = numpy.linalg.norm(molecule.coordinates[index1] - molecule.coordinates[index2])
                        if distance < nonbond_threshold_factor*(radii[index1] + radii[index2]):
                            return False
            return True

        for counter in xrange(max_tries):
            result = copy.deepcopy(self.molecule)
            for atom1, atom2 in self.pairs:
                delta = result.coordinates[atom1] - result.coordinates[atom2]
                try:
                    half_atoms = self.get_half(atom1, atom2)
                except GraphError:
                    continue
                if half_atoms is not None:
                    # Random bond stretch
                    translation = delta*bond_fraction*numpy.random.uniform(-1, 1)
                    for half_atom in half_atoms:
                        result.coordinates[half_atom] += translation
                    # Random dihedral rotation
                    direction = delta / numpy.linalg.norm(delta)
                    R1 = rotation_around_center(
                        result.coordinates[atom1],
                        numpy.random.uniform(-dihedral_rotation, dihedral_rotation),
                        direction,
                    )
                    for half_atom in half_atoms:
                        result.coordinates[half_atom] = R1.vector_apply(result.coordinates[half_atom])
                    # Random bending angle rotation 1
                    axis = random_orthonormal(direction)
                    R2 = rotation_around_center(
                        result.coordinates[atom1],
                        numpy.random.uniform(-bending_rotation, +bending_rotation),
                        axis,
                    )
                    for half_atom in half_atoms:
                        result.coordinates[half_atom] = R2.vector_apply(result.coordinates[half_atom])
                    # Random bending angle rotation 2
                    axis = random_orthonormal(direction)
                    R3 = rotation_around_center(
                        result.coordinates[atom2],
                        numpy.random.uniform(-bending_rotation, +bending_rotation),
                        axis,
                    )
                    for half_atom in half_atoms:
                        result.coordinates[half_atom] = R3.vector_apply(result.coordinates[half_atom])

            result.coordinates -= result.coordinates.mean(axis=0)

            if not check(result):
                continue

            return result

# molecular criteria


class Anything(object):
    def __call__(self, id):
        return True


class MolecularCriterion(object):
    def init_graph(self, graph):
        self.graph = graph
        graph.init_index()


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
        return self.graph.molecule.numbers[self.graph.index[atom]] == self.number


class HasNumNeighbors(MolecularCriterion):
    def __init__(self, number):
        self.number = number

    def __call__(self, atom):
        return len(self.graph.neighbors[self.graph.index[atom]]) == self.number


class HasNeighborNumbers(MolecularCriterion):
    def __init__(self, numbers):
        self.numbers = list(numbers)
        self.numbers.sort()

    def __call__(self, atom):
        neighbors = self.graph.neighbors[atom]
        if not len(neighbors) == len(self.numbers):
            return
        neighbors = [self.graph.molecule.numbers[self.graph.index[neighbor]] for neighbor in neighbors]
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
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([frozenset([0, 1])], [0, 1])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class BendingAngleMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
        ], [0, 1, 2])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class DihedralAngleMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([1, 2]),
            frozenset([2, 3]),
        ], [0, 1, 2, 3])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class OutOfPlaneMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
        subgraph = Graph([
            frozenset([0, 1]),
            frozenset([0, 2]),
            frozenset([0, 3]),
        ], [0, 1, 2, 3])
        MolecularSubgraphMatchDefinition.__init__(self, subgraph, criteria_sets, node_tags)


class TetraMatchDefinition(MolecularSubgraphMatchDefinition):
    def __init__(self, criteria_sets, node_tags={}):
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
    mgs1 = [graph1.subgraph(group) for group in graph1.get_nodes_per_independent_graph()]
    mgs2 = [graph2.subgraph(group) for group in graph2.get_nodes_per_independent_graph()]

    if len(mgs1) != len(mgs2):
        #print "not same number of molecules"
        return

    matches = []

    while len(mgs1) > 0:
        subgraph1 = mgs1.pop()
        atom_criteria = dict((index, HasAtomNumber(number)) for index, number in zip(subgraph1.nodes, subgraph1.molecule.numbers))
        md = MolecularExactMatchDefinition(subgraph1, [CriteriaSet("foo", atom_criteria)])
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


