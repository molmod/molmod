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

    def get_canonical(self):
        """Returns a canonical representation of the molecular graph.

        This means that two graphs that only differ in the order of the nodes
        and their specific numbers, will result in the same canonical representation.

        The result is a tuple of two lists. The first list contains the atom
        numbers. The second list contains two-tuples that describe the bonds.
        """

        # To make this work, we need a unique starting point. If necessary,
        # symmetries will be used to make this more efficient.
        self.init_distances()
        self.init_neighbors()
        self.init_nodes()

        # pick the most central atom(s)
        longest_distances = self.distances.sum(axis=0)
        starting_atoms = (longest_distances == longest_distances.min()).nonzero()[0]
        #print "starting_atoms", starting_atoms

        if len(starting_atoms) > 1:
            # pick the heaviest atom(s)
            atom_numbers = self.numbers[starting_atoms]
            starting_atoms = starting_atoms[atom_numbers == atom_numbers.max()]
            #print "starting_atoms", starting_atoms

        #if len(starting_atoms) > 1:
        #    self.init_symmetries()
        #    # pick the atom(s) with the least number of symmetric equivalents
        #    num_equivalents = numpy.array([len(self.equivalent_nodes[i]) for i in starting_atoms])
        #    starting_atoms = starting_atoms[num_equivalents == num_equivalents.min()]
        #    print "starting_atoms", starting_atoms

            # if symmetric equivalents exist, pick from each group just one atom
        #    groups = set([])
        #    for i in starting_atoms:
        #        groups.add(frozenset(self.equivalent_nodes[i]))

        #    starting_atoms = numpy.array([min(group) for group in groups])
        #    print "starting_atoms", starting_atoms

        # It might very well be that the number of starting atoms is larger than one.
        # We will use each one to start a canonical representation. We will sorte
        # the representations with a consistent compare function. The first from
        # the sorted list of representations is finally taken

        def canonical_from(starting_atom):
            todo_atoms = set(self.nodes)
            todo_atoms.discard(starting_atom)
            boundary_atoms = set([starting_atom])

            # make a tree-representation of the molecule that starts
            # from the starting atom. This tree representation does
            # not contain atom indexes and is therefore independent
            # of any atom order.
            tree = AtomTreeNode(self.numbers[starting_atom], starting_atom)
            done_atoms = {starting_atom: tree}
            while len(todo_atoms) > 0:
                new_atoms = {}
                for b in boundary_atoms:
                    for n in self.neighbors[b]:
                        if n in todo_atoms and n not in new_atoms:
                            new_atoms[n] = AtomTreeNode(self.numbers[n], n)
                for n, node in new_atoms.iteritems():
                    for neighbor in self.neighbors[n]:
                        neighbor_node = done_atoms.get(neighbor)
                        if neighbor_node is not None:
                            neighbor_node.add_child(node)
                for n, node in new_atoms.iteritems():
                    for neighbor in self.neighbors[n]:
                        neighbor_node = new_atoms.get(neighbor)
                        if neighbor_node is not None:
                            neighbor_node.add_sibling(node)
                done_atoms.update(new_atoms)
                boundary_atoms = set(new_atoms)
                del new_atoms
                todo_atoms -= boundary_atoms

            # Now we order the lists inside the three based on
            # the topology of the tree.
            tree.sort()

            # We iterate over the sorted tree and hence get canonically
            # ordered atoms. They only thing that determines the order
            # is the starting atom.
            collected_atoms = set([starting_atom])
            sorted_atoms = {}
            for node in tree.yield_children():
                i = node.index
                if i not in collected_atoms:
                    l = sorted_atoms.setdefault(node.depth,[])
                    l.append(i) # make sure we respect the depth order
                    collected_atoms.add(i)

            sorted_atoms = sum((l for depth, l in sorted(sorted_atoms.iteritems())), [starting_atom])
            # Now get the atom pairs
            sorted_index = dict((old,new) for new,old in enumerate(sorted_atoms))
            atom_numbers = list(self.numbers[i] for i in sorted_atoms)
            bond_pairs = []
            for i,j in self.pairs:
                i,j = sorted_index[i], sorted_index[j]
                if i > j:
                    i,j = j,i
                bond_pairs.append((i,j))
            bond_pairs.sort()
            return atom_numbers, bond_pairs#, numpy.array(sorted_atoms)

        solutions = []
        for starting_atom in starting_atoms:
            solutions.append(canonical_from(starting_atom))
        solutions.sort()

        return solutions[0]



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


class AtomTreeNode(object):
    """Defines a tree structure of atoms.

    This class is used for the generation of a canonical graph.
    """
    def __init__(self, number, index):
        self.number = number
        self.index = index
        self.children = []
        self.parents = []
        self.siblings = []
        #self.num_parents = 0
        #self.num_siblings = 0
        self.depth = 0
        self.tag = 0

    num_siblings = property(lambda self: len(self.siblings))
    num_parents = property(lambda self: len(self.parents))

    def add_child(self, other_node):
        self.children.append(other_node)
        other_node.parents.append(self)
        other_node.depth = self.depth+1

    def add_sibling(self, other_node):
        self.siblings.append(other_node)

    def sort(self):
        for child in self.children:
            child.sort()
        self.children.sort()
        self.siblings.sort()

    def __cmp__(self, other):
        result = -cmp(
            (self.number, self.num_parents, self.num_siblings, len(self.children)),
            (other.number, other.num_parents, other.num_siblings, len(other.children))
        )
        if result != 0:
            return result
        for self_child, other_child in zip(self.children, other.children):
            result = self_child.__cmp__(other_child)
            if result != 0:
                return result
        if self.num_parents > 1:
            #not (self in other.siblings):
            result = cmp(self.index, other.index)
            if result != 0:
                return result
        for self_sibling, other_sibling in zip(self.siblings, other.siblings):
            if self_sibling.index > self.index:
                self_tag = (self.index, self_sibling.index)
            else:
                self_tag = (self_sibling.index, self.index)
            if other_sibling.index > other.index:
                other_tag = (other.index, other_sibling.index)
            else:
                other_tag = (other_sibling.index, other.index)
            result = cmp(self_tag, other_tag)
            if result != 0:
                return result

        return 0

    def yield_children(self):
        for child in self.children:
            yield child
        for child in self.children:
            for v in child.yield_children():
                yield v


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





