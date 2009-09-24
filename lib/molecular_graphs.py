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
"""Extension of the graphs module with molecular features"""


from molmod.graphs import cached, Graph, SubgraphPattern
from molmod.binning import PairSearch

import numpy


__all__ = [
    "MolecularGraph",
    "HasAtomNumber", "HasNumNeighbors", "HasNeighborNumbers", "HasNeighbors", "BondLongerThan",
    "atom_criteria", "BondPattern", "BendingAnglePattern", "DihedralAnglePattern",
    "OutOfPlanePattern", "TetraPattern", "NRingPattern",
]


class MolecularGraph(Graph):
    """Describes a molecular graph: connectivity, atom numbers and bond orders.

       This class inherits all features from the Graph class and adds methods
       and attributes that are specific from molecular graphs. Instances are
       immutable, so if you want to modify a molecular graph, just instantiate a
       new object with modified connectivity, numbers and orders. The advantage
       is that various graph analysis and properties can be cached.
    """

    @classmethod
    def from_geometry(cls, molecule, unit_cell=None, do_orders=False):
        """Construct a MolecularGraph object based on interatomic distances

           All short distances are computed with the binning module and compared
           with a database of bond lengths. Based on this comparison, bonded
           atoms are detected.
        """
        from molmod.bonds import bonds

        pair_search = PairSearch(
            molecule.coordinates,
            bonds.max_length*bonds.bond_tolerance,
            unit_cell
        )

        orders = []
        lengths = []
        pairs = []

        for i0, i1, delta, distance in pair_search:
            bond_order = bonds.bonded(molecule.numbers[i0], molecule.numbers[i1], distance)
            if bond_order is not None:
                if do_orders:
                    orders.append(bond_order)
                lengths.append(distance)
                pairs.append((i0,i1))

        if do_orders:
            result = cls(pairs, molecule.numbers, orders)
        else:
            result = cls(pairs, molecule.numbers)
        result.bond_lengths = numpy.array(lengths)

        return result

    @classmethod
    def from_blob(cls, s):
        """Construct a molecular graph from the blob representation"""
        atom_str, pair_str = s.split()
        numbers = numpy.array([int(s) for s in atom_str.split(",")])
        pairs = []
        orders = []
        for s in pair_str.split(","):
            i, j, o = (int(w) for w in s.split("_"))
            pairs.append((i, j))
            orders.append(o)
        return cls(pairs, numbers, numpy.array(orders))

    def __init__(self, pairs, numbers, orders=None):
        """Initialize a molecular graph

           Arguments:
             pairs -- See base class (Graph) documentation
             numbers -- consecutive atom numbers
             orders -- bond orders

           When the nature of an atom or a bond is unclear ambiguous, set the
           corresponding integer to zero. This means the nature of the atom or
           bond is unspecified. When the bond orders are not given, they are all
           set to  zero.

           If you want to use 'special' atom types, use negative numbers. The
           same for bond orders. e.g. a nice choice for the bond order of a
           hybrid bond is -1.
        """
        if orders is None:
            orders = numpy.ones(len(pairs), int)
        elif len(orders) != len(pairs):
            raise ValueError("The number of (bond) orders must be equal to the number of pairs")
        Graph.__init__(self, pairs, len(numbers))
        self._init_attributes({"numbers": numbers, "orders": orders}, {})

    def __mul__(self, repeat):
        """Construct a graph that repeats this graph a number of times

           Arguments:
             repeat -- The number of repetitions.
        """
        if not isinstance(repeat, int):
            raise TypeError("Can only multiply a graph with an integer")
        # copy pairs
        new_pairs = []
        for i in xrange(repeat):
            for node1, node2 in self.pairs:
                new_pairs.append(frozenset([node1+i*self.num_nodes, node2+i*self.num_nodes]))
        # copy numbers
        new_numbers = numpy.zeros((repeat, len(self.numbers)), int)
        new_numbers[:] = self.numbers
        new_numbers = new_numbers.ravel()
        # copy orders
        new_orders = numpy.zeros((repeat, len(self.orders)), int)
        new_orders[:] = self.orders
        new_orders = new_orders.ravel()
        return MolecularGraph(new_pairs, new_numbers, new_orders)

    __rmul__ = __mul__

    @cached
    def blob(self):
        """Create a compact text representation of the graph"""
        atom_str = ",".join(str(number) for number in self.numbers)
        pair_str = ",".join("%i_%i_%i" % (i, j, o) for (i, j), o in zip(self.pairs, self.orders))
        return "%s %s" % (atom_str, pair_str)

    def get_node_string(self, i):
        """Return a string based on the atom number"""
        number = self.numbers[i]
        if number == 0:
            return Graph.get_node_string(self, i)
        else:
            # pad with zeros to make sure that string sort is identical to number sort
            return "%03i" % number

    def get_pair_string(self, i):
        """Return a string based on the bond order"""
        order = self.orders[i]
        if order == 0:
            return Graph.get_pair_string(self, i)
        else:
            # pad with zeros to make sure that string sort is identical to number sort
            return "%03i" % order

    def get_subgraph(self, subnodes, normalize=False):
        """Creates a subgraph of the current graph

           See help(Graph.get_subgraph) for more information.
        """
        graph = Graph.get_subgraph(self, subnodes, normalize)
        if normalize:
            new_numbers = self.numbers[graph._old_node_indexes] # nodes do change
        else:
            new_numbers = self.numbers # nodes don't change!
        new_orders = self.orders[graph._old_pair_indexes]
        result = MolecularGraph(graph.pairs, new_numbers, new_orders)
        if normalize:
            result._old_node_indexes = graph._old_node_indexes
        result._old_pair_indexes = graph._old_pair_indexes
        return result

    def add_hydrogens(self, formal_charges=None):
        """Returns a molecular graph where hydrogens are added explicitely

           When the bond order is unknown, it assumes bond order one. If the
           graph  has an attribute formal_charges, this routine will take it
           into account when counting the number of hydrogens to be added. The
           returned graph will also have a formal_charges attribute.

           This routine only adds hydrogen atoms for a limited set of atoms from
           the periodic system: B, C, N, O, F, Al, Si, P, S, Cl, Br.
        """

        new_pairs = list(self.pairs)
        counter = self.num_nodes
        for i in xrange(self.num_nodes):
            num_elec = self.numbers[i]
            if formal_charges is not None:
                num_elec -= formal_charges[i]
            if num_elec >= 5 and num_elec <= 9:
                num_hydrogen = num_elec - 10 + 8
            elif num_elec >= 13 and num_elec <= 17:
                num_hydrogen = num_elec - 18 + 8
            elif num_elec == 35:
                num_hydrogen = 1
            else:
                continue
            if num_hydrogen > 4:
                num_hydrogen = 8-num_hydrogen
            for n in self.neighbors[i]:
                bo = self.orders[self.pair_index[frozenset([i, n])]]
                if bo <= 0:
                    bo = 1
                num_hydrogen -= bo
            for j in xrange(num_hydrogen):
                new_pairs.append((i, counter))
                counter += 1
        new_numbers = numpy.zeros(counter, int)
        new_numbers[:self.num_nodes] = self.numbers
        new_numbers[self.num_nodes:] = 1
        new_orders = numpy.zeros(len(new_pairs), int)
        new_orders[:self.num_pairs] = self.orders
        new_orders[self.num_pairs:] = 1
        result = MolecularGraph(new_pairs, new_numbers, new_orders)
        return result




# basic criteria for molecular patterns

class HasAtomNumber(object):
    """Criterion for the atom number of a vertex"""

    def __init__(self, number):
        """Initialize a HasAtomNumber object

           Arguments:
             number  --  the expected atom number
        """
        self.number = number

    def __call__(self, index, graph):
        """Return True only if the atom number is correct

           Arguments:
             index  --  the index of the vertex/edge on which the criterion is
                        applied
             graph  --  the graph on which the criterion is tested
        """
        return graph.numbers[index] == self.number


class HasNumNeighbors(object):
    """Criterion for the number of neighboring vertexes"""

    def __init__(self, count):
        """Initialize a HasNumNeighbors object

           Arguments:
             count  --  the expected number of neighbors
        """
        self.count = count

    def __call__(self, index, graph):
        """Return True only if the number of neighbors is correct

           Arguments:
             index  --  the index of the vertex/edge on which the criterion is
                        applied
             graph  --  the graph on which the criterion is tested
        """
        return len(graph.neighbors[index]) == self.count


class HasNeighborNumbers(object):
    """Criterion for the atom numbers of the neighbor vertexes"""

    def __init__(self, *numbers):
        """Initialize a HasNeighborNumbers object

           Arguments:
             *numbers  --  a list with atom numbers
        """
        self.numbers = list(numbers)
        self.numbers.sort()

    def __call__(self, index, graph):
        """Return True only if each neighbor can be linked with an atom number

           Arguments:
             index  --  the index of the vertex/edge on which the criterion is
                        applied
             graph  --  the graph on which the criterion is tested
        """
        neighbors = graph.neighbors[index]
        if not len(neighbors) == len(self.numbers):
            return
        neighbor_numbers = [graph.numbers[neighbor] for neighbor in neighbors]
        neighbor_numbers.sort()
        return neighbor_numbers == self.numbers


class HasNeighbors(object):
    """Tests if the neighbors of a vertex match the given criteria"""
    def __init__(self, *neighbor_criteria):
        """Initialize a HasNeighbors object

           Arguments:
             *neighbor_criteria  --  a list of criteria objects
        """
        self.neighbor_criteria = list(neighbor_criteria)

    def __call__(self, index, graph):
        """Return True only if each neighbor can be linked with a positive criterion

           Arguments:
             index  --  the index of the vertex/edge on which the criterion is
                        applied
             graph  --  the graph on which the criterion is tested
        """

        def all_permutations(l):
            """Iterate over all permutations"""
            if len(l) == 1:
                yield l
                return
            for i in xrange(len(l)):
                for sub in all_permutations(l[:i]+l[i+1:]):
                    yield [l[i]] + sub

        neighbors = graph.neighbors[index]
        if not len(neighbors) == len(self.neighbor_criteria):
            return
        # consider all permutations. If one matches, return True
        for perm_neighbors in all_permutations(list(neighbors)):
            ok = True
            for neighbor, crit in zip(perm_neighbors, self.neighbor_criteria):
                if not crit(neighbor, graph):
                    ok = False
                    break
            if ok:
                return True
        return False


class BondLongerThan(object):
    """A vertex criterion to select bonds longer than a given threshold"""
    def __init__(self, length):
        """Initialize a BondLongerThan object

           This criterion assumes that the molecular graph has an attribute
           self.bond_lengths

           Argument:
             length -- the minimum length of the bond
        """
        self.length = length

    def __call__(self, index, graph):
        """Return True only if the bond is longer than the threshold

           Arguments:
             index  --  the index of the vertex/edge on which the criterion is
                        applied
             graph  --  the graph on which the criterion is tested
        """
        return graph.bond_lengths[index] > self.length


def atom_criteria(*params):
    """An auxiliary function to construct a dictionary of Criteria"""
    result = {}
    for index, param in enumerate(params):
        if param is None:
            continue
        elif isinstance(param, int):
            result[index] = HasAtomNumber(param)
        else:
            result[index] = param
    return result


# common patterns for molecular structures


class BondPattern(SubgraphPattern):
    """Pattern for two consecutive vertices"""
    def __init__(self, criteria_sets=None, node_tags=None):
        """Initialize a BondPattern object

           Arguments: see SubgraphPattern.__init__
        """
        if node_tags is None:
            node_tags = {}
        subgraph = Graph([(0, 1)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class BendingAnglePattern(SubgraphPattern):
    """Pattern for three consecutive vertices"""
    def __init__(self, criteria_sets=None, node_tags=None):
        """Initialize a BendingAnglePattern object

           Arguments: see SubgraphPattern.__init__
        """
        if node_tags is None:
            node_tags = {}
        subgraph = Graph([(0, 1), (1, 2)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class DihedralAnglePattern(SubgraphPattern):
    """Pattern for four consecutive vertices"""
    def __init__(self, criteria_sets=None, node_tags=None):
        """Initialize a DihedralAnglePattern object

           Arguments: see SubgraphPattern.__init__
        """
        if node_tags is None:
            node_tags = {}
        subgraph = Graph([(0, 1), (1, 2), (2, 3)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class OutOfPlanePattern(SubgraphPattern):
    """Pattern for a central vertex connected to three other vertices"""
    def __init__(self, criteria_sets=None, node_tags=None):
        """Initialize a TetraPattern object

           Arguments: see SubgraphPattern.__init__
        """
        if node_tags is None:
            node_tags = {}
        subgraph = Graph([(0, 1), (0, 2), (0, 3)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class TetraPattern(SubgraphPattern):
    """Pattern for a central vertex connected to four other vertices"""
    def __init__(self, criteria_sets=None, node_tags=None):
        """Initialize a TetraPattern object

           Arguments: see SubgraphPattern.__init__
        """
        if node_tags is None:
            node_tags = {}
        subgraph = Graph([(0, 1), (0, 2), (0, 3), (0, 4)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class NRingPattern(SubgraphPattern):
    """Pattern for strong rings with a fixed size"""

    def __init__(self, size, criteria_sets=None, node_tags=None, strong=False):
        """Initialize a NRingPattern object

           Argument:
             size  --  the size of the ring
        """
        if node_tags is None:
            node_tags = {}
        self.size = size
        self.strong = strong
        subgraph = Graph([(i, (i+1)%size) for i in xrange(size)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)

    def check_next_match(self, match, new_relations):
        """Check if the (onset for a) match can be a valid (part of a) ring"""
        if not SubgraphPattern.check_next_match(self, match, new_relations):
            return False
        if self.strong:
            # can this ever become a strong ring?
            node1_start = match.forward[self.subgraph.central_node]
            for node1 in new_relations.itervalues():
                paths = list(self.graph.iter_shortest_paths(node1, node1_start))
                if self.size % 2 == 0 and len(match) == self.size:
                    if len(paths) != 2:
                        #print "NRingPattern.check_next_match: not strong a.1"
                        return False
                    for path in paths:
                        if len(path) != len(match)/2+1:
                            #print "NRingPattern.check_next_match: not strong a.2"
                            return False
                else:
                    if len(paths) != 1:
                        #print "NRingPattern.check_next_match: not strong b.1"
                        return False
                    if len(paths[0]) != (len(match)+1)/2:
                        #print "NRingPattern.check_next_match: not strong b.2"
                        return False
            #print "RingPattern.check_next_match: no remarks"
        return True

    def complete(self, match):
        """Check the completeness of the ring match"""
        if not SubgraphPattern.complete(self, match):
            return False
        if self.strong:
            # If the ring is not strong, return False
            if self.size % 2 == 0:
                # even ring
                for i in xrange(self.size/2):
                    node1_start = match.forward[i]
                    node1_stop = match.forward[i+self.size/2]
                    paths = list(self.graph.iter_shortest_paths(node1_start, node1_stop))
                    if len(paths) != 2:
                        #print "Even ring must have two paths between opposite nodes"
                        return False
                    for path in paths:
                        if len(path) != self.size/2+1:
                            #print "Paths between opposite nodes must half the size of the ring+1"
                            return False
            else:
                # odd ring
                for i in xrange(self.size/2+1):
                    node1_start = match.forward[i]
                    node1_stop = match.forward[i+self.size/2]
                    paths = list(self.graph.iter_shortest_paths(node1_start, node1_stop))
                    if len(paths) > 1:
                        return False
                    if len(paths[0]) != self.size/2+1:
                        return False
                    node1_stop = match.forward[i+self.size/2+1]
                    paths = list(self.graph.iter_shortest_paths(node1_start, node1_stop))
                    if len(paths) > 1:
                        return False
                    if len(paths[0]) != self.size/2+1:
                        return False
        return True

