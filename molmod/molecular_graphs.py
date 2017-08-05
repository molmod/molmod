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
"""Extension of the graphs module with molecular features"""


from __future__ import division

from builtins import range
import numpy as np

from molmod.graphs import cached, Graph, CustomPattern
from molmod.utils import ReadOnlyAttribute
from molmod.binning import PairSearchIntra


__all__ = [
    "MolecularGraph",
    "HasAtomNumber", "HasNumNeighbors", "HasNeighborNumbers", "HasNeighbors",
    "BondLongerThan", "atom_criteria",
    "BondPattern", "BendingAnglePattern", "DihedralAnglePattern",
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
    def _check_numbers(self, numbers):
        """the number of vertices and atomic numbers must be the same"""
        if len(numbers) != self.num_vertices:
            raise TypeError("The number of vertices must match the length of "
                "the atomic numbers array.")

    def _check_orders(self, orders):
        """the number of edges and bond orders must be the same"""
        if len(orders) != self.num_edges:
            raise TypeError("The number of edges must match the length of "
                "the bond orders array.")

    def _check_symbols(self, symbols):
        """the size must be the same as the length of the array numbers and all elements must be strings"""
        if len(symbols) != self.size:
            raise TypeError("The number of symbols in the graph does not "
                "match the length of the atomic numbers array.")
        for symbol in symbols:
            if not isinstance(symbol, str):
                raise TypeError("All symbols must be strings.")

    numbers = ReadOnlyAttribute(np.ndarray, none=False, check=_check_numbers,
        npdim=1, npdtype=int, doc="the atomic numbers associated with the "
        "vertices")
    orders = ReadOnlyAttribute(np.ndarray, none=False, check=_check_orders,
        npdim=1, npdtype=float, doc="the bond orders associated with the edges")
    symbols = ReadOnlyAttribute(tuple, _check_symbols, doc="symbols for the "
        "atoms, which can be element names for force-field atom types")

    @classmethod
    def from_geometry(cls, molecule, do_orders=False, scaling=1.0):
        """Construct a MolecularGraph object based on interatomic distances

           All short distances are computed with the binning module and compared
           with a database of bond lengths. Based on this comparison, bonded
           atoms are detected.

           Argument:
            | ``molecule``  --  The molecule to derive the graph from

           Optional arguments:
            | ``do_orders``  --  set to True to estimate the bond order
            | ``scaling``  --  scale the threshold for the connectivity. increase
                               this to 1.5 in case of transition states when a
                               fully connected topology is required.
        """
        from molmod.bonds import bonds

        unit_cell = molecule.unit_cell
        pair_search = PairSearchIntra(
            molecule.coordinates,
            bonds.max_length*bonds.bond_tolerance,
            unit_cell
        )

        orders = []
        lengths = []
        edges = []

        for i0, i1, delta, distance in pair_search:
            bond_order = bonds.bonded(molecule.numbers[i0], molecule.numbers[i1], distance/scaling)
            if bond_order is not None:
                if do_orders:
                    orders.append(bond_order)
                lengths.append(distance)
                edges.append((i0,i1))

        if do_orders:
            result = cls(edges, molecule.numbers, orders, symbols=molecule.symbols)
        else:
            result = cls(edges, molecule.numbers, symbols=molecule.symbols)

        # run a check on all neighbors. if two bonds point in a direction that
        # differs only by 45 deg. the longest of the two is discarded. the
        # double loop over the neighbors is done such that the longest bonds
        # are eliminated first
        slated_for_removal = set([])
        threshold = 0.5**0.5
        for c, ns in result.neighbors.items():
            lengths_ns = []
            for n in ns:
                delta = molecule.coordinates[n] - molecule.coordinates[c]
                if unit_cell is not None:
                    delta = unit_cell.shortest_vector(delta)
                length = np.linalg.norm(delta)
                lengths_ns.append([length, delta, n])
            lengths_ns.sort(reverse=True, key=(lambda r: r[0]))
            for i0, (length0, delta0, n0) in enumerate(lengths_ns):
                for i1, (length1, delta1, n1) in enumerate(lengths_ns[:i0]):
                    if length1 == 0.0:
                        continue
                    cosine = np.dot(delta0, delta1)/length0/length1
                    if cosine > threshold:
                        # length1 > length0
                        slated_for_removal.add((c,n1))
                        lengths_ns[i1][0] = 0.0
        # construct a mask
        mask = np.ones(len(edges), bool)
        for i0, i1 in slated_for_removal:
            edge_index = result.edge_index.get(frozenset([i0,i1]))
            if edge_index is None:
                raise ValueError('Could not find edge that has to be removed: %i %i' % (i0, i1))
            mask[edge_index] = False
        # actual removal
        edges = [edges[i] for i in range(len(edges)) if mask[i]]
        if do_orders:
            bond_order = [bond_order[i] for i in range(len(bond_order)) if mask[i]]
            result = cls(edges, molecule.numbers, orders)
        else:
            result = cls(edges, molecule.numbers)

        lengths = [lengths[i] for i in range(len(lengths)) if mask[i]]
        result.bond_lengths = np.array(lengths)

        return result

    @classmethod
    def from_blob(cls, s):
        """Construct a molecular graph from the blob representation"""
        atom_str, edge_str = s.split()
        numbers = np.array([int(s) for s in atom_str.split(",")])
        edges = []
        orders = []
        for s in edge_str.split(","):
            i, j, o = (int(w) for w in s.split("_"))
            edges.append((i, j))
            orders.append(o)
        return cls(edges, numbers, np.array(orders))

    def __init__(self, edges, numbers, orders=None, symbols=None, num_vertices=None):
        """
           Arguments:
            | ``edges``  --  See base class (Graph) documentation
            | ``numbers``  --  consecutive atom numbers

           Optional arguments:
            | ``orders``  --  bond orders
            | ``symbols``  --  atomic symbols
            | ``num_vertices``  --  must be the same as the number of atoms or None

           When the nature of an atom or a bond is unclear ambiguous, set the
           corresponding integer to zero. This means the nature of the atom or
           bond is unspecified. When the bond orders are not given, they are all
           set to  zero.

           If you want to use 'special' atom types, use negative numbers. The
           same for bond orders. e.g. a nice choice for the bond order of a
           hybrid bond is -1.
        """
        if num_vertices is not None and num_vertices != len(numbers):
            raise ValueError('The number of vertices must be the same as the number of atoms.')
        if orders is None:
            orders = np.ones(len(edges), float)
        Graph.__init__(self, edges, len(numbers))
        self.numbers = numbers
        self.orders = orders
        self.symbols = symbols

    def __mul__(self, repeat):
        """Construct a graph that repeats this graph a number of times

           Arguments:
            | ``repeat`` -- The number of repetitions.
        """
        if not isinstance(repeat, int):
            raise TypeError("Can only multiply a graph with an integer")
        # copy edges
        new_edges = []
        for i in range(repeat):
            for vertex1, vertex2 in self.edges:
                new_edges.append(frozenset([vertex1+i*self.num_vertices, vertex2+i*self.num_vertices]))
        # copy numbers
        new_numbers = np.zeros((repeat, len(self.numbers)), int)
        new_numbers[:] = self.numbers
        new_numbers = new_numbers.ravel()
        # copy orders
        new_orders = np.zeros((repeat, len(self.orders)), int)
        new_orders[:] = self.orders
        new_orders = new_orders.ravel()
        # copy symbols
        if self.symbols is not None:
            new_symbols = self.symbols*repeat
        else:
            new_symbols = None
        return MolecularGraph(new_edges, new_numbers, new_orders, new_symbols)

    __rmul__ = __mul__

    @cached
    def blob(self):
        """A compact text representation of the graph"""
        atom_str = ",".join(str(number) for number in self.numbers)
        edge_str = ",".join("%i_%i_%i" % (i, j, o) for (i, j), o in zip(self.edges, self.orders))
        return "%s %s" % (atom_str, edge_str)

    def get_vertex_string(self, i):
        """Return a string based on the atom number"""
        number = self.numbers[i]
        if number == 0:
            return Graph.get_vertex_string(self, i)
        else:
            # pad with zeros to make sure that string sort is identical to number sort
            return "%03i" % number

    def get_edge_string(self, i):
        """Return a string based on the bond order"""
        order = self.orders[i]
        if order == 0:
            return Graph.get_edge_string(self, i)
        else:
            # pad with zeros to make sure that string sort is identical to number sort
            return "%03i" % order

    def get_subgraph(self, subvertices, normalize=False):
        """Creates a subgraph of the current graph

           See :meth:`molmod.graphs.Graph.get_subgraph` for more information.
        """
        graph = Graph.get_subgraph(self, subvertices, normalize)
        if normalize:
            new_numbers = self.numbers[graph._old_vertex_indexes] # vertices do change
        else:
            new_numbers = self.numbers # vertices don't change!
        if self.symbols is None:
            new_symbols = None
        elif normalize:
            new_symbols = tuple(self.symbols[i] for i in graph._old_vertex_indexes)
        else:
            new_symbols = self.symbols
        new_orders = self.orders[graph._old_edge_indexes]
        result = MolecularGraph(graph.edges, new_numbers, new_orders, new_symbols)
        if normalize:
            result._old_vertex_indexes = graph._old_vertex_indexes
        result._old_edge_indexes = graph._old_edge_indexes
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

        new_edges = list(self.edges)
        counter = self.num_vertices
        for i in range(self.num_vertices):
            num_elec = self.numbers[i]
            if formal_charges is not None:
                num_elec -= int(formal_charges[i])
            if num_elec >= 5 and num_elec <= 9:
                num_hydrogen = num_elec - 10 + 8
            elif num_elec >= 13 and num_elec <= 17:
                num_hydrogen = num_elec - 18 + 8
            elif num_elec == 35:
                num_hydrogen = 1
            else:
                continue
            if num_hydrogen > 4:
                num_hydrogen = 8 - num_hydrogen
            for n in self.neighbors[i]:
                bo = self.orders[self.edge_index[frozenset([i, n])]]
                if bo <= 0:
                    bo = 1
                num_hydrogen -= int(bo)
            for j in range(num_hydrogen):
                new_edges.append((i, counter))
                counter += 1
        new_numbers = np.zeros(counter, int)
        new_numbers[:self.num_vertices] = self.numbers
        new_numbers[self.num_vertices:] = 1
        new_orders = np.zeros(len(new_edges), int)
        new_orders[:self.num_edges] = self.orders
        new_orders[self.num_edges:] = 1
        result = MolecularGraph(new_edges, new_numbers, new_orders)
        return result




# basic criteria for molecular patterns

class HasAtomNumber(object):
    """Criterion for the atom number of a vertex"""

    def __init__(self, number):
        """
           Arguments:
            | ``number``  --  the expected atom number
        """
        self.number = number

    def __call__(self, index, graph):
        """Return True only if the atom number is correct

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion is
                         applied
            | ``graph``  --  the graph on which the criterion is tested
        """
        return graph.numbers[index] == self.number


class HasNumNeighbors(object):
    """Criterion for the number of neighboring vertexes"""

    def __init__(self, count):
        """
           Arguments:
            | ``count``  --  the expected number of neighbors
        """
        self.count = count

    def __call__(self, index, graph):
        """Return True only if the number of neighbors is correct

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion is
                             applied
            | ``graph``  --  the graph on which the criterion is tested
        """
        return len(graph.neighbors[index]) == self.count


class HasNeighborNumbers(object):
    """Criterion for the atom numbers of the neighbor vertexes"""

    def __init__(self, *numbers):
        """
           Arguments:
            | ``*numbers``  --  a list with atom numbers
        """
        for number in numbers:
            if not isinstance(number, int):
                raise TypeError("All arguments must be integers, found a %s" % type(number))
        self.numbers = list(numbers)
        self.numbers.sort()

    def __call__(self, index, graph):
        """Return True only if each neighbor can be linked with an atom number

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion is
                             applied
            | ``graph``  --  the graph on which the criterion is tested
        """
        neighbors = graph.neighbors[index]
        if not len(neighbors) == len(self.numbers):
            return False
        neighbor_numbers = sorted([graph.numbers[neighbor] for neighbor in neighbors])
        return neighbor_numbers == self.numbers


class HasNeighbors(object):
    """Tests if the neighbors of a vertex match the given criteria"""
    def __init__(self, *neighbor_criteria):
        """
           Arguments:
            | ``*neighbor_criteria``  --  a list of criteria objects
        """
        for criterion in neighbor_criteria:
            if not hasattr(criterion, "__call__"):
                raise TypeError("All arguments must be callable")
        self.neighbor_criteria = list(neighbor_criteria)

    def __call__(self, index, graph):
        """Return True only if each neighbor can be linked with a positive criterion

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion is
                             applied
            | ``graph``  --  the graph on which the criterion is tested
        """

        def all_permutations(l):
            """Iterate over all permutations"""
            if len(l) == 1:
                yield l
                return
            for i in range(len(l)):
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
        """
           This criterion assumes that the molecular graph has an attribute
           self.bond_lengths

           Argument:
            | ``length`` -- the minimum length of the bond
        """
        self.length = length

    def __call__(self, index, graph):
        """Return True only if the bond is longer than the threshold

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion is
                             applied
            | ``graph``  --  the graph on which the criterion is tested
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


class BondPattern(CustomPattern):
    """Pattern for two consecutive vertices"""
    def __init__(self, criteria_sets=None, vertex_tags=None):
        """
           Arguments: see :class:`molmod.graphs.CustomPattern`
        """
        if vertex_tags is None:
            vertex_tags = {}
        pattern_graph = Graph([(0, 1)])
        CustomPattern.__init__(self, pattern_graph, criteria_sets, vertex_tags)


class BendingAnglePattern(CustomPattern):
    """Pattern for three consecutive vertices"""
    def __init__(self, criteria_sets=None, vertex_tags=None):
        """
           Arguments: see :class:`molmod.graphs.CustomPattern`
        """
        if vertex_tags is None:
            vertex_tags = {}
        pattern_graph = Graph([(0, 1), (1, 2)])
        CustomPattern.__init__(self, pattern_graph, criteria_sets, vertex_tags)


class DihedralAnglePattern(CustomPattern):
    """Pattern for four consecutive vertices"""
    def __init__(self, criteria_sets=None, vertex_tags=None):
        """
           Arguments: see :class:`molmod.graphs.CustomPattern`
        """
        if vertex_tags is None:
            vertex_tags = {}
        pattern_graph = Graph([(0, 1), (1, 2), (2, 3)])
        CustomPattern.__init__(self, pattern_graph, criteria_sets, vertex_tags)


class OutOfPlanePattern(CustomPattern):
    """Pattern for a central vertex connected to three other vertices"""
    def __init__(self, criteria_sets=None, vertex_tags=None):
        """
           Arguments: see :class:`molmod.graphs.CustomPattern`
        """
        if vertex_tags is None:
            vertex_tags = {}
        pattern_graph = Graph([(0, 1), (0, 2), (0, 3)])
        CustomPattern.__init__(self, pattern_graph, criteria_sets, vertex_tags)


class TetraPattern(CustomPattern):
    """Pattern for a central vertex connected to four other vertices"""
    def __init__(self, criteria_sets=None, vertex_tags=None):
        """
           Arguments: see :class:`molmod.graphs.CustomPattern`
        """
        if vertex_tags is None:
            vertex_tags = {}
        pattern_graph = Graph([(0, 1), (0, 2), (0, 3), (0, 4)])
        CustomPattern.__init__(self, pattern_graph, criteria_sets, vertex_tags)


class NRingPattern(CustomPattern):
    """Pattern for strong rings with a fixed size"""

    def __init__(self, size, criteria_sets=None, vertex_tags=None, strong=False):
        """
           Argument:
            | ``size``  --  the size of the ring
        """
        if vertex_tags is None:
            vertex_tags = {}
        self.size = size
        self.strong = strong
        pattern_graph = Graph([(i, (i+1)%size) for i in range(size)])
        CustomPattern.__init__(self, pattern_graph, criteria_sets, vertex_tags)

    def check_next_match(self, match, new_relations, subject_graph, one_match):
        """Check if the (onset for a) match can be a valid (part of a) ring"""
        if not CustomPattern.check_next_match(self, match, new_relations, subject_graph, one_match):
            return False
        if self.strong:
            # can this ever become a strong ring?
            vertex1_start = match.forward[self.pattern_graph.central_vertex]
            for vertex1 in new_relations.values():
                paths = list(subject_graph.iter_shortest_paths(vertex1, vertex1_start))
                if self.size % 2 == 0 and len(match) == self.size:
                    if len(paths) != 2:
                        #print "NRingPattern.check_next_match: not strong a.1"
                        return False
                    for path in paths:
                        if len(path) != len(match)//2+1:
                            #print "NRingPattern.check_next_match: not strong a.2"
                            return False
                else:
                    if len(paths) != 1:
                        #print "NRingPattern.check_next_match: not strong b.1"
                        return False
                    if len(paths[0]) != (len(match)+1)//2:
                        #print "NRingPattern.check_next_match: not strong b.2"
                        return False
            #print "RingPattern.check_next_match: no remarks"
        return True

    def complete(self, match, subject_graph):
        """Check the completeness of the ring match"""
        if not CustomPattern.complete(self, match, subject_graph):
            return False
        if self.strong:
            # If the ring is not strong, return False
            if self.size % 2 == 0:
                # even ring
                for i in range(self.size//2):
                    vertex1_start = match.forward[i]
                    vertex1_stop = match.forward[i+self.size//2]
                    paths = list(subject_graph.iter_shortest_paths(vertex1_start, vertex1_stop))
                    if len(paths) != 2:
                        #print "Even ring must have two paths between opposite vertices"
                        return False
                    for path in paths:
                        if len(path) != self.size//2+1:
                            #print "Paths between opposite vertices must half the size of the ring+1"
                            return False
            else:
                # odd ring
                for i in range(self.size//2+1):
                    vertex1_start = match.forward[i]
                    vertex1_stop = match.forward[i+self.size//2]
                    paths = list(subject_graph.iter_shortest_paths(vertex1_start, vertex1_stop))
                    if len(paths) > 1:
                        return False
                    if len(paths[0]) != self.size//2+1:
                        return False
                    vertex1_stop = match.forward[i+self.size//2+1]
                    paths = list(subject_graph.iter_shortest_paths(vertex1_start, vertex1_stop))
                    if len(paths) > 1:
                        return False
                    if len(paths[0]) != self.size//2+1:
                        return False
        return True
