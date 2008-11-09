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
from molmod.graphs import cached, Graph, SubgraphPattern, EqualPattern, Match, OneToOne, GraphSearch, CriteriaSet
from molmod.binning import IntraAnalyseNeighboringObjects, PositionedObject, SparseBinnedObjects
from molmod.units import angstrom

from molmod.ext import ff_dm_quad, ff_dm_reci, ff_bond_quad, ff_bond_hyper

import numpy, copy


__all__ = [
    "MolecularGraph", "ToyFF",
    "HasAtomNumber", "HasNumNeighbors", "HasNeighborNumbers", "BondLongerThan",
    "atom_criteria", "BondPattern", "BendingAnglePattern", "DihedralAnglePattern",
    "OutOfPlanePattern", "TetraPattern", "NRingPattern",
]


class MolecularGraph(Graph):
    """Describes a molecular graph: connectivity, atom numbers and bond orders.

    This class inherits all features from the Graph class and adds methods and
    attributes that are specific from molecular graphs. Instances are immutable,
    so if you want to modify a molecular graph, just instantiate a new object
    with modified connectivity, numbers and orders. The advantage is that
    various graph analysis and properties can be cached.
    """
    @classmethod
    def from_geometry(cls, molecule, unit_cell=None):
        from molmod.data.bonds import bonds

        def yield_positioned_atoms():
            for index in xrange(molecule.size):
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
            (frozenset([positioned.id for positioned in key]), data)
            for key, data
            in IntraAnalyseNeighboringObjects(binned_atoms, compare_function)(unit_cell)
        )
        pairs = tuple(pair for pair, data in bond_data)
        orders = numpy.array([data[0] for pair, data in bond_data], dtype=int)
        lengths = numpy.array([data[1] for pair, data in bond_data], dtype=float)

        result = cls(pairs, molecule.numbers, orders)
        result.bond_lengths = lengths
        return result

    @classmethod
    def from_blob(cls, s):
        """Construct a molecular graph from the blob representation created with get_blob"""
        atom_str, pair_str = s.split()
        numbers = numpy.array([int(s) for s in atom_str.split(",")])
        pairs = []
        orders = []
        for s in pair_str.split(","):
            i,j,o = (int(w) for w in s.split("_"))
            pairs.append((i,j))
            orders.append(o)
        return cls(pairs,numbers,numpy.array(orders))

    def __init__(self, pairs, numbers, orders=None):
        """Initialize a molecular graph

        Arguments:
          pairs -- See base class (Graph) documentation
          numbers -- consecutive atom numbers
          orders -- bond orders

        When the nature of an atom or a bond is unclear ambiguous, set the
        corresponding integer to zero. This means the nature of the atom or bond
        is unspecified. When the bond orders are not given, they are all set to
        zero.

        If you want to use 'special' atom types, use negative numbers. The same
        for bond orders. e.g. a nice choice for the bond order of a hybrid bond
        is -1.
        """
        if orders is None:
            orders = numpy.zeros(len(pairs), dtype=int)
        elif len(orders) != len(pairs):
            raise ValueError("The number of (bond) orders must be equal to the number of pairs")
        Graph.__init__(self, pairs, len(numbers))
        self._numbers = numbers
        self._numbers.setflags(write=False)
        self._orders = orders
        self._orders.setflags(write=False)

    numbers = property(lambda self: self._numbers)
    orders = property(lambda self: self._orders)

    def __mul__(self, repeat):
        result = Graph.__mul__(self, repeat)
        result.__class__ = MolecularGraph
        # copy numbers
        numbers = numpy.zeros((repeat, len(self.numbers)), int)
        numbers[:] = self.numbers
        result._numbers = numbers.ravel()
        result._numbers.setflags(write=False)
        # copy orders
        orders = numpy.zeros((repeat, len(self.orders)), int)
        orders[:] = self.orders
        result._orders = orders.ravel()
        result._orders.setflags(write=False)
        return result

    __rmul__ = __mul__

    @cached
    def blob(self):
        """Create a compact text representation of the graph."""
        atom_str = ",".join(str(number) for number in self.numbers)
        pair_str = ",".join("%i_%i_%i" % (i,j,o) for (i,j),o in zip(self.pairs,self.orders))
        return "%s %s" % (atom_str, pair_str)

    def get_node_string(self, i):
        number = self.numbers[i]
        if number == 0:
            return Graph.get_node_string(self, i)
        else:
            # pad with zeros to make sure that string sort is identical to number sort
            return "%03i" % number

    def get_pair_string(self, i):
        order = self.orders[i]
        if order == 0:
            return Graph.get_pair_string(self, i)
        else:
            # pad with zeros to make sure that string sort is identical to number sort
            return "%03i" % order

    def get_subgraph(self, subnodes, normalize=False):
        """Creates a subgraph of the current graph.

        See help(Graph.get_subgraph) for more information.
        """
        result = Graph.get_subgraph(self, subnodes, normalize)
        result.__class__ = MolecularGraph
        if normalize:
            result._numbers = self.numbers[result.old_node_indexes] # nodes do change
            result._numbers.setflags(write=False)
        else:
            result._numbers = self.numbers # nodes don't change!
        result._orders = self.orders[result.old_pair_indexes]
        result._orders.setflags(write=False)
        return result

    def guess_geometry(self):
        N = len(self.numbers)
        from molmod.minimizer import Minimizer, NewtonGLineSearch

        ff = ToyFF(self)
        x_init = numpy.random.normal(0,1,N*3)

        #  level 1 geometry optimization
        ff.dm_quad = 1.0
        minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-6, 1e-10, 2*N, 500, do_gradient=True, verbose=False)
        x_init = minimizer.x

        #  level 2 geometry optimization
        ff.dm_quad = 0.0
        ff.dm_reci = 0.1
        ff.bond_quad = 1.0
        minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-4, 1e-10, 2*N, 500, do_gradient=True, verbose=False)
        x_init = minimizer.x

        #  level 3 geometry optimization
        ff.bond_quad = 0.0
        ff.bond_hyper = 1.0
        ff.span_quad = 10.0
        ff.dm_reci = 0.5
        minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-5, 1e-10, 2*N, 500, do_gradient=True, verbose=False)
        x_init = minimizer.x

        x_opt = x_init

        mol = Molecule(self.numbers, x_opt.reshape((N,3)))
        return mol

    def add_hydrogens(self, formal_charges=None):
        """Returns a molecular graph where hydrogens are added explicitely.

        When the bond order is unknown, it assumes bond order one. If the graph
        has an attribute formal_charges, this routine will take it into account
        when counting the number of hydrogens to be added. The returned graph
        will also have a formal_charges attribute.

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
                bo = self.orders[self.pair_index[frozenset([i,n])]]
                if bo <= 0:
                    bo = 1
                num_hydrogen -= bo
            for j in xrange(num_hydrogen):
                new_pairs.append((i,counter))
                counter += 1
        new_numbers = numpy.zeros(counter, int)
        new_numbers[:self.num_nodes] = self.numbers
        new_numbers[self.num_nodes:] = 1
        new_orders = numpy.zeros(len(new_pairs), int)
        new_orders[:self.num_pairs] = self.orders
        new_orders[self.num_pairs:] = 1
        result = MolecularGraph(new_pairs, new_numbers, new_orders)
        return result

# initial geometry based on graph

class ToyFF(object):
    def __init__(self, graph):
        from molmod.data.bonds import bonds

        self.dm0 = graph.distances.astype(numpy.int32)
        self.numbers = graph.numbers.astype(numpy.int32)

        num_bonds = len(graph.pairs)
        bond_pairs = []
        bond_lengths = []
        for counter, (i,j) in enumerate(graph.pairs):
            bond_pairs.append((i,j))
            bond_lengths.append(bonds.get_length(graph.numbers[i],graph.numbers[j]))
        self.bond_pairs = numpy.array(bond_pairs, numpy.int32)
        self.bond_lengths = numpy.array(bond_lengths, float)

        span_pairs = []
        span_lengths = []
        for i, neighbors in graph.neighbors.iteritems():
            number_i = graph.numbers[i]
            if (number_i >= 5 and number_i <=8):
                order = len(neighbors) + abs(number_i-6)
            elif number_i >= 13 and number_i <= 16:
                order = len(neighbors) + abs(number_i-14)
            else:
                order = -1
            if order < 2 or order > 6:
                angle = numpy.pi/180.0*115.0
            elif order == 2:
                angle = numpy.pi
            elif order == 3:
                angle = numpy.pi/180.0*125.0
            elif order == 4:
                angle = numpy.pi/180.0*109.0
            elif order == 5:
                angle = numpy.pi/180.0*100.0
            elif order == 6:
                angle = numpy.pi/180.0*90.0
            for j in neighbors:
                number_j = graph.numbers[j]
                for k in neighbors:
                    if j<k and not frozenset([j,k]) in graph.pairs:
                        number_k = graph.numbers[k]
                        dj = bonds.get_length(number_i,number_j)
                        dk = bonds.get_length(number_i,number_k)
                        d = numpy.sqrt(dj**2+dk**2-2*dj*dk*numpy.cos(angle))
                        span_pairs.append((j,k))
                        span_lengths.append(d)
        self.span_pairs = numpy.array(span_pairs, numpy.int32)
        self.span_lengths = numpy.array(span_lengths, float)

        self.dm_quad = 0.0
        self.dm_reci = 0.0
        self.bond_quad = 0.0
        self.span_quad = 0.0
        self.bond_hyper = 0.0

    def __call__(self, x, do_gradient=False):
        x = x.reshape((-1,3))
        result = 0.0

        gradient = numpy.zeros(x.shape, float)
        if self.dm_quad > 0.0:
            result += ff_dm_quad(self.numbers, x, self.dm0, 1.5*angstrom, self.dm_quad, gradient)
        if self.dm_reci:
            result += ff_dm_reci(self.numbers, x, self.dm0, self.dm_reci, gradient)
        if self.bond_quad:
            result += ff_bond_quad(x, self.bond_pairs, self.bond_lengths, self.bond_quad, gradient)
        if self.span_quad:
            result += ff_bond_quad(x, self.span_pairs, self.span_lengths, self.span_quad, gradient)
        if self.bond_hyper:
            result += ff_bond_hyper(x, self.bond_pairs, self.bond_lengths, 10.0, self.bond_hyper, gradient)
        if do_gradient:
            return result, gradient.ravel()
        else:
            return result


# basic criteria for molecular patterns

class HasAtomNumber(object):
    def __init__(self, number):
        self.number = number

    def __call__(self, atom, graph):
        return graph.numbers[atom] == self.number


class HasNumNeighbors(object):
    def __init__(self, count):
        self.count = count

    def __call__(self, atom, graph):
        return len(graph.neighbors[atom]) == self.count


class HasNeighborNumbers(object):
    def __init__(self, numbers):
        self.numbers = list(numbers)
        self.numbers.sort()

    def __call__(self, atom, graph):
        neighbors = graph.neighbors[atom]
        if not len(neighbors) == len(self.numbers):
            return
        neighbor_numbers = [graph.numbers[neighbor] for neighbor in neighbors]
        neighbor_numbers.sort()
        return neighbor_numbers == self.numbers


class BondLongerThan(object):
    def __init__(self, length):
        self.length = length

    def __call__(self, pair_index, graph):
        return graph.bond_lengths[pair_index] > self.length


def atom_criteria(*params):
    """An auxiliary function to construct a dictionary of Criteria geared
    towards molecular patterns."""
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
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([(0, 1)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class BendingAnglePattern(SubgraphPattern):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([(0, 1), (1, 2)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class DihedralAnglePattern(SubgraphPattern):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([(0, 1), (1, 2), (2, 3)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class OutOfPlanePattern(SubgraphPattern):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([(0, 1), (0, 2), (0, 3)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class TetraPattern(SubgraphPattern):
    def __init__(self, criteria_sets=None, node_tags={}):
        subgraph = Graph([(0, 1), (0, 2), (0, 3), (0, 4)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)


class NRingPattern(SubgraphPattern):
    def __init__(self, size, criteria_sets=None, node_tags={}, strong=False):
        self.size = size
        self.strong = strong
        subgraph = Graph([(i,(i+1)%size) for i in xrange(size)])
        SubgraphPattern.__init__(self, subgraph, criteria_sets, node_tags)

    def check_next_match(self, match, new_relations):
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
        if not SubgraphPattern.complete(self, match):
            return False
        if self.strong:
            # If the ring is not strong, return False
            if self.size%2 == 0:
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
