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
#


"""
This module contains tools for creating and analyzing graphs. For some historic
reasons, vertices are called nodes and edges are called pairs.

Some distinctive features include:
* Iterating over all shortest paths between two nodes (Dijkstra)
    http://en.wikipedia.org/wiki/Dijkstra's_algorithm
* Iterating over nodes or pairs using the Breadth First convention
    http://en.wikipedia.org/wiki/Breadth-first_search
* The all pairs shortest path matrix (Floyd-Warshall)
    http://en.wikipedia.org/wiki/Floyd-Warshall_algorithm
* Symmetry analysis of graphs (automorphisms)
    The Graph class can generate a list of permutations between nodes that
    map the graph onto itself. This can be used to generate (and test) all
    possible geometric symmetries in a molecule.
    http://en.wikipedia.org/wiki/Graph_automorphism
* Scanning a graph for patterns
    The GraphSearch is a generic class that can scan a graph for certain
    patterns, e.g. given subgraphs, strong rings, isomorphisms,
    automorphisms, ... The subgraph patter can deal with (multiple sets of)
    additional conditions that must be satisfied, such as "give me all dihedral
    angles where the central atoms are carbons" without duplicates.

The central class in this module is 'Graph'. It caches most of the analysis
results, which implies that the graph structure can not be changed once the
object is created. If you feel the need to do this, just construct a new graph
object.
"""


import copy, numpy


__all__ = [
    "OneToOneError", "OneToOne", "GraphError", "cached", "Graph",
    "Match", "PatternError", "Pattern",
    "CriteriaSet", "Anything", "CritOr", "CritAnd", "CritNodeString", "CritPairString",
    "SubgraphPatternError", "SubgraphPattern", "EqualPattern", "EgoMatch",
    "EgoPattern", "RingPattern", "GraphSearch",
]

class OneToOneError(Exception):
    pass


class OneToOne(object):
    """
    Implements a discrete bijection between source and destination elements.

    The implementation is based on a consistent set of forward and reverse
    relations, stored in dictionaries.
    """

    def __init__(self, pairs=[]):
        self.forward = {}
        self.reverse = {}
        self.add_relations(pairs)

    def __len__(self):
        return len(self.forward)

    def __str__(self):
        result = "|"
        for source, destination in self.forward.iteritems():
            result += " %s -> %s |" % (source, destination)
        return result

    def __copy__(self):
        class EmptyClass(object):
            pass
        result = EmptyClass()
        result.__class__ = self.__class__
        result.__dict__ = self.__dict__.copy()
        result.forward = self.forward.copy()
        result.reverse = self.reverse.copy()
        return result

    def __mul__(self, other):
        """Return the result of the 'after' operator."""
        result = OneToOne()
        for source, mid in other.forward.iteritems():
            destination = self.forward[mid]
            result.forward[source] = destination
            result.reverse[destination] = source
        return result

    def add_relation(self, source, destination):
        if self.in_sources(source):
            if self.forward[source] != destination:
                raise OneToOneError("Source is already in use. Destination does not match.")
            else:
                raise OneToOneError("Source-Destination relation already exists.")
        elif self.in_destinations(destination):
            raise OneToOneError("Destination is already in use. Source does not match.")
        else:
            self.forward[source] = destination
            self.reverse[destination] = source

    def add_relations(self, pairs):
        for source, destination in pairs:
            self.add_relation(source, destination)

    def get_destination(self, source):
        return self.forward[source]

    def get_source(self, destination):
        return self.reverse[destination]

    def in_destinations(self, destination):
        return destination in self.reverse

    def in_sources(self, source):
        return source in self.forward

    def inverse(self):
        """Returns the inverse bijection."""
        result = self.__class__()
        result.forward = copy.copy(self.reverse)
        result.reverse = copy.copy(self.forward)
        return result


class GraphError(Exception):
    pass


class cached(object):
    """A decorator that will turn a method into a caching descriptor.

    When attribute is requested for the firs time, the original method will
    be called and its return value is cached. Subsequent access to the attribute
    will just return the cached value.
    """
    def __init__(self, fn):
        self.fn = fn
        self.attribute_name = "_cache_%s" % fn.__name__
        self.__doc__ = fn.__doc__

    def __get__(self, instance, owner):
        value = getattr(instance, self.attribute_name, self)
        if value is self:
            #print "COMPUTING %s" % self.attribute_name
            value = self.fn(instance)
            setattr(instance, self.attribute_name, value)
        return value



class Graph(object):
    """The Graph class implements an undirected graph. All pairs have equal weight.

    The pairs attribute is the most elementary and is always available. All
    other attributes are optional and can be generated with their respective
    init_* methods.

    Graphs are meant to be immutable objects. Once created they are not
    supposed to be modified. If you want an adapted graph, create a new object.
    The main reason is that all graph analysis work is cached in this object.
    When the pairs change, the cache becomes invalid and should be erased. The
    latter is not supported as it is easier to create just a new graph object
    with the modified pairs.

    If you want to associate data with nodes or pairs, add custom dictionary or
    list attributes to the graph object. This is the way to work attributes that
    are not applicable to all nodes or pairs:

    node_property = {node1: blablabla, node2: blablabla, ...}
    pair_property = {frozenset([node1,node2]): blablabla, ..}

    If a certain property applies to all nodes or pairs, it is sometimes more
    practical to work with lists are numpy arrays that have the same ordering
    as the nodes or the pairs:

    node_property = numpy.array([6, 6, 1, 1, 1, 1], int) # atom numbers for ethene
    pair_property = numpy.array([2, 1, 1, 1, 1], int)    # bond orders for ethene

    """

    def __init__(self, pairs, num_nodes=None):
        """Initialize a Graph object.

        Arguments:
          pairs -- tuple(frozenset([node1, node2]), ...)
          num_nodes -- number of nodes

        node1 to nodeN must be integers from 0 to N-1. If nodes above the
        highest node value are not connected by pairs, use the num_nodes
        argument to tell what the total number of nodes is.

        If the pairs argument does not have the correct format, it will be
        converted.

        """
        self._pairs = []
        for pair in pairs:
            if len(pair) != 2:
                raise TypeError("The pairs must be a iterable with 2 elements")

                raise TypeError("The pairs must be frozen sets of exactly 2 elements each.")
            i,j = pair
            if i==j:
                raise ValueError("A pair must contain two different values.")
            if not (isinstance(i, int) and isinstance(j, int)):
                raise TypeError("The pairs must contain integers.")
            if i < 0 or j < 0:
                raise TypeError("The pairs must contain positive integers.")
            self._pairs.append(frozenset([i,j]))
        self._pairs = tuple(self._pairs)

        if len(self._pairs) == 0:
            self._num_nodes = 0
        else:
            self._num_nodes = max(max(a,b) for a,b in self._pairs)+1
        if num_nodes is not None:
            if not isinstance(num_nodes, int):
                raise TypeError("The optional argument num_nodes must be an integer when given.")
            if num_nodes < self._num_nodes:
                raise ValueError("num_nodes must be equal or larger to the number of nodes deduced from the pair list.")
            self._num_nodes = num_nodes

    pairs = property(lambda self: self._pairs)
    num_nodes = property(lambda self: self._num_nodes)
    num_pairs = property(lambda self: len(self._pairs))

    def __mul__(self, repeat):
        """Construct a graph that repeats this graph a number of times.

        Arguments:
          repeat -- The number of repetitions.
        """
        if not isinstance(repeat, int):
            raise TypeError("Can only multiply a graph with an integer")
        new_pairs = []
        for i in xrange(repeat):
            for node1, node2 in self.pairs:
                new_pairs.append(frozenset([node1+i*self.num_nodes,node2+i*self.num_nodes]))
        result = Graph(new_pairs, self.num_nodes*repeat)
        return result

    __rmul__ = __mul__

    def __str__(self):
        return " ".join("%i-%i" % tuple(sorted([i,j])) for i,j in self.pairs)

    # functions that should be implemented by derived classes

    def get_node_string(self, i):
        """Returns a string that fully characterizes the nature of the node.

        In case of an ordinary graph, all nodes are the same.
        """
        return ""

    def get_pair_string(self, i):
        """Returns a string that fully characterizes the nature of a pair.

        In case of an ordinary graph, all pairs are the same.
        """
        return ""

    # cached attributes:

    @cached
    def _work_nodes(self):  # for internal usage only
        return numpy.zeros(self._num_nodes, int)

    @cached
    def pair_index(self):
        """construct a map to look up the index of a pair"""
        return dict((pair,index) for index,pair in enumerate(self.pairs))

    @cached
    def neighbors(self):
        """Setup a dictionary with neighbors.

        The dictionary will have the following form:
          {nodeX: (nodeY1, nodeY2, ...), ...}
        This means that nodeX and nodeY1 have a pair etc. This also implies
        that the following elements are part of the dictionary:
          {nodeY1: (nodeX, ...), nodeY2: (nodeX, ...), ...}
        """
        neighbors = dict((node, set([])) for node in xrange(self.num_nodes))
        for a, b in self.pairs:
            neighbors[a].add(b)
            neighbors[b].add(a)
        return neighbors

    @cached
    def distances(self):
        """Construct the matrix with the all-pairs shortest path."""
        from molmod.ext import graphs_floyd_warshall
        distances = numpy.zeros((self.num_nodes, self.num_nodes), numpy.int32)
        #distances[:] = -1 # set all -1, which is just a very big integer
        #distances.ravel()[::len(distances)+1] = 0 # set diagonal to zero
        for i,j in self.pairs: # set pairs to one
            distances[i,j] = 1
            distances[j,i] = 1
        graphs_floyd_warshall(distances)
        return distances

    @cached
    def max_distance(self):
        """The maximum value in the distances matrix."""
        if self.distances.shape[0] == 0:
            return 0
        else:
            return self.distances.max()

    @cached
    def central_nodes(self):
        """Define the nodes that have the lowest maximum distance to any other node."""
        max_distances = self.distances.max(0)
        max_distances_min = max_distances[max_distances>0].min()
        return (max_distances==max_distances_min).nonzero()[0]

    @cached
    def central_node(self):
        """Define the node that has the lowest maximum distance to any other node.

        This definition does not lead to a unique solution. One arbitrary solution
        is selected.
        """
        return self.central_nodes[0]

    @cached
    def independent_nodes(self):
        """Generate lists of nodes are only interconnected within each list.

        This means that there is no path from a node in one list to a node
        in another list. In case of a molecular graph, this would yield the
        atoms that belong to individual molecules.
        """
        candidates = set(xrange(self.num_nodes))

        result = []
        while len(candidates) > 0:
            pivot = candidates.pop()
            group = [node for node, distance in self.iter_breadth_first(pivot)]
            candidates.difference_update(group)

            # this sort makes sure that the order of the nodes is respected
            group.sort()
            result.append(group)
        return result

    @cached
    def fingerprint(self):
        """A total graph fingerprint that is invariant under permutation if the
        node indexes. The chance that two different (molecular) graphs yield
        the same fingerprint is small but not zero. (See unit tests.)"""
        if self.num_nodes == 0:
            return numpy.zeros(20, numpy.ubyte)
        else:
            return sum(self.node_fingerprints)

    @cached
    def node_fingerprints(self):
        """A per-node fingerprint that is invariant under permutation if the
        node indexes. Atoms that are symmetrically equivalent will get the same
        finger print, e.g. the hydrogens in methane would get the same
        fingerprint."""
        return self.get_node_fingerprints(
            [self.get_node_string(i) for i in xrange(self.num_nodes)],
            [self.get_pair_string(i) for i in xrange(self.num_pairs)],
        )

    @cached
    def equivalent_nodes(self):
        """A dictionary with symmetrically equivalent nodes."""
        level1 = {}
        for i, row in enumerate(self.node_fingerprints):
            key = str(buffer(row))
            l = level1.get(key)
            if l is None:
                l = set([i])
                level1[key] = l
            else:
                l.add(i)
        level2 = {}
        for key, nodes in level1.iteritems():
            for node in nodes:
                level2[node] = nodes
        return level2

    @cached
    def symmetries(self):
        """Graph symmetries are permutations that map the graph onto itself."""

        symmetry_cycles = set([])
        symmetries = set([])
        for match in GraphSearch(EgoPattern())(self):
            match.cycles = match.get_closed_cycles()
            assert match.cycles not in symmetry_cycles, "Duplicates in EgoMatch"
            symmetry_cycles.add(match.cycles)
            symmetries.add(match)
        return symmetries

    @cached
    def symmetry_cycles(self):
        """Construct the cycle representations of the graph symmetries"""
        result = set([])
        for symmetry in self.symmetries:
            result.add(symmetry.cycles)
        return result

    @cached
    def canonical_order(self):
        """The nodes in a canonical or normalized order.

        This routine will return a list of nodes in an order that does not
        depend on the initial order, but only depends on the connectivity and
        the return values of the function self.get_node_string.

        Only the nodes that are involved in pairs will be included. The result
        can be given as first argument to self.get_subgraph, with reduce=True
        as second argument. This will return a complete canonical graph.

        The routine is designed not to use symmetry relations that are obtained
        with the GraphSearch routine. We also tried to create an ordering that
        feels like natural, i.e. starting in the center and pushing nodes with
        few equivalents to the front. If necessary, the nature of the nodes and
        their bonds to atoms closer to the center will also play a role, but
        only as a last resort.
        """
        # A) find an appropriate starting node.
        # Here we take a central node that has a minimal number of symmetrical
        # equivalents, 'the highest atom number', and the highest fingerprint.
        # Note that the symmetrical equivalents are computed from the node
        # fingerprints, i.e. without the GraphSearch.
        starting_node = max(
            (
                -len(self.equivalent_nodes[node]),
                self.get_node_string(node),
                str(buffer(self.node_fingerprints[node])),
                node
            ) for node in self.central_nodes
        )[-1]

        # B) sort all nodes based on
        #      1) distance from central node
        #      2) number of equivalent nodes
        #      3) node string, (higher atom numbers come first)
        #      4) fingerprint
        #      5) node index
        # The last field is only included to collect the result of the sort.
        # The fingerprint on itself would be sufficient, but the three first are
        # there to have a naturally appealing result.
        l = [
            [
                -distance,
                -len(self.equivalent_nodes[node]),
                self.get_node_string(node),
                str(buffer(self.node_fingerprints[node])),
                node
            ] for node, distance in self.iter_breadth_first(starting_node)
            if len(self.neighbors[node]) > 0
        ]
        l.sort(reverse=True)

        # C) The order of some nodes is still not completely set. e.g. consider
        # the case of allene. The four hydrogen atoms are equivalent, but one
        # can have two different orders: make geminiles consecutive or don't. It
        # is more trikcy than one would think at first sight. In the case of
        # allene, geminility could easily solve the problem. Consider a big flat
        # rotationally symmetric molecule (order 2). The first five shells are
        # order 4 and one would just give a random order to four segemnts in the
        # first shell. Only when one reaches the outer part that has order two,
        # it turns out that the arbitrary choices in the inner shell play a
        # role. So it does not help to look at relations with nodes at inner or
        # current shells only. One has to consider the whole picture.
        # (unit testing reveals troubles like these)

        # I need some sleep now. The code below checks for potential fuzz and
        # will raise an error if the ordering is not fully determined yet. One
        # day, I'll need this code more than I do now, and I'll fix things up.
        # I know how to do this, but I don't care enough right now.
        # -- Toon
        for i in xrange(1, len(l)):
            if l[i][:-1] == l[i-1][:-1]:
                raise NotImplementedError("Volunteers, please finish this routine!")

        # D) Return only the node indexes.
        return [record[-1] for record in l]

    # other usefull graph functions

    def iter_breadth_first(self, start=None, do_paths=False, do_duplicates=False):
        """Iterate over the nodes with the breadth first algorithm.

        See http://en.wikipedia.org/wiki/Breadth-first_search for more info.
        If not start node is given, the central node is taken.

        By default, the distance to the starting node is also computed. If the
        path to the starting node should be computed instead, set path to True.

        When duplicate is True, then nodes that can be reached trhough different
        paths of equal length, will be iterated twice. This tipically only makes
        sencse when path==True.
        """
        if start is None:
            start = self.central_node
        else:
            try:
                start = int(start)
            except ValueError:
                raise TypeError("First argument (start) must be an integer.")
            if start < 0 or start >= self.num_nodes:
                raise ValueError("start must be in the range [0,%i[" % self.num_nodes)
        from collections import deque
        work = self._work_nodes
        work[:] = -1
        work[start] = 0
        if do_paths:
            result = (start, 0, (start,))
        else:
            result = (start, 0)
        yield result
        todo = deque([result])
        while len(todo) > 0:
            if do_paths:
                parent, parent_length, parent_path = todo.popleft()
            else:
                parent, parent_length = todo.popleft()
            current_length = parent_length + 1
            for current in self.neighbors[parent]:
                visited = work[current]
                if visited == -1 or (do_duplicates and visited == current_length):
                    work[current] = current_length
                    if do_paths:
                        current_path = parent_path + (current,)
                        result = (current, current_length, current_path)
                    else:
                        result = (current, current_length)
                    #print "iter_breadth_first", result
                    yield result
                    todo.append(result)

    def iter_shortest_paths(self, a, b):
        """Iterate over all the shortest paths between node a and b."""
        max_len = None
        for node, length, path in self.iter_breadth_first(a, do_paths=True, do_duplicates=True):
            if max_len is not None:
                if length > max_len:
                    return
            if node == b:
                max_len = length
                yield path

    def iter_breadth_first_pairs(self, start=None):
        """Iterate over the pairs with the breadth first convention.

        We need this for the pattern matching algorithms, but a quick look at
        wikipedia did not result in a known and named algorithm.

        The pairs are yield one by one, together with the distance of the pair
        from the starting node and a flag that indicates whether the yielded
        pair connects two nodes that are at the same distance from the starting
        node. If that flag is False, the distance from the starting node to
        pair[0] is equal to the distance variable and the distance from pair[1]
        to the starting node is equal to distance+1. One item has the following
        format: ((i,j), distance, flag)
        """
        if start is None:
            start = self.central_node
        else:
            try:
                start = int(start)
            except ValueError:
                raise TypeError("First argument (start) must be an integer.")
            if start < 0 or start >= self.num_nodes:
                raise ValueError("start must be in the range [0,%i[" % self.num_nodes)
        from collections import deque
        work = self._work_nodes
        work[:] = -1
        work[start] = 0
        todo = deque([start])
        while len(todo) > 0:
            parent = todo.popleft()
            distance = work[parent]
            for current in self.neighbors[parent]:
                if work[current] == -1:
                    yield (parent,current), distance, False
                    work[current] = distance+1
                    todo.append(current)
                elif work[current] == distance and current>parent:
                    # second equation in elif avoids duplicates
                    yield (parent,current), distance, True
                elif work[current] == distance+1:
                    yield (parent,current), distance, False

    def get_subgraph(self, subnodes, normalize=False):
        """Constructs a subgraph of the current graph

        Arguments:
          subnodes -- The nodes that should be retained.
          normalize -- Whether or not the nodes should renumbered and reduced to
              the given set of subnodes. When True, also the pairs are sorted.
              It the end, this means that new order of the pairs does not depend
              on the original order, but only on the order of the argument
              subnodes.
              This option is False by default. When False, only pairs will be
              discarded, but the retained data remain unchanged. Also the
              parameter num_nodes is not affected.

        The returned graph will have an attribute old_pair_indexes that relates
        the positions of the new and the old pairs as follows:
          self.pairs[result.old_pair_indexes[i]] = result.pairs[i]
        In derived classes, the following should be supported:
          self.pair_property[result.old_pair_indexes[i]] = result.pair_property[i]

        When normalize==True, also the nodes are affected and the derived classes
        should make sure that the following works:
          self.node_property[result.old_node_indexes[i]] = result.node_property[i]
        The attribute old_node_indexes is only constructed when normalize==True.
        """
        if normalize:
            revorder = dict((j,i) for i,j in enumerate(subnodes))
            new_pairs = []
            old_pair_indexes = []
            for counter, (i,j) in enumerate(self.pairs):
                new_i = revorder.get(i)
                if new_i is None:
                    continue
                new_j = revorder.get(j)
                if new_j is None:
                    continue
                new_pairs.append((new_i, new_j))
                old_pair_indexes.append(counter)
            # sort the pairs
            order = range(len(new_pairs))
            order.sort( key=(lambda i: tuple(sorted(new_pairs[i]))) ) # argsort in pure python
            new_pairs = [new_pairs[i] for i in order]
            old_pair_indexes = [old_pair_indexes[i] for i in order]

            result = Graph(new_pairs)
            result.old_node_indexes = numpy.array(subnodes, dtype=int)
            #result.new_node_indexes = revorder
            result.old_pair_indexes = numpy.array(old_pair_indexes, dtype=int)
        else:
            subnodes = set(subnodes)
            old_pair_indexes = numpy.array([i for i, pair in enumerate(self.pairs) if pair.issubset(subnodes)], dtype=int)
            new_pairs = tuple(self.pairs[i] for i in old_pair_indexes)
            result = Graph(new_pairs, self.num_nodes)
            result.old_pair_indexes = old_pair_indexes
            # no need for old and new node_indexes because they remain the same.
        return result

    def get_node_fingerprints(self, node_strings, pair_strings, num_iter=None):
        import sha
        str2array = lambda x: numpy.frombuffer(x, numpy.ubyte)
        hashrow = lambda x: str2array(sha.new(x.data).digest())
        # initialization
        result = numpy.zeros((self.num_nodes, 20), numpy.ubyte)
        for i in xrange(self.num_nodes):
            result[i] = hashrow(str2array(node_strings[i]))
        for i in xrange(self.num_pairs):
            a,b = self.pairs[i]
            tmp = hashrow(str2array(pair_strings[i]))
            result[a] += tmp
            result[b] += tmp
        work = result.copy()
        # iterations
        if num_iter is None:
            num_iter = self.max_distance
        for i in xrange(num_iter):
            for a,b in self.pairs:
                work[a] += result[b]
                work[b] += result[a]
            #for a in xrange(self.num_nodes):
            #    for b in xrange(self.num_nodes):
            #        work[a] += hashrow(result[b]*self.distances[a,b])
            for a in xrange(self.num_nodes):
                result[a] = hashrow(work[a])
        return result

    def get_halfs(self, node1, node2):
        """Compute the two parts of the graph separated by the pair (node1, node2)

        If this is not possible (due to loops connecting both ends), a GraphError
        is raised.

        Returns the nodes in both halfs.
        """
        node1_new = set(self.neighbors[node1])
        if node2 not in node1_new:
            raise GraphError("Node1 and node2 must be connected.")
        node1_new.discard(node2)
        node1_part = set([node1])

        while len(node1_new) > 0:
            pivot = node1_new.pop()
            if pivot == node2:
                raise GraphError("The graph can not be separated in two halfs by disconnecting node1 and node2.")
            pivot_neighbors = set(self.neighbors[pivot])
            pivot_neighbors -= node1_part
            node1_new |= pivot_neighbors
            node1_part.add(pivot)

        # find node_b_part: easy, is just the rest
        node2_part = set(xrange(self.num_nodes)) - node1_part

        return node1_part, node2_part

    def get_halfs_double(self, node_a1, node_b1, node_a2, node_b2):
        """Compute the two parts separated by (node_a1, node_b1) and (node_a2, node_b2)

        Raise a GraphError when (node_a1, node_b1) and (node_a2, node_b2) do not
        separate the graph in two disconnected parts. The pairs must be
        neigbors. If not a GraphError is raised. The for nodes must not coincide
        or a GraphError is raised.

        Returns the nodes of the two halfs and the four 'hinge' nodes in the
        correct order, i.e. both node_a1 and node_a2 are in the first half and
        both node_b1 and node_b2 are in the second half.
        """
        if node_a1 not in self.neighbors[node_b1]:
            raise GraphError("Node_a1 must be a neighbor of node_b1.")
        if node_a2 not in self.neighbors[node_b2]:
            raise GraphError("Node_a2 must be a neighbor of node_b2.")

        # find node_a_part (and possibly switch node_a2, node_b2)
        node_a_new = set(self.neighbors[node_a1])
        node_a_new.discard(node_b1)
        if node_a1 == node_b2:
            # we now that we have to swap node_a2 and node_b2. The algo
            # below will fail otherwise in this 'exotic' case.
            node_a2, node_b2 = node_b2, node_a2
            #node_a_new.discard(node_a2) # in case there is overlap
        if node_a1 == node_a2: node_a_new.discard(node_b2) # in case there is overlap
        node_a_part = set([node_a1])

        touched = False # True if (the switched) node_a2 has been reached.
        while len(node_a_new) > 0:
            pivot = node_a_new.pop()
            if pivot == node_b1:
                raise GraphError("The graph can not be separated in two halfs. node_b1 reached by node_a1.")
            node_a_part.add(pivot)
            pivot_neighbors = set(self.neighbors[pivot]) # create a new set that we can modify
            pivot_neighbors -= node_a_part
            if pivot == node_a2 or pivot == node_b2:
                if pivot == node_b2:
                    if touched:
                        raise GraphError("The graph can not be separated in two halfs. node_b2 reached by node_a1.")
                    else:
                        # put them in the correct order
                        node_a2, node_b2 = node_b2, node_a2
                pivot_neighbors.discard(node_b2)
                touched = True
            node_a_new |= pivot_neighbors

        if node_a2 not in node_a_part:
            raise GraphError("The graph can not be separated in two halfs. node_a1 can not reach node_a2 trough node_a_part")

        # find node_b_part: easy, is just the rest ...
        #node_b_part = set(xrange(self.num_nodes)) - node_a_part

        # ... but we also want that there is a path in node_b_part from node_b1 to node_b2
        if node_b1 == node_b2:
            closed = True
        else:
            node_b_new = set(self.neighbors[node_b1])
            node_b_new.discard(node_a1)
            node_b_part = set([node_b1])

            closed = False
            while len(node_b_new) > 0:
                pivot = node_b_new.pop()
                if pivot == node_b2:
                    closed = True
                    break
                pivot_neighbors = set(self.neighbors[pivot])
                pivot_neighbors -= node_b_part
                node_b_new |= pivot_neighbors
                node_b_part.add(pivot)

        if not closed:
            raise GraphError("The graph can not be separated in two halfs. node_b1 can not reach node_b2 trough node_b_part")

        # finaly compute the real node_b_part, the former loop might break
        # early for efficiency.
        node_b_part = set(xrange(self.num_nodes)) - node_a_part

        # done!
        return node_a_part, node_b_part, (node_a1, node_b1, node_a2, node_b2)

    def full_match(self, other):
        """Given another graph with the same connectivity, find the mapping
        between node indexes in self and other.

        This also works on disconnected graphs. Derived classes should just
        implement get_node_string and get_pair_string to make this method
        aware of the different nature of certain nodes. In case molecules,
        this would make the algorithm sensitive to atom numbers etc.
        """
        graphs0 = [self.get_subgraph(group,normalize=True) for group in self.independent_nodes]
        # we need normalize subgraphs because these graphs are used as patterns.
        graphs1 = [other.get_subgraph(group) for group in other.independent_nodes]

        if len(graphs0) != len(graphs1):
            return

        matches = []

        for graph0 in graphs0:
            pattern = EqualPattern(graph0)
            found_match = False
            for i, graph1 in enumerate(graphs1):
                local_matches = list(GraphSearch(pattern)(graph1,one_match=True))
                if len(local_matches) == 1:
                    match = local_matches[0]
                    # we need to restore the relation between the normalize graph0
                    # and its original indexes
                    old_to_new = OneToOne(((j,i) for i,j in enumerate(graph0.old_node_indexes)))
                    matches.append(match * old_to_new)
                    del graphs1[i]
                    found_match = True
                    break
            if not found_match:
                return

        result = OneToOne()
        for match in matches:
            result.add_relations(match.forward.iteritems())
        return result


# Pattern matching


class Match(OneToOne):
    def __init__(self, node0, node1):
        OneToOne.__init__(self, [(node0,node1)])
        self.previous_ends1 = set([node1])

    def get_new_pairs(self, graph):
        result = []
        #print "Match.get_new_pairs self.previous_ends1", self.previous_ends1
        for node in self.previous_ends1:
            for neighbor in graph.neighbors[node]:
                if neighbor not in self.reverse:
                    result.append((node, neighbor))
        return result

    def copy_with_new_relations(self, new_relations):
        result = self.__copy__()
        result.add_relations(new_relations.iteritems())
        result.previous_ends1 = set(new_relations.itervalues())
        return result


class PatternError(Exception):
    pass


class Pattern(object):
    sub = True # This means that matching nodes must not have equal number of neighbors

    """Base class for a pattern in a graph.

    Note the following conventions:
      * A pattern can always be represented by a graph (or a set of graphs) and
        some additional conditions. This graph is the so called 'PATTERN GRAPH'.
        For technical reasons, this pattern graph is not always constructed
        explicitly. Variables related to this graph often get suffix '0'. Note
        that a pattern graph is always fully connected.
      * The graph in which we search for the pattern, is called the 'SUBJECT
        GRAPH'. Variables related to this graph often get suffix '1'.
    """
    MatchClass = Match

    def init_graph(self, graph, one_match):
        """Checks whether this pattern makes sense for the given subject graph."""
        self.graph = graph
        self.one_match = one_match

    def iter_initial_relations(self):
        """Yields the initial relations (source,destination) to start with.

        The function iterates of single relations between a pattern node and a
        subject node. In practice it is sufficient to select one node in the
        pattern and then associate it with each (reasonable) node in the subject
        graph.
        """
        raise NotImplementedError

    def get_new_pairs(self, level):
        raise NotImplementedError

    def check_symmetry(self, new_relations, current_match, next_match):
        "Check wether the new_relations correspond the reference case of all possible symetric equivalents"
        return True

    def compare(self, node0, node1):
        """
        Test if node0 and node1 can be equal. False positives are allowed, but
        the less false positives, the more efficient the GraphSearch will be.
        """
        return True

    def check_next_match(self, match, new_relations):
        """Does this match object make sense for the current pattern.

        Return False if some symmetry or other considerations are not satisfied.
        The checks in this function are typically only possible by considering
        the whole instead of looking just at a few nodes/pairs/relations.
        """
        return True

    def complete(self, match):
        """Returns True if not more additional relations are required."""
        return True

    def iter_final_matches(self, match):
        yield match


class CriteriaSet(object):
    """A set of criteria that can be associated with a subgraph pattern."""
    def __init__(self, thing_criteria=None, relation_criteria=None, global_criteria=None, **kwargs):
        if thing_criteria is None:
            self.thing_criteria = {}
        else:
            self.thing_criteria = thing_criteria
        if relation_criteria is None:
            self.relation_criteria = {}
        else:
            self.relation_criteria = relation_criteria
        if global_criteria is None:
            self.global_criteria = {}
        else:
            self.global_criteria = global_criteria
        self.info = kwargs

    def test_match(self, match, subgraph, graph):
        for node0, c in self.thing_criteria.iteritems():
            node1 = match.forward[node0]
            if not c(node1, graph): return False
        for pair0_index, c in self.relation_criteria.iteritems():
            node0a, node0b = subgraph.pairs[pair0_index]
            pair1_index = graph.pair_index[frozenset([
                match.forward[node0a],
                match.forward[node0b],
            ])]
            if not c(pair1_index, graph): return False
        for c in self.global_criteria:
            if not c(match, graph): return False
        return True

# few basic example criteria

class Anything(object):
    def __call__(self, index, graph):
        return True


class CritOr(object):
    def __init__(self, *criteria):
        self.criteria = criteria

    def __call__(self, index, graph):
        for c in self.criteria:
            if c(index, graph):
                return True
        return False


class CritAnd(object):
    def __init__(self, *criteria):
        self.criteria = criteria

    def __call__(self, index, graph):
        for c in self.criteria:
            if not c(index, graph):
                return False
        return True


class CritNodeString(object):
    def __init__(self, reference):
        self.reference = reference

    def __call__(self, index, graph):
        s0 = self.reference.get_node_string(index)
        s1 = graph.get_node_string(index)
        return s1=="" or s2=="" or s1==s2 # an aspecific node acts as a wildcard


class CritPairString(object):
    def __init__(self, reference):
        self.reference = reference

    def __call__(self, index, graph):
        s0 = self.reference.get_pair_string(index)
        s1 = graph.get_pair_string(index)
        return s1=="" or s2=="" or s1==s2 # an aspecific node acts as a wildcard


# pattern and match stuff

class SubgraphPatternError(PatternError):
    pass


class SubgraphPattern(Pattern):
    def __init__(self, subgraph, criteria_sets=None, node_tags={}):
        """Initialise a subgraph pattern.

        Arguments:
          subgraph -- the pattern that has to be found in the subject graph.
          criteria_sets -- Criteria sets associate additional conditions with
              nodes and pairs, and can also introduce global match conditions.
          node_tags -- Node tags can reduce the symmetry of the subgraph
              pattern. An example case where this is usefull: Consider atoms
              0,1,2 that are bonded in this order. We want to compute the
              distance from atom 2 to the line (0,1). In this case the
              matches (0->a, 1->b, 2->c) and (0->c, 1->b, 2->a) correspond to
              different internal coordinates. We want the graph search to return
              the two solutions. In order to do this, set node_tags={0:0,1:0,2:1}.
              This means that node 0 and 1 are equivalent, but that node 2 has a
              different nature. In the case of a bending angle, only one match
              like (0->a, 1->b, 2->c) is sufficient and we do not want to reduce
              the symmetry of the subgraph. In this case, one should not use
              node_tags at all.
        """
        self.criteria_sets = criteria_sets
        self.node_tags = node_tags
        # get the essential information from the subgraph:
        self.set_subgraph(subgraph)
        Pattern.__init__(self)

    def set_subgraph(self, subgraph):
        if subgraph is not None and len(subgraph.independent_nodes) != 1:
            raise SubgraphPatternError("A subgraph pattern must not be a disconnected graph.")
        self.subgraph = subgraph
        # A) the levels for the incremental pattern matching
        self.level_pairs = {}
        self.level_constraints = {}
        if subgraph is not None:
            for pair, distance, constraint in self.subgraph.iter_breadth_first_pairs():
                if constraint:
                    l = self.level_constraints.setdefault(distance-1, [])
                else:
                    l = self.level_pairs.setdefault(distance, [])
                l.append(pair)
        #print "level_pairs", self.level_pairs
        #print "level_constraints", self.level_constraints
        # B) The comparisons the should me checked when one wants to avoid
        # symmetrically duplicate pattern matches
        self.duplicate_checks = set([])
        if not (subgraph is None or self.criteria_sets is None):
            for cycles in subgraph.symmetry_cycles:
                if len(cycles) > 0:
                    self.duplicate_checks.add((cycles[0][0], cycles[0][1]))


    def iter_initial_relations(self):
        node0 = self.subgraph.central_node
        for node1 in xrange(self.graph.num_nodes):
            if self.compare(node0, node1):
                yield node0, node1

    def get_new_pairs(self, level):
        return self.level_pairs.get(level, []), self.level_constraints.get(level, [])

    def check_next_match(self, match, new_relations):
        # only returns true for ecaxtly one set of new_relations from all the
        # ones that are symmetrically equivalent
        if not (self.criteria_sets is None or self.one_match):
            for check in self.duplicate_checks:
                node_a = new_relations.get(check[0])
                node_b = new_relations.get(check[1])
                if node_a is None and node_b is None:
                    continue # if this pair is completely absent in the new
                    # relations, it is either completely in the match or it
                    # is to be matched. So it is either already checked for
                    # symmetry duplicates, or it will be check in future.
                if node_a is None:
                    # maybe node_a is in the match and node_b is the only one
                    # in the new relations. try to get node_a from the match.
                    node_a = match.forward.get(check[0])
                    if node_a is None:
                        # ok, node_a is to be found, don't care about it right
                        # now. it will be checked in future calls.
                        continue
                elif node_b is None:
                    # maybe node_b is in the match and node_a is the only one
                    # in the new relations. try to get node_b from the match.
                    node_b = match.forward.get(check[1])
                    if node_b is None:
                        # ok, node_b is to be found, don't care about it right
                        # now. it will be checked in future calls.
                        continue
                if node_a > node_b:
                    # Why does this work? The answer is not so easy to explain,
                    # and certainly not easy to find. if node_a > node_b, it
                    # means that there is a symmetry operation that leads to
                    # an equivalent match where node_b < node_a. The latter
                    # match is prefered for as much pairs (node_a,node_b) as
                    # possible without rejecting all possible matches. The real
                    # difficulty is to construct a proper list of
                    # (node_a,node_b) pairs that will reject all but one matches.
                    # I conjecture that this list contains all the first two
                    # nodes from each normalized symmetry cycle of the sub graph.
                    # We need a math guy to do the proof. -- Toon
                    return False
            return True
        return True

    def complete(self, match):
        return len(match) == self.subgraph.num_nodes

    def iter_final_matches(self, canonical_match):
        if self.criteria_sets is None or self.one_match:
            yield canonical_match
        else:
            for criteria_set in self.criteria_sets:
                satisfied_match_tags = set([])
                for symmetry in self.subgraph.symmetries:
                    final_match = canonical_match * symmetry
                    #print final_match
                    if criteria_set.test_match(final_match, self.subgraph, self.graph):
                        match_tags = tuple(
                            self.node_tags.get(symmetry.forward[node0])
                            for node0
                            in xrange(self.subgraph.num_nodes)
                        )
                        if match_tags not in satisfied_match_tags:
                            final_match.__dict__.update(criteria_set.info)
                            yield final_match
                            satisfied_match_tags.add(match_tags)


class EqualPattern(SubgraphPattern):
    sub = False

    def __init__(self, subgraph):
        # Don't allow criteria sets and node_tags
        SubgraphPattern.__init__(self, subgraph)

    def iter_initial_relations(self):
        if self.subgraph.num_pairs != self.graph.num_pairs:
            return # don't even try
        node0 = self.subgraph.central_node
        for node1 in self.graph.central_nodes:
            if self.compare(node0, node1):
                yield node0, node1

    def compare(self, node0, node1):
        return (self.subgraph.node_fingerprints[node0] == self.graph.node_fingerprints[node1]).all()


class EgoMatch(Match):
    def get_closed(self):
        """Return wether this permutation is closed."""
        for source in self.forward:
            if source not in self.reverse:
                return False
        return True

    def get_closed_cycles(self):
        """Return the closed cycles corresponding to this permutation.

        The cycle will be normalized to facilitate the elimination of
        duplicates. The following is guaranteed:
          1) If this permutation is represented by disconnected cycles, the
             cycles will be sorted by the lowest index they contain.
          2) Each cycle starts with its lowest index. (unique starting point)
          3) Singletons are discarded. (because they are boring)

        The following normalization step is commented out. It does NOT apply.
          4) The second element in a cycle is lower than the last element in a
             cycle. (unique direction). This manipulation is only caried out
             in the end when all cycles are known. It is applied on the first
             cycle and all subsequent cycles undergo the same transformation so
             that the interpretation of the cycle notation is not altered.
        """
        # A) construct all the cycles
        closed_cycles = []
        todo = set(self.forward.keys())
        current_node = None
        while len(todo) > 0:
            if current_node == None:
                current_node = todo.pop()
                current_cycle = []
            else:
                todo.discard(current_node)
            current_cycle.append(current_node)
            next_node = self.get_destination(current_node)
            if next_node == current_cycle[0]:
                if len(current_cycle) > 1:
                    # bring the lowest element in front
                    pivot = numpy.argmin(current_cycle)
                    current_cycle = current_cycle[pivot:] + current_cycle[:pivot]
                    closed_cycles.append(current_cycle)
                current_node = None
            else:
                current_node = next_node
        # B) normalize the cycle representation
        closed_cycles.sort() # a normal sort is sufficient because only the
                             # first item of each cycle is considered

        #if len(closed_cycles) > 0 and closed_cycles[0][1] > closed_cycles[0][-1]:
        #    # swap all cycles and bring the lowest in front
        #    for cycle in closed_cycles:
        #        cycle.reverse()
        #        cycle.insert(0,cycle.pop(-1))

        # transform the structure into a tuple of tuples
        closed_cycles = tuple(tuple(cycle) for cycle in closed_cycles)
        return closed_cycles



class EgoPattern(EqualPattern):
    MatchClass = EgoMatch

    def __init__(self):
        EqualPattern.__init__(self, None)

    def init_graph(self, graph, one_match):
        self.set_subgraph(graph)
        EqualPattern.init_graph(self, graph, one_match)


class RingPattern(Pattern):
    def __init__(self, max_size):
        if max_size < 3:
            raise ValueError("Ring sizes must be at least 3.")
        self.max_size = max_size

    def init_graph(self, graph, one_match):
        Pattern.init_graph(self, graph, one_match)

    def iter_initial_relations(self):
        node0 = 0
        for node1 in xrange(self.graph.num_nodes):
            yield node0, node1

    def get_new_pairs(self, level):
        if level == 0:
            pairs0 = [(0,1),(0,2)]
        elif level >= (self.max_size-1)/2:
            pairs0 = []
        else:
            l2 = level*2
            pairs0 = [(l2-1,l2+1),(l2,l2+2)]
        return pairs0, []

    def check_next_match(self, match, new_relations):
        # avoid duplicate rings (order of traversal)
        if len(match) == 3:
            if match.forward[1] < match.forward[2]:
                #print "RingPattern.check_next_match: duplicate order"
                return False
        # avoid duplicate rings (starting point)
        for node1 in new_relations.itervalues():
            if node1 < match.forward[0]:
                #print "RingPattern.check_next_match: duplicate start"
                return False
        # can this ever become a strong ring?
        for node1 in new_relations.itervalues():
            paths = list(self.graph.iter_shortest_paths(node1, match.forward[0]))
            if len(paths) != 1:
                #print "RingPattern.check_next_match: not strong 1"
                return False
            if len(paths[0]) != (len(match)+1)/2:
                #print "RingPattern.check_next_match: not strong 2"
                return False
        #print "RingPattern.check_next_match: no remarks"
        return True

    def complete(self, match):
        size = len(match)
        # check whether we have an odd strong ring
        if match.forward[size-1] in self.graph.neighbors[match.forward[size-2]]:
            # we have an odd closed cycle. check if this is a strong ring
            order = range(0,size,2) + range(1, size-1, 2)[::-1]
            ok = True
            for i in xrange(len(order)/2):
                if len(list(self.graph.iter_shortest_paths(match.forward[order[i]], match.forward[order[(i+size/2)%size]]))) > 1:
                    ok = False
                    break
                if len(list(self.graph.iter_shortest_paths(match.forward[order[i]], match.forward[order[(i+size/2+1)%size]]))) > 1:
                    ok = False
                    break
            if ok:
                match.ring_nodes = tuple(match.forward[i] for i in order)
                #print "RingPattern.complete: found odd ring"
                return True
            #print "RingPattern.complete: no odd ring"
        # check whether we have an even strong ring
        paths = list(self.graph.iter_shortest_paths(match.forward[size-1], match.forward[size-2]))
        #print "RingPattern.complete: even paths", paths
        if (size > 3 and len(paths) == 1 and len(paths[0]) == 3) or \
           (size == 3 and len(paths) == 2 and len(paths[0]) == 3):
            path = paths[0]
            if size == 3 and path[1] == match.forward[0]:
                path = paths[1]
            # we have an even closed cycle. check if this is a strong ring
            match.add_relation(size, path[1])
            size += 1
            order = range(0,size,2) + range(size-1, 0, -2)
            ok = True
            for i in xrange(len(order)/2):
                if len(list(self.graph.iter_shortest_paths(match.forward[order[i]], match.forward[order[(i+size/2)%size]]))) != 2:
                    ok = False
                    break
            if not ok:
                node1 = match.forward[size-1]
                del match.forward[size-1]
                del match.reverse[node1]
                #print "RingPattern.complete: no even ring"
            else:
                match.ring_nodes = tuple(match.forward[i] for i in order)
                #print "RingPattern.complete: found even ring"
            return ok
        #print "RingPattern.complete: not at all"
        return False


class GraphSearch(object):
    def __init__(self, pattern, debug=False):
        self.pattern = pattern
        self.debug = debug

    def __call__(self, graph, one_match=False):
        """Yields all match of self.pattern in the given graph.

        Arguments:
            graph  --  The graph in which the matches according to
                       self.matchdefinition have to be found.
            one_match --  If True, only one match will be returned. This
                          allows certain optimizations.
        """
        self.pattern.init_graph(graph, one_match)
        # Matches are grown iteratively.
        for node0, node1 in self.pattern.iter_initial_relations():
            init_match = self.pattern.MatchClass(node0, node1)
            # init_match cotains only one source -> dest relation. starting from
            # this initial match, the function iter_matches extends the match
            # in all possible ways and yields the completed matches
            for canonical_match in self._iter_matches(init_match, graph):
                # Some patterns my exclude symmetrically equivalent matches as
                # to aviod dupplicates. with such a 'canonical' solution,
                # the pattern is allowed to generate just those symmatrical
                # duplicates of interest.
                for final_match in self.pattern.iter_final_matches(canonical_match):
                    self.print_debug("final_match: %s" % final_match)
                    yield final_match
                    if one_match: return

    def print_debug(self, text, indent=0):
        if self.debug:
            if indent > 0:
                print " "*self.debug, text
            self.debug += indent
            if indent <= 0:
                print " "*self.debug, text

    def _iter_candidate_groups(self, init_match, pairs0, pairs1):
        # collect all end nodes0 and end nodes1 that belong to the same group.
        sources = {}
        for start_node0, end_node0 in pairs0:
            l = sources.setdefault(start_node0, [])
            l.append(end_node0)
        dests = {}
        for start_node1, end_node1 in pairs1:
            start_node0 = init_match.reverse[start_node1]
            l = dests.setdefault(start_node0, [])
            l.append(end_node1)
        for start_node0, end_nodes0 in sources.iteritems():
            end_nodes1 = dests.get(start_node0, [])
            yield end_nodes0, end_nodes1


    def _iter_new_relations(self, init_match, graph, pairs0, constraints0, pairs1):
        # Count the number of unique pairs0[i][1] values. This is also be
        # the number of new relations.
        num_new_relations = len(set(j for i,j in pairs0))

        def combine_small(relations, num):
            if len(relations) == 0:
                return
            for i, pivot in enumerate(relations):
                if num == 1:
                    yield (pivot,)
                else:
                    compatible_relations = list(
                        item for item in relations[:i]
                        if pivot[0]!=item[0] and pivot[1]!=item[1]
                    )
                    for tail in combine_small(compatible_relations, num-1):
                        yield (pivot,) + tail

        # generate candidate relations
        candidate_relations = []
        for end_nodes0, end_nodes1 in self._iter_candidate_groups(init_match, pairs0, pairs1):
            if len(end_nodes0) > len(end_nodes1):
                return # this can never work, the subject graph is 'too small'
            elif not self.pattern.sub and len(end_nodes0) != len(end_nodes1):
                return # an exact match is sought, this can never work
            l = []
            for end_node0 in end_nodes0:
                for end_node1 in end_nodes1:
                    if self.pattern.compare(end_node0,end_node1):
                        l.append((end_node0,end_node1))
            # len(end_nodes0) = the total number of relations that must be made in this group
            if len(l) > 0:
                # turn l into a list of sets of internally compatible candidate
                # relations in this group
                l = list(combine_small(l, len(end_nodes0)))
                candidate_relations.append(l)
        if len(candidate_relations) == 0:
            return
        self.print_debug("candidate_relations: %s" % candidate_relations)

        def combine_big(pos=0):
            # pos is an index in candidate_relations
            crs = candidate_relations[pos]
            if pos == len(candidate_relations)-1:
                for relations in crs:
                    yield relations
            else:
                for tail in combine_big(pos+1):
                    for relations in crs:
                        yield relations + tail

        # final loop
        for new_relations in combine_big():
            new_relations = set(new_relations)
            self.print_debug("new_relations: %s" % (new_relations,))
            # check the total number of new relations
            if len(new_relations) != num_new_relations:
                continue
            # check sanity of relations
            forward = dict(new_relations)
            if len(forward) != num_new_relations:
                continue
            reverse = dict((j,i) for i,j in new_relations)
            if len(reverse) != num_new_relations:
                continue
            # check the constraints
            for a0,b0 in constraints0:
                if forward[a0] not in graph.neighbors[forward[b0]]:
                    forward = None
                    break
            if forward is None:
                continue
            yield forward

    def _iter_matches(self, input_match, graph, level=0):
        self.print_debug("ENTERING _ITER_MATCHES", 1)
        self.print_debug("input_match: %s" % input_match)
        # A) collect the new pairs to extend the match.
        # Note that the pairs are ordered. pair[0] is always in the match.
        # pair[1] is never in the match. The constraints contain information
        # about the end points of pairs0. It is a list of two-tupples where
        # (a,b) means that a and b must be bonded.
        pairs0, constraints0 = self.pattern.get_new_pairs(level)
        pairs1 = input_match.get_new_pairs(graph)
        self.print_debug("pairs0: %s" % pairs0)
        self.print_debug("constraints0: %s" % constraints0)
        self.print_debug("pairs1: %s" % pairs1)

        # B) iterate over the sets of new relations: [(node0[i],node1[j]), ...]
        # that contain all endpoints of pairs0, that satisfy the constraints0
        # and where (node0[i],node1[j]) only occurs if these are end points
        # of a pair0 and pair1 whose starting points are already in init_match.
        # These conditions are implemented in an iterator as to separate concerns.
        # This iterator also calls the routines that check whether node1[j] also
        # satisfies additional conditions inherent node node0[i].
        for new_relations in self._iter_new_relations(input_match, graph, pairs0, constraints0, pairs1):
            # for each set of new_relations, construct a next_match and recurse
            next_match = input_match.copy_with_new_relations(new_relations)
            if not self.pattern.check_next_match(next_match, new_relations):
                continue
            if self.pattern.complete(next_match):
                yield next_match
            else:
                for match in self._iter_matches(next_match, graph, level+1):
                    yield match
        self.print_debug("LEAVING _ITER_MATCHES", -1)


