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
#
"""
This module contains tools for creating and analyzing graphs. The main
features are:

A) symmetry analysis of graphs
The Graph class will generate a list of permutations between nodes (during
it's initialization) that map the graph onto it'self.

B) scanning a graph for subgraphs
the yield_matches method iterates over all subgraphs congruent to the given
graph that are completely contained within self. This method requires the
symmetries of the graph to be known, in order to avoid duplicates.
"""


from clusters import ClusterFactory

import copy, numpy


__all__ = ["OneToOne", "Permutation", "Graph", "MatchGenerator"]


class OneToOneError(Exception):
    pass


class OneToOne(object):
    """
    Implements a discrete bijection between source and destination elements.

    The implementation is based on a consistent set of forward and backward
    relations, stored in dictionaries.
    """

    def __init__(self, pairs=[]):
        self.forward = {}
        self.backward = {}
        for source, destination in pairs:
            self.add_relation(source, destination)

    def __len__(self):
        return len(self.forward)

    def __str__(self):
        result = "|"
        for source, destination in self.forward.iteritems():
            result += " %s -> %s |" % (source, destination)
        return result

    def __copy__(self):
        result = self.__class__()
        result.forward = copy.copy(self.forward)
        result.backward = copy.copy(self.backward)
        return result

    def __mul__(self, other):
        """Return the result of the 'after' operator."""
        result = self.__class__()
        for source, mid in other.forward.iteritems():
            destination = self.forward[mid]
            result.forward[source] = destination
            result.backward[destination] = source
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
            self.backward[destination] = source

    def get_destination(self, source):
        return self.forward[source]

    def get_source(self, destination):
        return self.backward[destination]

    def in_destinations(self, destination):
        return destination in self.backward

    def in_sources(self, source):
        return source in self.forward

    def inverse(self):
        """Returns the inverse bijection."""
        result = self.__class__()
        result.forward = copy.copy(self.backward)
        result.backward = copy.copy(self.forward)
        return result


class Permutation(OneToOne):
    """
    A special type of bijection where both source and destination elements
    are part of the same set. A permutation is closed when each destination
    is also a source.
    """

    def get_closed(self):
        """Return wether this permutation is closed."""
        for source in self.forward:
            if source not in self.backward:
                return False
        return True

    def get_closed_cycles(self):
        """Return the closed cycles that form this permutation."""
        closed_cycles = []
        sources = self.forward.keys()
        current_source = None
        current_cycle = []
        while len(sources) > 0:
            if current_source == None:
                current_source = sources[0]
                current_cycle = []
            current_cycle.append(current_source)
            sources.remove(current_source)
            current_source = self.get_destination(current_source)
            if current_source == current_cycle[0]:
                if len(current_cycle) > 1:
                    closed_cycles.append(tuple(current_cycle))
                current_source = None
        return tuple(closed_cycles)


def all_relations(list1, list2, worth_trying=None):
    """
    Yields all possible relation sets between elements from list1 and list2.

    This is a helper function for the Graph class.

    Arguments:
    list1 -- list of source items
    list2 -- list of destination items
    worth_trying -- a function that takes a relation (tuple) as argument and
                    returns False if the relation is not usefull.
    """
    if len(list1) == 0 or len(list2) == 0:
        yield []
    else:
        for index2 in xrange(len(list2)):
            pair = (list1[0], list2[index2])
            if worth_trying == None or worth_trying(pair):
                for relation_set in all_relations(list1[1:], list2[:index2] + list2[index2+1:], worth_trying):
                    yield [pair] + relation_set


class Graph(object):
    """
    A Graph object contains two typical pythonic (not oriented) graph
    representations: pairs and neigbors.
    """

    def __init__(self, pairs, ordered_nodes=None):
        """
        Initialize a Graph object.

        Arguments:
        ordered_nodes -- [node1, node2, ...]
        pairs -- [frozenset([node1, node2]), ...]

        """
        self.nodes = ordered_nodes
        self.pairs = pairs
        self.init_index()
        self.init_neighbors()
        self.init_trees_and_shells()
        self.init_distances()
        self.init_central_node()

    def init_index(self):
        tmp = set([])
        for a, b in self.pairs:
            tmp.add(a)
            tmp.add(b)
        if self.nodes is None:
            self.nodes = list(tmp)
            self.nodes.sort()
        else:
            assert tmp.issubset(self.nodes)
            assert len(tmp) == len(self.nodes)
        self.index = dict((node, index) for index, node in enumerate(self.nodes))

    def init_neighbors(self):
        self.neighbors = dict((node, set([])) for node in self.nodes)
        for a, b in self.pairs:
            self.neighbors[a].add(b)
            self.neighbors[b].add(a)

    def init_trees_and_shells(self):
        self.trees = {}
        self.shells = {}
        self.shell_sizes = {}
        for central_node in self.nodes:
            tree = dict((node, set([])) for node in self.nodes)
            shells = []
            shell_sizes = []
            excludes = set([central_node])
            previous_shell = set([central_node])
            while len(previous_shell) > 0:
                shells.append(previous_shell)
                shell_sizes.append(len(previous_shell))
                new_shell = set([])
                for previous_node in previous_shell:
                    new_neighbors = set([
                        neighbor
                        for neighbor in self.neighbors[previous_node]
                        if neighbor not in excludes
                    ])
                    for new_neighbor in new_neighbors:
                        tree[new_neighbor].add(previous_node)
                    new_shell.update(new_neighbors)
                excludes.update(new_shell)
                previous_shell = new_shell
            self.shells[central_node] = shells
            self.shell_sizes[central_node] = numpy.array(shell_sizes, int)
            self.trees[central_node] = tree

    def init_distances(self):
        self.distances = numpy.zeros((len(self.nodes), len(self.nodes)), int)
        for node, shells in self.shells.iteritems():
            for distance, shell in enumerate(shells):
                for shell_node in shell:
                    self.distances[self.index[node], self.index[shell_node]] = distance

    def init_central_node(self):
        self.central_node = self.nodes[self.distances.max(0).argmin()]

    def can_overlap(self, central_node, node_a, node_b):
        tree = self.trees[central_node]
        # sanity checks
        if node_a == central_node or node_b == central_node:
            raise GraphError("can_overlap only makes sense of node_a != central_node and node_b != central_node")
        distance_a = self.distances[self.index[node_a], self.index[central_node]]
        distance_b = self.distances[self.index[node_b], self.index[central_node]]
        if distance_a == 0 or distance_b == 0:
            raise GraphError("can_overlap is only applicable when node_a and node_b are connected to central_node")
        if distance_a != distance_b:
            raise GraphError("can_overlap only works if distance(central_node, node_a) == distance(central_node, node_b)")
        if node_a == node_b:
            return True
        # the real algorithm
        return self._unsafe_can_overlap(tree, node_a, node_b, distance_a)

    def _unsafe_can_overlap(self, tree, node_a, node_b, distance):
        print "_"*distance, "_unsafe_can_overlap", distance
        if distance > 0:
            down_a = tree[node_a]
            down_b = tree[node_b]
            if len(down_a.intersection(down_b)) > 0:
                return True
            else:
                for new_node_a in down_a:
                    for new_node_b in down_b:
                        if self._unsafe_can_overlap(tree, node_a, node_b, distance-1):
                            return True
        else:
            return False


class MatchFilterError(Exception):
    pass


class MatchFilter(object):
    def init_graphs(self, subgraph, graph):
        "Checks initialy whether it makes sense to match both graphs."
        self.subgraph = subgraph
        self.graph = graph

    def new_pools(self, node0, node1, partial_match, shell_index):
        "Creates pool0 for subgraph and pool1 for graph that contain the nodes for the new relations."
        # The default behavior is to allow all nodes that have not been used yet.
        pool0 = set(self.subgraph.nodes).difference(partial_match.forward.keys())
        pool1 = set(self.graph.nodes).difference(partial_match.backward.keys())
        return pool0, pool1

    def useful_relation_group(self, relation_group):
        "Checks if the relation_group is worth looking at"
        return not (len(relation_group[0]) == 0 or len(relation_group[1]) == 0)

    def valid_relation_group(self, relation_group):
        "Checks if the relation_group is valid"
        return True

    def check_symmetry(self, new_relations):
        "Check wether the new_relations correspond the reference case of all possible symetric equivalents"
        return True

    def compare(self, node0, node1):
        """
        Test if node0 and node1 can be equal. False positives are allowed, but
        the less false positives, the more efficient the MatchGenerator will be.
        """
        return True

    def complete(self, node0, node1, match, shell_index):
        return True


class DistanceMatchFilter(MatchFilter):
    def new_pools(self, node0, node1, partial_match, shell_index):
        # It is often a good approximation and also in many cases a 'physical'
        # requirement that distance(node0, new0) == distance(node1, new1) for
        # every relation (new0, new1) in the match.
        if len(self.subgraph.shells[node0]) <= shell_index:
            pool0 = set([])
        else:
            pool0 = self.subgraph.shells[node0][shell_index]
        if len(self.graph.shells[node1]) <= shell_index:
            pool1 = set([])
        else:
            pool1 = self.graph.shells[node1][shell_index]
        return pool0, pool1



class ExactMatchFilter(DistanceMatchFilter):
    def init_graphs(self, subgraph, graph):
        if len(subgraph.nodes) != len(graph.nodes):
            raise MatchFilterError("It does not make sense to find an exact match between two graphs if the number of nodes is different")
        if len(subgraph.pairs) != len(graph.pairs):
            raise MatchFilterError("It does not make sense to find an exact match between two graphs if the number of relations is different")
        DistanceMatchFilter.init_graphs(self, subgraph, graph)

    def valid_relation_group(self, relation_group):
        return len(relation_group[0]) == len(relation_group[1])

    def check_symmetry(self, new_relations):
        return True

    def compare(self, node0, node1):
        sizes0 = self.subgraph.shell_sizes[node0]
        sizes1 = self.graph.shell_sizes[node1]
        return (
            (sizes0.shape == sizes1.shape) and
            (self.subgraph.shell_sizes[node0] == self.graph.shell_sizes[node1]).all()
        )

    def complete(self, node0, node1, match, shell_index):
        return len(self.subgraph.nodes) == len(match.forward)


class MatchGenerator(object):
    def __init__(self, match_filter, subgraph, graph, debug=False):
        self.match_filter = match_filter
        self.subgraph = subgraph
        self.graph = graph
        self.debug = debug

        match_filter.init_graphs(subgraph, graph)

    def __call__(self):
        node0 = self.subgraph.central_node
        for node1 in self.graph.nodes:
            if self.match_filter.compare(node0, node1):
                init_relations = [(node0, node1)]
                init_match = OneToOne(init_relations)
                for match in self.yield_matches(node0, node1, init_match, 1, init_relations):
                    yield match

    def combine(self, group0, group1):
        if len(group0) == 0 or len(group1) == 0:
            yield []
        else:
            for index1 in xrange(len(group1)):
                group1_copy = copy.copy(group1)
                pop1 = group1_copy.pop(index1)
                for tail in self.combine(group0[1:], group1_copy):
                    yield [(group0[0], pop1)] + tail

    def combine_relations(self, relation_groups):
        if len(relation_groups) == 0:
            yield []
        else:
            group0, group1 = relation_groups[0]
            for new_relations_rest in self.combine_relations(relation_groups[1:]):
                for new_relations_0 in self.combine(group0, group1):
                    yield new_relations_0 + new_relations_rest

    def yield_matches(self, node0, node1, input_match, shell_index, previous_relations):
        if self.debug: prefix = " "*shell_index
        if self.debug: print prefix, "ENTERING SHELL %i" % shell_index

        pool0, pool1 = self.match_filter.new_pools(node0, node1, input_match, shell_index)
        if len(pool0) == 0 or len(pool1) == 0:
            if self.debug: print prefix, "LEAVING SHELL (empty pool(s)) %i" % shell_index
            return

        cf = ClusterFactory()
        for previous_node0, previous_node1 in previous_relations:
            neighbors0 = pool0.intersection(self.subgraph.neighbors[previous_node0])
            neighbors1 = pool1.intersection(self.graph.neighbors[previous_node1])
            for neighbor0 in neighbors0:
                for neighbor1 in neighbors1:
                    if self.match_filter.compare(neighbor0, neighbor1):
                        cf.add_members([(0, neighbor0), (1, neighbor1)])

        clusters = cf.get_clusters()
        relation_groups = []
        for cluster in clusters:
            relation_group = [
                [neighbor for gid, neighbor in cluster.members if gid==0], # gid = graph index
                [neighbor for gid, neighbor in cluster.members if gid==1]
            ]
            if not self.match_filter.useful_relation_group(relation_group):
                continue
            if not self.match_filter.valid_relation_group(relation_group):
                if self.debug: print prefix, "LEAVING SHELL (invalid relation_groups) %i" % shell_index
                return
            relation_groups.append(relation_group)
        if len(relation_groups) == 0:
            if self.debug: print prefix, "LEAVING SHELL (zero relation_groups) %i" % shell_index
            return

        for new_relations in self.combine_relations(relation_groups):
            if self.match_filter.check_symmetry(new_relations):
                next_match = copy.copy(input_match)
                for new_node0, new_node1 in new_relations:
                    next_match.add_relation(new_node0, new_node1)
                if self.match_filter.complete(node0, node1, next_match, shell_index):
                    yield next_match
                else:
                    for match in self.yield_matches(node0, node1, next_match, shell_index+1, new_relations):
                        yield match

        if self.debug: print prefix, "LEAVING SHELL %i" % shell_index
