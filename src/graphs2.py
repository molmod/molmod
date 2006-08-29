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
        result.backward = self.backward.copy()
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

    def add_relations(self, pairs):
        for source, destination in pairs:
            self.add_relation(source, destination)

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

    def get_distance(self, node_a, node_b):
        return self.distances[self.index[node_a], self.index[node_b]]

    def yield_shortest_paths(self, node_a, node_b):
        if node_a == node_b:
            yield []
        else:
            for down_b in self.trees[node_a][node_b]:
                for path in self.yield_shortest_paths(node_a, down_b):
                    yield [down_b] + path


class Match(OneToOne):
    def __init__(self, init_relations):
        OneToOne.__init__(self, init_relations)
        self.previous_relations = init_relations

    def copy_with_new_relations(self, new_relations):
        result = self.__copy__()
        result.add_relations(new_relations)
        result.previous_relations = new_relations
        return result


class MatchDefinitionError(Exception):
    pass


class MatchDefinition(object):
    MatchClass = Match

    def init_graph(self, graph):
        "Checks initialy whether it makes sense to match both graphs."
        self.graph = graph

    def init_matches(self):
        "Yields the initial matches to start with."
        return

    def new_pools(self, partial_match):
        "Returns pool0 for subgraph and pool1 for graph that contain the nodes for the new relations."
        # derived classes should return both pool0 and pool1. Here only pool1
        # is implemented:
        return set(self.graph.nodes).difference(partial_match.backward.keys())

    def subgraph_neighbors(self, node0):
        return None

    def graph_neighbors(self, node1):
        return self.graph.neighbors[node1]

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

    def complete(self, match):
        return True


class SubGraphMatchDefinition(MatchDefinition):
    def __init__(self, subgraph):
        self.subgraph = subgraph

    def init_matches(self):
        node0 = self.subgraph.central_node
        for node1 in self.graph.nodes:
            yield self.MatchClass([(node0, node1)])

    def new_pools(self, partial_match):
        "Creates pool0 for subgraph and pool1 for graph that contain the nodes for the new relations."
        pool1 = MatchDefinition.new_pools(self, partial_match)
        # The default behavior is to allow all nodes that have not been used yet.
        pool0 = set(self.subgraph.nodes).difference(partial_match.forward.keys())
        return pool0, pool1

    def subgraph_neighbors(self, node0):
        return self.subgraph.neighbors[node0]


class ExactMatch(Match):
    def __init__(self, init_relations):
        assert len(init_relations) == 1, "An exact match filter can only start with one relation!"
        self.node0, self.node1 = init_relations[0]
        self.shell_index = 1
        Match.__init__(self, init_relations)

    def copy_with_new_relations(self, new_relations):
        result = Match.copy_with_new_relations(self, new_relations)
        result.shell_index += 1
        return result


class ExactMatchDefinition(SubGraphMatchDefinition):
    MatchClass = ExactMatch

    def init_graph(self, graph):
        if len(self.subgraph.nodes) != len(graph.nodes):
            raise MatchDefinitionError("It does not make sense to find an exact match between two graphs if the number of nodes is different")
        if len(self.subgraph.pairs) != len(graph.pairs):
            raise MatchDefinitionError("It does not make sense to find an exact match between two graphs if the number of relations is different")
        MatchDefinition.init_graph(self, graph)

    def new_pools(self, partial_match):
        # for an exact match, the matching nodes come from matching shells.
        shell_index = partial_match.shell_index
        node0 = partial_match.node0
        node1 = partial_match.node1
        if len(self.subgraph.shells[node0]) <= shell_index:
            pool0 = set([])
        else:
            pool0 = self.subgraph.shells[node0][shell_index]
        if len(self.graph.shells[node1]) <= shell_index:
            pool1 = set([])
        else:
            pool1 = self.graph.shells[node1][shell_index]
        return pool0, pool1

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

    def complete(self, match):
        return len(self.subgraph.nodes) == len(match.forward)


class EgoMatchDefinition(ExactMatchDefinition):
    def __init__(self):
        ExactMatchDefinition.__init__(self, None)

    def init_graph(self, graph):
        self.subgraph = graph
        ExactMatchDefinition.init_graph(self, graph)


class RingMatchError(Exception):
    pass


class RingMatch(Match):
    def __init__(self, init_relations):
        assert len(init_relations) == 1, "A ring match filter can only start with one relation!"
        self.increasing = True
        self.first_node = init_relations[0][1]
        self.size = None
        self.length = 0
        self.cache = {}
        Match.__init__(self, init_relations)

    def cache_node(self, node, increasing=False, decreasing=False):
        self.cache[node] = (increasing, decreasing)

    def check_oposite_even(self, node, graph):
        if self.size is None:
            max_distance = self.length
        else:
            max_distance = self.size/2
        oposite_node = self.forward[self.length + 1 - max_distance]
        #print node, oposite_node
        #print max_distance, graph.get_distance(oposite_node, node)
        if not graph.get_distance(oposite_node, node) == max_distance:
            return False
        counter = 0
        for path in graph.yield_shortest_paths(node, oposite_node):
            counter += 1
            if counter > 2:
                return False
        return True

    def check_oposite_odd(self, node, graph):
        if self.size is None:
            max_distance = self.length
        else:
            max_distance = self.size/2
        oposite_node1 = self.forward[self.length + 1 - max_distance]
        oposite_node2 = self.forward[self.length + 1 - max_distance - 1]
        #print node, oposite_node1, oposite_node2
        #print max_distance, graph.get_distance(oposite_node1, node), graph.get_distance(oposite_node2, node)
        if not (graph.get_distance(oposite_node1, node) == max_distance and
                graph.get_distance(oposite_node2, node) == max_distance):
            return False
        counter = 0
        for path in graph.yield_shortest_paths(node, oposite_node1):
            counter += 1
            if counter > 1:
                return False
        counter = 0
        for path in graph.yield_shortest_paths(node, oposite_node2):
            counter += 1
            if counter > 1:
                return False
        return True

    def check_distances(self, node, graph, max_size):
        new_distance = graph.get_distance(node, self.first_node)
        last_distance = graph.get_distance(self.forward[len(self)-1], self.first_node)
        #print "new-, last-distance", new_distance, last_distance
        if self.increasing:
            if new_distance > max_size:
                return False
            if new_distance == len(self):
                self.cache_node(node, increasing=True)
                return True
            elif new_distance == len(self) - 1:
                if self.check_oposite_odd(node, graph):
                    self.cache_node(node)
                    return True
                else:
                    return False
            elif new_distance == len(self) - 2:
                if self.check_oposite_even(node, graph):
                    self.cache_node(node, decreasing=True)
                    return True
                else:
                    return False
            else:
                raise RingMatchError("Distances can only vary by one at a time. There must be an error in the graph.distances matrix.")
        else:
            # do distance checks
            if new_distance != last_distance - 1:
                return False
            if self.size % 2 == 0:
                if self.check_oposite_even(node, graph):
                    self.cache_node(node, decreasing=True)
                    return True
                else:
                    return False
            else:
                if self.check_oposite_odd(node, graph):
                    self.cache_node(node, decreasing=True)
                    return True
                else:
                    return False

    def copy_with_new_relations(self, new_relations):
        assert len(new_relations) == 1, "Only one new relation at a time."
        result = Match.copy_with_new_relations(self, new_relations)
        result.cache = {}
        result.length += 1

        node1 = new_relations[0][1]
        increasing, decreasing = self.cache[node1]
        if not increasing and self.size == None:
            #print "not increasing"
            result.increasing = False
            if decreasing:
                result.size = self.length * 2
            else:
                result.size = self.length * 2 + 1
            #print "SIZE", result.size
        return result



class RingMatchDefinition(MatchDefinition):
    MatchClass = RingMatch

    def __init__(self, max_size):
        self.max_size = max_size

    def init_matches(self):
        node0 = 0
        for node1 in self.graph.nodes:
            yield self.MatchClass([(node0, node1)])

    def new_pools(self, partial_match):
        first_node1 = partial_match.forward[0]
        last_node1 = partial_match.forward[len(partial_match)-1]
        #pool1 = MatchDefinition.new_pools(self, partial_match)
        pool1 = set([
            node1
            for node1
            in self.graph.neighbors[last_node1]
            if (
                node1 > partial_match.first_node and
                node1 not in partial_match.backward and
                partial_match.check_distances(node1, self.graph, self.max_size)
           )
        ])
        pool0 = set([len(partial_match)])
        return pool0, pool1

    def subgraph_neighbors(self, node0):
        return [node0 + 1]

    def compare(self, node0, node1):
        return True

    def complete(self, match):
        # TODO: a more efficient way of eliminating duplicate rings (that only
        # differ in the direction
        return (match.size == match.length + 1) and match.forward[1] < match.forward[match.length]


class MatchGenerator(object):
    def __init__(self, match_definition, graph, debug=False):
        self.match_definition = match_definition
        self.debug = debug

        match_definition.init_graph(graph)

    def __call__(self):
        for init_match in self.match_definition.init_matches():
            for match in self.yield_matches(init_match):
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

    def print_debug(self, text, indent=0):
        if self.debug:
            if indent > 0:
                print " "*self.debug, text
            self.debug += indent
            if indent <= 0:
                print " "*self.debug, text

    def yield_matches(self, input_match):
        self.print_debug("ENTERING YIELD_MATCHES", 1)
        self.print_debug("input_match: %s" % input_match)
        pool0, pool1 = self.match_definition.new_pools(input_match)
        if len(pool0) == 0 or len(pool1) == 0:
            self.print_debug("LEAVING YIELD_MATCHES (empty pools)", -1)
            return

        cf = ClusterFactory()
        for previous_node0, previous_node1 in input_match.previous_relations:
            neighbors0 = pool0.intersection(self.match_definition.subgraph_neighbors(previous_node0))
            neighbors1 = pool1.intersection(self.match_definition.graph_neighbors(previous_node1))
            for neighbor0 in neighbors0:
                for neighbor1 in neighbors1:
                    if self.match_definition.compare(neighbor0, neighbor1):
                        cf.add_members([(0, neighbor0), (1, neighbor1)])

        clusters = cf.get_clusters()
        relation_groups = []
        for cluster in clusters:
            relation_group = [
                [neighbor for gid, neighbor in cluster.members if gid==0], # gid = graph index
                [neighbor for gid, neighbor in cluster.members if gid==1]
            ]
            if not self.match_definition.useful_relation_group(relation_group):
                continue
            if not self.match_definition.valid_relation_group(relation_group):
                self.print_debug("LEAVING YIELD_MATCHES (invalid_relation_groups)", -1)
                return
            relation_groups.append(relation_group)
        if len(relation_groups) == 0:
            self.print_debug("LEAVING YIELD_MATCHES (no new relations)", -1)
            return

        self.print_debug("relation_groups: %s" % relation_groups)
        for new_relations in self.combine_relations(relation_groups):
            if self.match_definition.check_symmetry(new_relations):
                next_match = input_match.copy_with_new_relations(new_relations)
                if self.match_definition.complete(next_match):
                    yield next_match
                else:
                    for match in self.yield_matches(next_match):
                        yield match

        self.print_debug("LEAVING YIELD_MATCHES", -1)
