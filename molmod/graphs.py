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
#
"""Graph data structure, analysis and pattern matching.

   Some distinctive features include:

   * Iterating over all shortest paths between two vertices (Dijkstra). See
     http://en.wikipedia.org/wiki/Dijkstra's_algorithm for more info.
   * Iterating over vertices or edges using the Breadth First convention. See
     http://en.wikipedia.org/wiki/Breadth-first_search for more info.
   * The all pairs shortest path matrix (Floyd-Warshall). See
     http://en.wikipedia.org/wiki/Floyd-Warshall_algorithm for more info.
   * Symmetry analysis of graphs (automorphisms). The Graph class can generate a
     list of permutations between vertices that map the graph onto itself. This
     can be used to generate (and test) all possible geometric symmetries in a
     molecule. See http://en.wikipedia.org/wiki/Graph_automorphism for more
     info.
   * Scanning a graph for patterns.
     The GraphSearch is a generic class that can scan a graph for certain
     patterns, e.g. given pattern_graphs, strong rings, isomorphisms,
     automorphisms, ... The pattern_graph can deal with (multiple sets of)
     additional conditions that must be satisfied, such as "give me all
     dihedral angles where the central atoms are carbons" without duplicates.

   The central class in this module is 'Graph'. It caches most of the analysis
   results, which implies that the graph structure can not be changed once the
   object is created. If you feel the need to do this, just construct a new
   graph object.
"""


from __future__ import print_function, division

from builtins import range
import copy

import numpy as np

from molmod.utils import cached, ReadOnly, ReadOnlyAttribute


__all__ = [
    "GraphError", "Graph", "OneToOne", "Match", "Pattern",
    "CriteriaSet", "Anything", "CritOr", "CritAnd", "CritXor", "CritNot",
    "CustomPattern", "EqualPattern", "RingPattern", "GraphSearch",
]


class GraphError(Exception):
    """Raised when something goes wrong in one of the graph algorithms"""
    pass


class Graph(ReadOnly):
    """An undirected graph, where edges have equal weight

       The edges attribute is the most elementary and is always available. All
       other attributes are optional.

       Graphs are meant to be immutable objects. Once created they are not
       supposed to be modified. If you want an adapted graph, create a new
       object. The main reason is that all graph analysis work is cached in this
       object. When the edges change, the cache becomes invalid and should be
       erased. The latter is not supported as it is easier to create just a new
       graph object with the modified edges.

       If you want to associate data with vertices or edges, add custom
       dictionary or list attributes to the graph object. This is the way to
       work attributes that are not applicable to all vertices or edges:

       >>> graph.vertex_property = {vertex1: blablabla, vertex2: blablabla, ...}
       >>> graph.edge_property = {frozenset([vertex1, vertex2]): blablabla, ..}

       If a certain property applies to all vertices or edges, it is sometimes
       more practical to work with lists or ``numpy`` arrays that have the same
       ordering as the vertices or the edges:

       >>> # atom numbers of ethene
       >>> graph.vertex_property = np.array([6, 6, 1, 1, 1, 1], int)
       >>> # bond orders of ethene
       >>> graph.edge_property = np.array([2, 1, 1, 1, 1], int)
    """
    edges = ReadOnlyAttribute(tuple, none=False, doc="the incidence list")
    num_vertices = ReadOnlyAttribute(int, none=False, doc="the number of vertices")

    def __init__(self, edges, num_vertices=None):
        """
           Arguments:
            | ``edges`` -- ``tuple(frozenset([vertex1, vertex2]), ...)``
            | ``num_vertices`` -- number of vertices

           vertex1 to vertexN must be integers from 0 to N-1. If vertices above
           the highest vertex value are not connected by edges, use the
           num_vertices argument to tell what the total number of vertices is.

           If the edges argument does not have the correct format, it will be
           converted.
        """

        tmp = []
        for edge in edges:
            if len(edge) != 2:
                raise TypeError("The edges must be a iterable with 2 elements")
            i, j = edge
            i = int(i)
            j = int(j)
            if i == j:
                raise ValueError("A edge must contain two different values.")
            if i < 0 or j < 0:
                raise TypeError("The edges must contain positive integers.")
            tmp.append(frozenset([i, j]))
        edges = tuple(tmp)

        if len(edges) == 0:
            real_num_vertices = 0
        else:
            real_num_vertices = max(max(a, b) for a, b in edges)+1
        if num_vertices is not None:
            if not isinstance(num_vertices, int):
                raise TypeError("The optional argument num_vertices must be an "
                    "integer when given.")
            if num_vertices < real_num_vertices:
                raise ValueError("num_vertices must be equal or larger to the "
                    "number of vertices deduced from the edge list.")
            real_num_vertices = num_vertices

        self.edges = edges
        self.num_vertices = real_num_vertices

    num_edges = property(lambda self: len(self.edges),
        doc="*Read-only attribute:* the number of edges in the graph.")

    def __mul__(self, repeat):
        """Construct a graph that repeats this graph a number of times

           Arguments:
            | ``repeat`` -- The number of repetitions.
        """
        if not isinstance(repeat, int):
            raise TypeError("Can only multiply a graph with an integer")
        new_edges = []
        for i in range(repeat):
            for vertex1, vertex2 in self.edges:
                new_edges.append(frozenset([
                    vertex1+i*self.num_vertices,
                    vertex2+i*self.num_vertices
                ]))
        return Graph(new_edges, self.num_vertices*repeat)

    __rmul__ = __mul__

    def __str__(self):
        return " ".join("%i-%i" % tuple(sorted([i, j])) for i, j in self.edges)

    # functions that should be implemented by derived classes

    def get_vertex_string(self, i):
        """Returns a string that fully characterizes the nature of the vertex

           In case of an ordinary graph, all vertices are the same.
        """
        return ""

    def get_edge_string(self, i):
        """Returns a string that fully characterizes the nature of an edge

           In case of an ordinary graph, all edges are the same.
        """
        return ""

    # cached attributes:

    @cached
    def edge_index(self):
        """A map to look up the index of a edge"""
        return dict((edge, index) for index, edge in enumerate(self.edges))

    @cached
    def neighbors(self):
        """A dictionary with neighbors

           The dictionary will have the following form:
           ``{vertexX: (vertexY1, vertexY2, ...), ...}``
           This means that vertexX and vertexY1 are connected etc. This also
           implies that the following elements are part of the dictionary:
           ``{vertexY1: (vertexX, ...), vertexY2: (vertexX, ...), ...}``.
        """
        neighbors = dict(
            (vertex, []) for vertex
            in range(self.num_vertices)
        )
        for a, b in self.edges:
            neighbors[a].append(b)
            neighbors[b].append(a)
        # turn lists into frozensets
        neighbors = dict((key, frozenset(val)) for key, val in neighbors.items())
        return neighbors

    @cached
    def distances(self):
        """The matrix with the all-pairs shortest path lenghts"""
        from molmod.ext import graphs_floyd_warshall
        distances = np.zeros((self.num_vertices,)*2, dtype=int)
        #distances[:] = -1 # set all -1, which is just a very big integer
        #distances.ravel()[::len(distances)+1] = 0 # set diagonal to zero
        for i, j in self.edges: # set edges to one
            distances[i, j] = 1
            distances[j, i] = 1
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
    def central_vertices(self):
        """Vertices that have the lowest maximum distance to any other vertex"""
        max_distances = self.distances.max(0)
        max_distances_min = max_distances[max_distances > 0].min()
        return (max_distances == max_distances_min).nonzero()[0]

    @cached
    def central_vertex(self):
        """The vertex that has the lowest maximum distance to any other vertex

           This definition does not lead to a unique solution. One arbitrary
           solution is selected.
        """
        return self.central_vertices[0]

    @cached
    def independent_vertices(self):
        """Lists of vertices that are only interconnected within each list

           This means that there is no path from a vertex in one list to a
           vertex in another list. In case of a molecular graph, this would
           yield the atoms that belong to individual molecules.
        """
        candidates = set(range(self.num_vertices))

        result = []
        while len(candidates) > 0:
            pivot = candidates.pop()
            group = [
                vertex for vertex, distance
                in self.iter_breadth_first(pivot)
            ]
            candidates.difference_update(group)

            # this sort makes sure that the order of the vertices is respected
            group.sort()
            result.append(group)
        return result

    @cached
    def fingerprint(self):
        """A total graph fingerprint

           The result is invariant under permutation of the vertex indexes. The
           chance that two different (molecular) graphs yield the same
           fingerprint is small but not zero. (See unit tests.)"""
        if self.num_vertices == 0:
            return np.zeros(20, np.ubyte)
        else:
            return sum(self.vertex_fingerprints)

    @cached
    def vertex_fingerprints(self):
        """A fingerprint for each vertex

           The result is invariant under permutation of the vertex indexes.
           Vertices that are symmetrically equivalent will get the same
           fingerprint, e.g. the hydrogens in methane would get the same
           fingerprint.
        """
        return self.get_vertex_fingerprints(
            [self.get_vertex_string(i) for i in range(self.num_vertices)],
            [self.get_edge_string(i) for i in range(self.num_edges)],
        )

    @cached
    def equivalent_vertices(self):
        """A dictionary with symmetrically equivalent vertices."""
        level1 = {}
        for i, row in enumerate(self.vertex_fingerprints):
            key = row.tobytes()
            l = level1.get(key)
            if l is None:
                l = set([i])
                level1[key] = l
            else:
                l.add(i)
        level2 = {}
        for key, vertices in level1.items():
            for vertex in vertices:
                level2[vertex] = vertices
        return level2

    @cached
    def symmetries(self):
        """Graph symmetries (permutations) that map the graph onto itself."""

        symmetry_cycles = set([])
        symmetries = set([])
        for match in GraphSearch(EqualPattern(self))(self):
            match.cycles = match.get_closed_cycles()
            if match.cycles in symmetry_cycles:
                raise RuntimeError("Duplicates in EqualMatch")
            symmetry_cycles.add(match.cycles)
            symmetries.add(match)
        return symmetries

    @cached
    def symmetry_cycles(self):
        """The cycle representations of the graph symmetries"""
        result = set([])
        for symmetry in self.symmetries:
            result.add(symmetry.cycles)
        return result

    @cached
    def canonical_order(self):
        """The vertices in a canonical or normalized order.

           This routine will return a list of vertices in an order that does not
           depend on the initial order, but only depends on the connectivity and
           the return values of the function self.get_vertex_string.

           Only the vertices that are involved in edges will be included. The
           result can be given as first argument to self.get_subgraph, with
           reduce=True as second argument. This will return a complete canonical
           graph.

           The routine is designed not to use symmetry relations that are
           obtained with the GraphSearch routine. We also tried to create an
           ordering that feels like natural, i.e. starting in the center and
           pushing vertices with few equivalents to the front. If necessary, the
           nature of the vertices and  their bonds to atoms closer to the center
           will also play a role, but only as a last resort.
        """
        # A) find an appropriate starting vertex.
        # Here we take a central vertex that has a minimal number of symmetrical
        # equivalents, 'the highest atom number', and the highest fingerprint.
        # Note that the symmetrical equivalents are computed from the vertex
        # fingerprints, i.e. without the GraphSearch.
        starting_vertex = max(
            (
                -len(self.equivalent_vertices[vertex]),
                self.get_vertex_string(vertex),
                self.vertex_fingerprints[vertex].tobytes(),
                vertex
            ) for vertex in self.central_vertices
        )[-1]

        # B) sort all vertices based on
        #      1) distance from central vertex
        #      2) number of equivalent vertices
        #      3) vertex string, (higher atom numbers come first)
        #      4) fingerprint
        #      5) vertex index
        # The last field is only included to collect the result of the sort.
        # The fingerprint on itself would be sufficient, but the three first are
        # there to have a naturally appealing result.
        l = [
            [
                -distance,
                -len(self.equivalent_vertices[vertex]),
                self.get_vertex_string(vertex),
                self.vertex_fingerprints[vertex].tobytes(),
                vertex
            ] for vertex, distance in self.iter_breadth_first(starting_vertex)
            if len(self.neighbors[vertex]) > 0
        ]
        l.sort(reverse=True)

        # C) The order of some vertices is still not completely set. e.g.
        # consider the case of allene. The four hydrogen atoms are equivalent,
        # but one can have two different orders: make geminiles consecutive or
        # don't. It is more trikcy than one would think at first sight. In the
        # case of allene, geminility could easily solve the problem. Consider a
        # big flat rotationally symmetric molecule (order 2). The first five
        # shells are order 4 and one would just give a random order to four
        # segemnts in the first shell. Only when one reaches the outer part that
        # has order two, it turns out that the arbitrary choices in the inner
        # shell play a role. So it does not help to look at relations with
        # vertices at inner or current shells only. One has to consider the
        # whole picture. (unit testing reveals troubles like these)

        # I need some sleep now. The code below checks for potential fuzz and
        # will raise an error if the ordering is not fully determined yet. One
        # day, I'll need this code more than I do now, and I'll fix things up.
        # I know how to do this, but I don't care enough right now.
        # -- Toon
        for i in range(1, len(l)):
            if l[i][:-1] == l[i-1][:-1]:
                raise NotImplementedError

        # D) Return only the vertex indexes.
        return [record[-1] for record in l]

    # other usefull graph functions

    def iter_breadth_first(self, start=None, do_paths=False, do_duplicates=False):
        """Iterate over the vertices with the breadth first algorithm.

           See http://en.wikipedia.org/wiki/Breadth-first_search for more info.
           If not start vertex is given, the central vertex is taken.

           By default, the distance to the starting vertex is also computed. If
           the path to the starting vertex should be computed instead, set path
           to True.

           When duplicate is True, then vertices that can be reached through
           different  paths of equal length, will be iterated twice. This
           typically only makes sense when path==True.
        """
        if start is None:
            start = self.central_vertex
        else:
            try:
                start = int(start)
            except ValueError:
                raise TypeError("First argument (start) must be an integer.")
            if start < 0 or start >= self.num_vertices:
                raise ValueError("start must be in the range [0, %i[" %
                                 self.num_vertices)
        from collections import deque
        work = np.zeros(self.num_vertices, int)
        work[:] = -1
        work[start] = 0
        if do_paths:
            result = (start, 0, (start, ))
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
                        current_path = parent_path + (current, )
                        result = (current, current_length, current_path)
                    else:
                        result = (current, current_length)
                    #print "iter_breadth_first", result
                    yield result
                    todo.append(result)

    def iter_shortest_paths(self, a, b):
        """Iterate over all the shortest paths between vertex a and b."""
        max_len = None
        ibf = self.iter_breadth_first(a, do_paths=True, do_duplicates=True)
        for vertex, length, path in ibf:
            if max_len is not None:
                if length > max_len:
                    return
            if vertex == b:
                max_len = length
                yield path

    def iter_breadth_first_edges(self, start=None):
        """Iterate over the edges with the breadth first convention.

           We need this for the pattern matching algorithms, but a quick look at
           Wikipedia did not result in a known and named algorithm.

           The edges are yielded one by one, together with the distance of the
           edge from the starting vertex and a flag that indicates whether the
           yielded edge connects two vertices that are at the same distance from
           the starting vertex. If that flag is False, the distance from the
           starting vertex to edge[0] is equal to the distance variable and the
           distance from edge[1] to the starting vertex is equal to distance+1.
           One item has the following format: ((i, j), distance, flag)
        """
        if start is None:
            start = self.central_vertex
        else:
            try:
                start = int(start)
            except ValueError:
                raise TypeError("First argument (start) must be an integer.")
            if start < 0 or start >= self.num_vertices:
                raise ValueError("start must be in the range [0, %i[" %
                                 self.num_vertices)
        from collections import deque
        work = np.zeros(self.num_vertices, int)
        work[:] = -1
        work[start] = 0
        todo = deque([start])
        while len(todo) > 0:
            parent = todo.popleft()
            distance = work[parent]
            for current in self.neighbors[parent]:
                if work[current] == -1:
                    yield (parent, current), distance, False
                    work[current] = distance+1
                    todo.append(current)
                elif work[current] == distance and current > parent:
                    # second equation in elif avoids duplicates
                    yield (parent, current), distance, True
                elif work[current] == distance+1:
                    yield (parent, current), distance, False

    def get_subgraph(self, subvertices, normalize=False):
        """Constructs a subgraph of the current graph

           Arguments:
            | ``subvertices`` -- The vertices that should be retained.
            | ``normalize`` -- Whether or not the vertices should renumbered and
                 reduced to the given set of subvertices. When True, also the
                 edges are sorted. It the end, this means that new order of the
                 edges does not depend on the original order, but only on the
                 order of the argument subvertices.
                 This option is False by default. When False, only edges will be
                 discarded, but the retained data remain unchanged. Also the
                 parameter num_vertices is not affected.

           The returned graph will have an attribute ``old_edge_indexes`` that
           relates the positions of the new and the old edges as follows::

             >>> self.edges[result._old_edge_indexes[i]] = result.edges[i]

           In derived classes, the following should be supported::

             >>> self.edge_property[result._old_edge_indexes[i]] = result.edge_property[i]

           When ``normalize==True``, also the vertices are affected and the
           derived classes should make sure that the following works::

             >>> self.vertex_property[result._old_vertex_indexes[i]] = result.vertex_property[i]

           The attribute ``old_vertex_indexes`` is only constructed when
           ``normalize==True``.
        """
        if normalize:
            revorder = dict((j, i) for i, j in enumerate(subvertices))
            new_edges = []
            old_edge_indexes = []
            for counter, (i, j) in enumerate(self.edges):
                new_i = revorder.get(i)
                if new_i is None:
                    continue
                new_j = revorder.get(j)
                if new_j is None:
                    continue
                new_edges.append((new_i, new_j))
                old_edge_indexes.append(counter)
            # sort the edges
            order = list(range(len(new_edges)))
            # argsort in pure python
            order.sort( key=(lambda i: tuple(sorted(new_edges[i]))) )
            new_edges = [new_edges[i] for i in order]
            old_edge_indexes = [old_edge_indexes[i] for i in order]

            result = Graph(new_edges, num_vertices=len(subvertices))
            result._old_vertex_indexes = np.array(subvertices, dtype=int)
            #result.new_vertex_indexes = revorder
            result._old_edge_indexes = np.array(old_edge_indexes, dtype=int)
        else:
            subvertices = set(subvertices)
            old_edge_indexes = np.array([
                i for i, edge in enumerate(self.edges)
                if edge.issubset(subvertices)
            ], dtype=int)
            new_edges = tuple(self.edges[i] for i in old_edge_indexes)
            result = Graph(new_edges, self.num_vertices)
            result._old_edge_indexes = old_edge_indexes
            # no need for old and new vertex_indexes because they remain the
            # same.
        return result

    def get_vertex_fingerprints(self, vertex_strings, edge_strings, num_iter=None):
        """Return an array with fingerprints for each vertex"""
        import hashlib
        def str2array(x):
            """convert a hash string to a numpy array of bytes"""
            if len(x) == 0:
                return np.zeros(0, np.ubyte)
            else:
                return np.fromstring(x, np.ubyte)
        hashrow = lambda x: np.fromstring(hashlib.sha1(x.data).digest(), np.ubyte)
        # initialization
        result = np.zeros((self.num_vertices, 20), np.ubyte)
        for i in range(self.num_vertices):
            result[i] = hashrow(str2array(vertex_strings[i]))
        for i in range(self.num_edges):
            a, b = self.edges[i]
            tmp = hashrow(str2array(edge_strings[i]))
            result[a] += tmp
            result[b] += tmp
        work = result.copy()
        # iterations
        if num_iter is None:
            num_iter = self.max_distance
        for i in range(num_iter):
            for a, b in self.edges:
                work[a] += result[b]
                work[b] += result[a]
            #for a in xrange(self.num_vertices):
            #    for b in xrange(self.num_vertices):
            #        work[a] += hashrow(result[b]*self.distances[a, b])
            for a in range(self.num_vertices):
                result[a] = hashrow(work[a])
        return result

    def get_halfs(self, vertex1, vertex2):
        """Split the graph in two halfs by cutting the edge: vertex1-vertex2

           If this is not possible (due to loops connecting both ends), a
           GraphError is raised.

           Returns the vertices in both halfs.
        """
        def grow(origin, other):
            frontier = set(self.neighbors[origin])
            frontier.discard(other)
            result = set([origin])
            while len(frontier) > 0:
                pivot = frontier.pop()
                if pivot == other:
                    raise GraphError("The graph can not be separated in two halfs "
                                     "by disconnecting vertex1 and vertex2.")
                pivot_neighbors = set(self.neighbors[pivot])
                pivot_neighbors -= result
                frontier |= pivot_neighbors
                result.add(pivot)
            return result

        vertex1_part = grow(vertex1, vertex2)
        vertex2_part = grow(vertex2, vertex1)
        return vertex1_part, vertex2_part

    def get_part(self, vertex_in, vertices_border):
        """List all vertices that are connected to vertex_in, but are not
           included in or 'behind' vertices_border.
        """
        vertices_new = set(self.neighbors[vertex_in])
        vertices_part = set([vertex_in])

        while len(vertices_new) > 0:
            pivot = vertices_new.pop()
            if pivot in vertices_border:
                continue
            vertices_part.add(pivot)
            pivot_neighbors = set(self.neighbors[pivot])
            pivot_neighbors -= vertices_part
            vertices_new |= pivot_neighbors

        return vertices_part

    def get_halfs_double(self, vertex_a1, vertex_b1, vertex_a2, vertex_b2):
        """Compute the two parts separated by ``(vertex_a1, vertex_b1)`` and ``(vertex_a2, vertex_b2)``

           Raise a GraphError when ``(vertex_a1, vertex_b1)`` and
           ``(vertex_a2, vertex_b2)`` do not separate the graph in two
           disconnected parts. The edges must be neighbors. If not a GraphError
           is raised. The for vertices must not coincide or a GraphError is
           raised.

           Returns the vertices of the two halfs and the four 'hinge' vertices
           in the correct order, i.e. both ``vertex_a1`` and ``vertex_a2`` are
           in the first half and both ``vertex_b1`` and ``vertex_b2`` are in the
           second half.
        """
        if vertex_a1 not in self.neighbors[vertex_b1]:
            raise GraphError("vertex_a1 must be a neighbor of vertex_b1.")
        if vertex_a2 not in self.neighbors[vertex_b2]:
            raise GraphError("vertex_a2 must be a neighbor of vertex_b2.")

        # find vertex_a_part (and possibly switch vertex_a2, vertex_b2)
        vertex_a_new = set(self.neighbors[vertex_a1])
        vertex_a_new.discard(vertex_b1)
        if vertex_a1 == vertex_b2:
            # we now that we have to swap vertex_a2 and vertex_b2. The algo
            # below will fail otherwise in this 'exotic' case.
            vertex_a2, vertex_b2 = vertex_b2, vertex_a2
            #vertex_a_new.discard(vertex_a2) # in case there is overlap
        if vertex_a1 == vertex_a2:
            vertex_a_new.discard(vertex_b2) # in case there is overlap
        vertex_a_part = set([vertex_a1])

        touched = False # True if (the switched) vertex_a2 has been reached.
        while len(vertex_a_new) > 0:
            pivot = vertex_a_new.pop()
            if pivot == vertex_b1:
                raise GraphError("The graph can not be separated in two halfs. "
                                 "vertex_b1 reached by vertex_a1.")
            vertex_a_part.add(pivot)
            # create a new set that we can modify
            pivot_neighbors = set(self.neighbors[pivot])
            pivot_neighbors -= vertex_a_part
            if pivot == vertex_a2 or pivot == vertex_b2:
                if pivot == vertex_b2:
                    if touched:
                        raise GraphError("The graph can not be separated in "
                                         "two halfs. vertex_b2 reached by "
                                         "vertex_a1.")
                    else:
                        # put them in the correct order
                        vertex_a2, vertex_b2 = vertex_b2, vertex_a2
                pivot_neighbors.discard(vertex_b2)
                touched = True
            vertex_a_new |= pivot_neighbors

        if vertex_a2 not in vertex_a_part:
            raise GraphError("The graph can not be separated in two halfs. "
                             "vertex_a1 can not reach vertex_a2 trough "
                             "vertex_a_part")

        # find vertex_b_part: easy, is just the rest ...
        #vertex_b_part = set(xrange(self.num_vertices)) - vertex_a_part

        # ... but we also want that there is a path in vertex_b_part from
        # vertex_b1 to vertex_b2
        if vertex_b1 == vertex_b2:
            closed = True
        else:
            vertex_b_new = set(self.neighbors[vertex_b1])
            vertex_b_new.discard(vertex_a1)
            vertex_b_part = set([vertex_b1])

            closed = False
            while len(vertex_b_new) > 0:
                pivot = vertex_b_new.pop()
                if pivot == vertex_b2:
                    closed = True
                    break
                pivot_neighbors = set(self.neighbors[pivot])
                pivot_neighbors -= vertex_b_part
                vertex_b_new |= pivot_neighbors
                vertex_b_part.add(pivot)

        if not closed:
            raise GraphError("The graph can not be separated in two halfs. "
                             "vertex_b1 can not reach vertex_b2 trough "
                             "vertex_b_part")

        # finaly compute the real vertex_b_part, the former loop might break
        # early for efficiency.
        vertex_b_part = set(range(self.num_vertices)) - vertex_a_part

        # done!
        return vertex_a_part, vertex_b_part, \
               (vertex_a1, vertex_b1, vertex_a2, vertex_b2)

    def full_match(self, other):
        """Find the mapping between vertex indexes in self and other.

           This also works on disconnected graphs. Derived classes should just
           implement get_vertex_string and get_edge_string to make this method
           aware of the different nature of certain vertices. In case molecules,
           this would make the algorithm sensitive to atom numbers etc.
        """
        # we need normalize subgraphs because these graphs are used as patterns.
        graphs0 = [
            self.get_subgraph(group, normalize=True)
            for group in self.independent_vertices
        ]
        graphs1 = [
            other.get_subgraph(group)
            for group in other.independent_vertices
        ]

        if len(graphs0) != len(graphs1):
            return

        matches = []

        for graph0 in graphs0:
            pattern = EqualPattern(graph0)
            found_match = False
            for i, graph1 in enumerate(graphs1):
                local_matches = list(GraphSearch(pattern)(graph1, one_match=True))
                if len(local_matches) == 1:
                    match = local_matches[0]
                    # we need to restore the relation between the normalized
                    # graph0 and its original indexes
                    old_to_new = OneToOne((
                        (j, i) for i, j
                        in enumerate(graph0._old_vertex_indexes)
                    ))
                    matches.append(match * old_to_new)
                    del graphs1[i]
                    found_match = True
                    break
            if not found_match:
                return

        result = OneToOne()
        for match in matches:
            result.add_relations(match.forward.items())
        return result


# Pattern matching


class OneToOne(object):
    """Implements a discrete bijection between source and destination elements

       The implementation is based on a consistent set of forward and reverse
       relations, stored in dictionaries.
    """

    def __init__(self, relations=None):
        """
           Argument:
            | ``relations``  --  initial relations for the bijection
        """
        self.forward = {}
        self.reverse = {}
        if relations is not None:
            self.add_relations(relations)

    def __len__(self):
        return len(self.forward)

    def __str__(self):
        result = "|"
        for source, destination in self.forward.items():
            result += " %s -> %s |" % (source, destination)
        return result

    def __mul__(self, other):
        """Return the result of the 'after' operator."""
        result = OneToOne()
        for source, mid in other.forward.items():
            destination = self.forward[mid]
            result.forward[source] = destination
            result.reverse[destination] = source
        return result

    def add_relation(self, source, destination):
        """Add new a relation to the bejection"""
        if self.in_sources(source):
            if self.forward[source] != destination:
                raise ValueError("Source is already in use. Destination does "
                                 "not match.")
            else:
                raise ValueError("Source-Destination relation already exists.")
        elif self.in_destinations(destination):
            raise ValueError("Destination is already in use. Source does not "
                             "match.")
        else:
            self.forward[source] = destination
            self.reverse[destination] = source

    def add_relations(self, relations):
        """Add multiple relations to a bijection"""
        for source, destination in relations:
            self.add_relation(source, destination)

    def get_destination(self, source):
        """Get the end point of a relation that start with 'source'"""
        return self.forward[source]

    def get_source(self, destination):
        """Get the starting point of a relation that ends with 'destination'"""
        return self.reverse[destination]

    def in_destinations(self, destination):
        """Test if the given destination is present"""
        return destination in self.reverse

    def in_sources(self, source):
        """Test if the given source is present"""
        return source in self.forward

    def inverse(self):
        """Returns the inverse bijection."""
        result = self.__class__()
        result.forward = copy.copy(self.reverse)
        result.reverse = copy.copy(self.forward)
        return result


class Match(OneToOne):
    """A match between a pattern and a graph"""
    @classmethod
    def from_first_relation(cls, vertex0, vertex1):
        """Intialize a fresh match based on the first relation"""
        result = cls([(vertex0, vertex1)])
        result.previous_ends1 = set([vertex1])
        return result

    def get_new_edges(self, subject_graph):
        """Get new edges from the subject graph for the graph search algorithm

           The Graph search algorithm extends the matches iteratively by adding
           matching vertices that are one edge further from the starting vertex
           at each iteration.
        """
        result = []
        #print "Match.get_new_edges self.previous_ends1", self.previous_ends1
        for vertex in self.previous_ends1:
            for neighbor in subject_graph.neighbors[vertex]:
                if neighbor not in self.reverse:
                    result.append((vertex, neighbor))
        return result

    def copy_with_new_relations(self, new_relations):
        """Create a new match object extended with new relations"""
        result = self.__class__(self.forward.items())
        result.add_relations(new_relations.items())
        result.previous_ends1 = set(new_relations.values())
        return result


class Pattern(object):
    """Base class for a pattern in a graph.

       Note the following conventions:

       * A pattern can always be represented by a graph (or a set of graphs)
         and some additional conditions. This graph is the so called 'PATTERN
         GRAPH'. For technical reasons, this pattern graph is not always
         constructed explicitly. Variables related to this graph often get
         suffix '0'. Note that a pattern graph is always fully connected.
       * The graph in which we search for the pattern, is called the 'SUBJECT
         GRAPH'. Variables related to this graph often get suffix '1'.
    """

    # This means that matching vertices must not have equal number of neighbors:
    sub = True
    MatchClass = Match

    def iter_initial_relations(self, subject_graph):
        """Iterate over all initial relations to start a Graph Search

           The function iterates of single relations between a pattern vertex
           and a subject vertex. In practice it is sufficient to select one
           vertex in the pattern and then associate it with each (reasonable)
           vertex in the subject graph.
        """
        raise NotImplementedError

    def get_new_edges(self, level):
        """Get new edges from the pattern graph for the graph search algorithm

           The level argument denotes the distance of the new edges from the
           starting vertex in the pattern graph.
        """
        raise NotImplementedError

    def check_symmetry(self, new_relations, current_match, next_match):
        """Off all symmetric ``new_relations``, only allow one"""
        return True

    def compare(self, vertex0, vertex1, subject_graph):
        """Test if ``vertex0`` and ``vertex1`` can be equal

           False positives are allowed, but the less false positives, the more
           efficient the :class:`GraphSearch` will be.
        """
        return True

    def check_next_match(self, match, new_relations, subject_graph, one_match):
        """Does this match object make sense for the current pattern

           Return False if some symmetry or other considerations are not
           satisfied. The checks in this function are typically only possible by
           considering the whole instead of looking just at a few
           vertices/edges/relations.
        """
        return True

    def complete(self, match, subject_graph):
        """Returns ``True`` if no more additional relations are required"""
        return True

    def iter_final_matches(self, match, subject_graph, one_match):
        """Just return the original match

           Derived classes can specialize here to make efficient use of
           symmetries
        """
        yield match


class CriteriaSet(object):
    """A set of criteria that can be associated with a custum pattern."""

    def __init__(self, vertex_criteria=None, edge_criteria=None, **kwargs):
        """
           Arguments:
            | ``vertex_criteria``  --  a dictionary with criteria for the
                                       vertices, ``key=vertex_index``,
                                       ``value=criterion``
            | ``edge_criteria``  --  a dictionary with criteria for the edges
                                     ``key=edge_index``, ``value=criterion``

           Any other keyword argument will be assigned as attribute to matches
           that fulfill the criteria of this set.
        """
        if vertex_criteria is None:
            self.vertex_criteria = {}
        else:
            self.vertex_criteria = vertex_criteria
        if edge_criteria is None:
            self.edge_criteria = {}
        else:
            self.edge_criteria = edge_criteria
        self.info = kwargs

    def test_match(self, match, pattern_graph, subject_graph):
        """Test if a match satisfies the criteria"""
        for vertex0, c in self.vertex_criteria.items():
            vertex1 = match.forward[vertex0]
            if not c(vertex1, subject_graph):
                return False
        for edge0_index, c in self.edge_criteria.items():
            vertex0a, vertex0b = pattern_graph.edges[edge0_index]
            edge1_index = subject_graph.edge_index[frozenset([
                match.forward[vertex0a],
                match.forward[vertex0b],
            ])]
            if not c(edge1_index, subject_graph):
                return False
        return True

# few basic example criteria

class Anything(object):
    """A criterion that always returns True"""
    def __call__(self, index, subject_graph):
        """Always returns True"""
        return True


class CritOr(object):
    """OR Operator for criteria objects"""

    def __init__(self, *criteria):
        """
           Argument:
            | ``criteria``  --  a list of criteria to apply the OR operation to.
        """
        self.criteria = criteria

    def __call__(self, index, graph):
        """Evaluates all the criteria and applies an OR opartion

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion
                             is applied
            | ``graph``  --  the graph on which the criterion is tested
        """
        for c in self.criteria:
            if c(index, graph):
                return True
        return False


class CritAnd(object):
    """AND Operator for criteria objects"""

    def __init__(self, *criteria):
        """
           Argument:
            | ``criteria``  --  a list of criteria to apply the AND operation to
        """
        self.criteria = criteria

    def __call__(self, index, graph):
        """Evaluates all the criteria and applies an AND opartion

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion
                             is applied
            | ``graph``  --  the graph on which the criterion is tested
        """
        for c in self.criteria:
            if not c(index, graph):
                return False
        return True


class CritXor(object):
    """XOR Operator for criteria objects"""

    def __init__(self, *criteria):
        """
           Argument:
            | ``criteria``  --  a list of criteria to apply the XOR operation to.
        """
        self.criteria = criteria

    def __call__(self, index, graph):
        """Evaluates all the criteria and applies a generalized XOR opartion

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion
                             is applied
            | ``graph``  --  the graph on which the criterion is tested

           when the XOR operation is applied to more than two criteria, True
           is only returned when an odd number of criteria return True.
        """
        count = 0
        for c in self.criteria:
            if c(index, graph):
                count += 1
        return (count % 2) == 1


class CritNot(object):
    """Inverion of another criterion"""

    def __init__(self, criterion):
        """
           Argument:
            | ``criterion``  --  another criterion object
        """
        self.criterion = criterion

    def __call__(self, index, graph):
        """Evaluates all the criterion and applies an inversion opartion

           Arguments:
            | ``index``  --  the index of the vertex/edge on which the criterion is
                             applied
            | ``graph``  --  the graph on which the criterion is tested
        """
        return not self.criterion(index, graph)


# pattern and match stuff


class CustomPattern(Pattern):
    """A pattern based on a given pattern graph

       The pattern graph can be complemented with additional criteria for the
       vertices and edges. Additionally the effective symmetry of the pattern
       graph can be reduced by tagging the vertices in the pattern graph with
       different labels.
    """
    def __init__(self, pattern_graph, criteria_sets=None, vertex_tags=None, start_vertex=None):
        """
        Arguments:
          | ``pattern_graph`` -- the pattern that has to be found in the subject
                                 graph.
          | ``criteria_sets`` -- Criteria sets associate additional conditions
                                 with vertices and edges, and can also introduce
                                 global match conditions.
          | ``vertex_tags`` -- vertex tags can reduce the symmetry of the
              pattern_graph pattern. An example case where this is useful:
              Consider atoms 0, 1, 2 that are bonded in this order. We want to
              compute the distance from atom 2 to the line (0, 1). In this case
              the matches (0->a, 1->b, 2->c) and (0->c, 1->b, 2->a) correspond
              to different internal coordinates. We want the graph search to
              return the two solutions. In order to do this, set
              ``vertex_tags={0:0, 1:0, 2:1}``. This means that vertex 0 and 1
              are equivalent, but that vertex 2 has a different nature. In the
              case of a bending angle, only one match like (0->a, 1->b, 2->c) is
              sufficient and we do not want to reduce the symmetry of the
              ``pattern_graph``. In this case, one should not use vertex_tags at
              all.
          | ``start_vertex``  --  The first vertex in the pattern graph that is
              linked with a vertex in the subject graph. A wise choice can
              improve the performance of a graph search. If not given, the
              central vertex is take as start_vertex
        """
        self.criteria_sets = criteria_sets
        if vertex_tags is None:
            vertex_tags = {}
        self.vertex_tags = vertex_tags
        # get the essential information from the pattern_graph:
        if start_vertex is None:
            self.start_vertex = pattern_graph.central_vertex
        else:
            self.start_vertex = start_vertex
        self._set_pattern_graph(pattern_graph)
        Pattern.__init__(self)

    def _set_pattern_graph(self, pattern_graph):
        """Initialize the pattern_graph"""
        self.pattern_graph = pattern_graph
        self.level_edges = {}
        self.level_constraints = {}
        self.duplicate_checks = set([])
        if pattern_graph is None:
            return
        if len(pattern_graph.independent_vertices) != 1:
            raise ValueError("A pattern_graph must not be a disconnected "
                             "graph.")
        # A) the levels for the incremental pattern matching
        ibfe = self.pattern_graph.iter_breadth_first_edges(self.start_vertex)
        for edge, distance, constraint in ibfe:
            if constraint:
                l = self.level_constraints.setdefault(distance-1, [])
            else:
                l = self.level_edges.setdefault(distance, [])
            l.append(edge)
        #print "level_edges", self.level_edges
        #print "level_constraints", self.level_constraints
        # B) The comparisons the should be checked when one wants to avoid
        # symmetrically duplicate pattern matches
        if self.criteria_sets is not None:
            for cycles in pattern_graph.symmetry_cycles:
                if len(cycles) > 0:
                    self.duplicate_checks.add((cycles[0][0], cycles[0][1]))


    def iter_initial_relations(self, subject_graph):
        """Iterate over all valid initial relations for a match"""
        vertex0 = self.start_vertex
        for vertex1 in range(subject_graph.num_vertices):
            if self.compare(vertex0, vertex1, subject_graph):
                yield vertex0, vertex1

    def get_new_edges(self, level):
        """Get new edges from the pattern graph for the graph search algorithm

           The level argument denotes the distance of the new edges from the
           starting vertex in the pattern graph.
        """
        return (
            self.level_edges.get(level, []),
            self.level_constraints.get(level, [])
        )

    def check_next_match(self, match, new_relations, subject_graph, one_match):
        """Check if the (onset for a) match can be a valid"""
        # only returns true for ecaxtly one set of new_relations from all the
        # ones that are symmetrically equivalent
        if not (self.criteria_sets is None or one_match):
            for check in self.duplicate_checks:
                vertex_a = new_relations.get(check[0])
                vertex_b = new_relations.get(check[1])
                if vertex_a is None and vertex_b is None:
                    continue # if this pair is completely absent in the new
                    # relations, it is either completely in the match or it
                    # is to be matched. So it is either already checked for
                    # symmetry duplicates, or it will be check in future.
                if vertex_a is None:
                    # maybe vertex_a is in the match and vertex_b is the only
                    # one in the new relations. try to get vertex_a from the
                    # match.
                    vertex_a = match.forward.get(check[0])
                    if vertex_a is None:
                        # ok, vertex_a is to be found, don't care about it right
                        # now. it will be checked in future calls.
                        continue
                elif vertex_b is None:
                    # maybe vertex_b is in the match and vertex_a is the only
                    # one in the new relations. try to get vertex_b from the
                    # match.
                    vertex_b = match.forward.get(check[1])
                    if vertex_b is None:
                        # ok, vertex_b is to be found, don't care about it right
                        # now. it will be checked in future calls.
                        continue
                if vertex_a > vertex_b:
                    # Why does this work? The answer is not so easy to explain,
                    # and certainly not easy to find. if vertex_a > vertex_b, it
                    # means that there is a symmetry operation that leads to
                    # an equivalent match where vertex_b < vertex_a. The latter
                    # match is preferred for as much pairs (vertex_a, vertex_b)
                    # as possible without rejecting all possible matches. The
                    # real difficulty is to construct a proper list of
                    # (vertex_a, vertex_b) pairs that will reject all but one
                    # matches. I conjecture that this list contains all the
                    # first two vertices from each normalized symmetry cycle of
                    # the pattern graph. We need a math guy to do the proof. -- Toon
                    return False
            return True
        return True

    def complete(self, match, subject_graph):
        """Return True of the match is complete"""
        return len(match) == self.pattern_graph.num_vertices

    def iter_final_matches(self, canonical_match, subject_graph, one_match):
        """Given a match, iterate over all related equivalent matches

           When criteria sets are defined, the iterator runs over all symmetric
           equivalent matches that fulfill one of the criteria sets. When not
           criteria sets are defined, the iterator only yields the input match.
        """
        if self.criteria_sets is None or one_match:
            yield canonical_match
        else:
            for criteria_set in self.criteria_sets:
                satisfied_match_tags = set([])
                for symmetry in self.pattern_graph.symmetries:
                    final_match = canonical_match * symmetry
                    #print final_match
                    if criteria_set.test_match(final_match, self.pattern_graph, subject_graph):
                        match_tags = tuple(
                            self.vertex_tags.get(symmetry.reverse[vertex0])
                            for vertex0
                            in range(self.pattern_graph.num_vertices)
                        )
                        if match_tags not in satisfied_match_tags:
                            final_match.__dict__.update(criteria_set.info)
                            yield final_match
                            satisfied_match_tags.add(match_tags)


class EqualMatch(Match):
    """A Match object with specialized functions for the EqualPattern"""
    def get_closed_cycles(self):
        """Return the closed cycles corresponding to this permutation

           The cycle will be normalized to facilitate the elimination of
           duplicates. The following is guaranteed:

           1) If this permutation is represented by disconnected cycles, the
              cycles will be sorted by the lowest index they contain.
           2) Each cycle starts with its lowest index. (unique starting point)
           3) Singletons are discarded. (because they are boring)
        """
        # A) construct all the cycles
        closed_cycles = []
        todo = set(self.forward.keys())
        if todo != set(self.forward.values()):
            raise GraphError("The subject and pattern graph must have the same "
                             "numbering.")
        current_vertex = None
        while len(todo) > 0:
            if current_vertex == None:
                current_vertex = todo.pop()
                current_cycle = []
            else:
                todo.discard(current_vertex)
            current_cycle.append(current_vertex)
            next_vertex = self.get_destination(current_vertex)
            if next_vertex == current_cycle[0]:
                if len(current_cycle) > 1:
                    # bring the lowest element in front
                    pivot = np.argmin(current_cycle)
                    current_cycle = current_cycle[pivot:] + \
                                    current_cycle[:pivot]
                    closed_cycles.append(current_cycle)
                current_vertex = None
            else:
                current_vertex = next_vertex
        # B) normalize the cycle representation
        closed_cycles.sort() # a normal sort is sufficient because only the
                             # first item of each cycle is considered

        # transform the structure into a tuple of tuples
        closed_cycles = tuple(tuple(cycle) for cycle in closed_cycles)
        return closed_cycles


class EqualPattern(CustomPattern):
    """Like CustomPattern, but the pattern has the same size as the subject graph"""
    sub = False
    MatchClass = EqualMatch

    def __init__(self, pattern_graph):
        """See :meth:`CustomPattern.__init__`"""
        # Don't allow criteria sets and vertex_tags. This limitation is due to
        # the compare method below. TODO: Is this a good idea?
        CustomPattern.__init__(self, pattern_graph)

    def iter_initial_relations(self, subject_graph):
        """Iterate over all valid initial relations for a match"""
        if self.pattern_graph.num_edges != subject_graph.num_edges:
            return # don't even try
        for pair in CustomPattern.iter_initial_relations(self, subject_graph):
            yield pair

    def compare(self, vertex0, vertex1, subject_graph):
        """Returns true when the two vertices are of the same kind"""
        return (
            self.pattern_graph.vertex_fingerprints[vertex0] ==
            subject_graph.vertex_fingerprints[vertex1]
        ).all()


class RingPattern(Pattern):
    """A pattern that matches strong rings up to a given size"""

    def __init__(self, max_size):
        """
           Argument:
            | ``max_size``  --  the maximum number of vertices in a ring
        """
        if max_size < 3:
            raise ValueError("Ring sizes must be at least 3.")
        self.max_size = max_size
        Pattern.__init__(self)

    def iter_initial_relations(self, subject_graph):
        """Iterate over all valid initial relations for a match"""
        vertex0 = 0
        for vertex1 in range(subject_graph.num_vertices):
            yield vertex0, vertex1

    def get_new_edges(self, level):
        """Get new edges from the pattern graph for the graph search algorithm

           The level argument denotes the distance of the new edges from the
           starting vertex in the pattern graph.
        """
        if level == 0:
            edges0 = [(0, 1), (0, 2)]
        elif level >= (self.max_size-1)//2:
            edges0 = []
        else:
            l2 = level*2
            edges0 = [(l2-1, l2+1), (l2, l2+2)]
        return edges0, []

    def check_next_match(self, match, new_relations, subject_graph, one_match):
        """Check if the (onset for a) match can be a valid (part of a) ring"""
        # avoid duplicate rings (order of traversal)
        if len(match) == 3:
            if match.forward[1] < match.forward[2]:
                #print "RingPattern.check_next_match: duplicate order", match.forward[1], match.forward[2]
                return False
        # avoid duplicate rings (starting point)
        for vertex1 in new_relations.values():
            if vertex1 < match.forward[0]:
                #print "RingPattern.check_next_match: duplicate start", vertex1, match.forward[0]
                return False
        # can this ever become a strong ring?
        for vertex1 in new_relations.values():
            paths = list(subject_graph.iter_shortest_paths(vertex1, match.forward[0]))
            if len(paths) != 1:
                #print "RingPattern.check_next_match: not strong 1"
                return False
            if len(paths[0]) != (len(match)+1)//2:
                #print "RingPattern.check_next_match: not strong 2"
                return False
        return True

    def complete(self, match, subject_graph):
        """Check the completeness of a ring match"""
        size = len(match)
        # check whether we have an odd strong ring
        if match.forward[size-1] in subject_graph.neighbors[match.forward[size-2]]:
            # we have an odd closed cycle. check if this is a strong ring
            order = list(range(0, size, 2)) + list(range(1, size-1, 2))[::-1]
            ok = True
            for i in range(len(order)//2):
                # Count the number of paths between two opposite points in the
                # ring. Since the ring has an odd number of vertices, each
                # vertex has two semi-opposite vertices.
                count = len(list(subject_graph.iter_shortest_paths(
                    match.forward[order[i]],
                    match.forward[order[(i+size//2)%size]]
                )))
                if count > 1:
                    ok = False
                    break
                count = len(list(subject_graph.iter_shortest_paths(
                    match.forward[order[i]],
                    match.forward[order[(i+size//2+1)%size]]
                )))
                if count > 1:
                    ok = False
                    break
            if ok:
                match.ring_vertices = tuple(match.forward[i] for i in order)
                #print "RingPattern.complete: found odd ring"
                return True
            #print "RingPattern.complete: no odd ring"
        # check whether we have an even strong ring
        paths = list(subject_graph.iter_shortest_paths(
            match.forward[size-1],
            match.forward[size-2]
        ))
        #print "RingPattern.complete: even paths", paths
        if (size > 3 and len(paths) == 1 and len(paths[0]) == 3) or \
           (size == 3 and len(paths) == 2 and len(paths[0]) == 3):
            path = paths[0]
            if size == 3 and path[1] == match.forward[0]:
                path = paths[1]
            # we have an even closed cycle. check if this is a strong ring
            match.add_relation(size, path[1])
            size += 1
            order = list(range(0, size, 2)) + list(range(size-1, 0, -2))
            ok = True
            for i in range(len(order)//2):
                count = len(list(subject_graph.iter_shortest_paths(
                    match.forward[order[i]],
                    match.forward[order[(i+size//2)%size]]
                )))
                if count != 2:
                    ok = False
                    break
            if ok:
                # also check if this does not violate the requirement for a
                # unique origin:
                if match.forward[size-1] < match.forward[0]:
                    ok = False
            if not ok:
                vertex1 = match.forward[size-1]
                del match.forward[size-1]
                del match.reverse[vertex1]
                size -= 1
                #print "RingPattern.complete: no even ring"
            else:
                match.ring_vertices = tuple(match.forward[i] for i in order)
                #print "RingPattern.complete: found even ring"
            return ok
        #print "RingPattern.complete: not at all"
        return False


class GraphSearch(object):
    """An algorithm that searches for all matches of a pattern in a graph

       Usage:

         >>> gs = GraphSearch(pattern)
         >>> for match in gs(graph):
         ...     print match.forward
    """

    def __init__(self, pattern, debug=False):
        """
           Arguments:
            | ``pattern``  --  A Pattern instance, describing the pattern to
                               look for
            | ``debug``  --  When true, debugging info is printed on screen
                             [default=False]
        """
        self.pattern = pattern
        self.debug = debug

    def __call__(self, subject_graph, one_match=False):
        """Iterator over all matches of self.pattern in the given graph.

           Arguments:
            | subject_graph  --  The subject_graph in which the matches
                                 according to self.pattern have to be found.
            | one_match --  If True, only one match will be returned. This
                            allows certain optimizations.
        """
        # Matches are grown iteratively.
        for vertex0, vertex1 in self.pattern.iter_initial_relations(subject_graph):
            init_match = self.pattern.MatchClass.from_first_relation(vertex0, vertex1)
            # init_match cotains only one source -> dest relation. starting from
            # this initial match, the function iter_matches extends the match
            # in all possible ways and yields the completed matches
            for canonical_match in self._iter_matches(init_match, subject_graph, one_match):
                # Some patterns my exclude symmetrically equivalent matches as
                # to aviod dupplicates. with such a 'canonical' solution,
                # the pattern is allowed to generate just those symmatrical
                # duplicates of interest.
                ifm = self.pattern.iter_final_matches(canonical_match, subject_graph, one_match)
                for final_match in ifm:
                    self.print_debug("final_match: %s" % final_match)
                    yield final_match
                    if one_match: return

    def print_debug(self, text, indent=0):
        """Only prints debug info on screen when self.debug == True."""
        if self.debug:
            if indent > 0:
                print(" "*self.debug, text)
            self.debug += indent
            if indent <= 0:
                print(" "*self.debug, text)

    def _iter_candidate_groups(self, init_match, edges0, edges1):
        """Divide the edges into groups"""
        # collect all end vertices0 and end vertices1 that belong to the same
        # group.
        sources = {}
        for start_vertex0, end_vertex0 in edges0:
            l = sources.setdefault(start_vertex0, [])
            l.append(end_vertex0)
        dests = {}
        for start_vertex1, end_vertex1 in edges1:
            start_vertex0 = init_match.reverse[start_vertex1]
            l = dests.setdefault(start_vertex0, [])
            l.append(end_vertex1)
        for start_vertex0, end_vertices0 in sources.items():
            end_vertices1 = dests.get(start_vertex0, [])
            yield end_vertices0, end_vertices1


    def _iter_new_relations(self, init_match, subject_graph, edges0, constraints0, edges1):
        """Given an onset for a match, iterate over all possible new key-value pairs"""
        # Count the number of unique edges0[i][1] values. This is also
        # the number of new relations.
        num_new_relations = len(set(j for i, j in edges0))

        def combine_small(relations, num):
            """iterate over all compatible combinations within one set of relations"""
            if len(relations) == 0:
                return
            for i, pivot in enumerate(relations):
                if num == 1:
                    yield (pivot, )
                else:
                    compatible_relations = list(
                        item for item in relations[:i]
                        if pivot[0]!=item[0] and pivot[1]!=item[1]
                    )
                    for tail in combine_small(compatible_relations, num-1):
                        yield (pivot, ) + tail

        # generate candidate relations
        candidate_relations = []
        icg = self._iter_candidate_groups(init_match, edges0, edges1)
        for end_vertices0, end_vertices1 in icg:
            if len(end_vertices0) > len(end_vertices1):
                return # this can never work, the subject graph is 'too small'
            elif not self.pattern.sub and \
                 len(end_vertices0) != len(end_vertices1):
                return # an exact match is sought, this can never work
            l = []
            for end_vertex0 in end_vertices0:
                for end_vertex1 in end_vertices1:
                    if self.pattern.compare(end_vertex0, end_vertex1, subject_graph):
                        l.append((end_vertex0, end_vertex1))
            # len(end_vertices0) = the total number of relations that must be
            # made in this group
            if len(l) > 0:
                # turn l into a list of sets of internally compatible candidate
                # relations in this group
                l = list(combine_small(l, len(end_vertices0)))
                candidate_relations.append(l)
        if len(candidate_relations) == 0:
            return
        self.print_debug("candidate_relations: %s" % candidate_relations)

        def combine_big(pos=0):
            """Iterate over all possible sets of relations"""
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
            self.print_debug("new_relations: %s" % (new_relations, ))
            # check the total number of new relations
            if len(new_relations) != num_new_relations:
                continue
            # check sanity of relations
            forward = dict(new_relations)
            if len(forward) != num_new_relations:
                continue
            reverse = dict((j, i) for i, j in new_relations)
            if len(reverse) != num_new_relations:
                continue
            # check the constraints
            for a0, b0 in constraints0:
                if forward[a0] not in subject_graph.neighbors[forward[b0]]:
                    forward = None
                    break
            if forward is None:
                continue
            yield forward

    def _iter_matches(self, input_match, subject_graph, one_match, level=0):
        """Given an onset for a match, iterate over all completions of that match

           This iterator works recursively. At each level the match is extended
           with a new set of relations based on vertices in the pattern graph
           that are at a distances 'level' from the starting vertex
        """
        self.print_debug("ENTERING _ITER_MATCHES", 1)
        self.print_debug("input_match: %s" % input_match)
        # A) collect the new edges in the pattern graph and the subject graph
        # to extend the match.
        #
        # Note that the edges are ordered. edge[0] is always in the match.
        # edge[1] is never in the match. The constraints contain information
        # about the end points of edges0. It is a list of two-tuples where
        # (a, b) means that a and b must be connected.
        #
        # Second note: suffix 0 indicates the pattern graph and suffix 1
        # is used for the subject graph.
        edges0, constraints0 = self.pattern.get_new_edges(level)
        edges1 = input_match.get_new_edges(subject_graph)
        self.print_debug("edges0: %s" % edges0)
        self.print_debug("constraints0: %s" % constraints0)
        self.print_debug("edges1: %s" % edges1)

        # B) iterate over the sets of new relations: [(vertex0[i], vertex1[j]),
        # ...] that contain all endpoints of edges0, that satisfy the
        # constraints0 and where (vertex0[i], vertex1[j]) only occurs if these
        # are end points of a edge0 and edge1 whose starting points are already
        # in init_match. These conditions are implemented in an iterator as to
        # separate concerns. This iterator also calls the routines that check
        # whether vertex1[j] also satisfies additional conditions inherent
        # vertex0[i].
        inr = self._iter_new_relations(input_match, subject_graph, edges0,
                                       constraints0, edges1)
        for new_relations in inr:
            # for each set of new_relations, construct a next_match and recurse
            next_match = input_match.copy_with_new_relations(new_relations)
            if not self.pattern.check_next_match(next_match, new_relations, subject_graph, one_match):
                continue
            if self.pattern.complete(next_match, subject_graph):
                yield next_match
            else:
                for match in self._iter_matches(next_match, subject_graph, one_match, level+1):
                    yield match
        self.print_debug("LEAVING_ITER_MATCHES", -1)
