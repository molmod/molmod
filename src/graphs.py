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

import copy

__all__ = ["OneToOne", "Permutation", "Graph",
           "SymmetricGraph", "Criterium",
           "MatchFilter", "MatchFilterParameterized", "MatchFilterMolecular"]


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
        """
        Extend the bijection with one new relation. Returns true of the
        new relation is accepted or when the new relation is already made.
        """
        if self.in_sources(source):
            if self.forward[source] == destination:
                return True
            else:
                return False
        elif self.in_destinations(destination):
            return False
        else:
            self.forward[source] = destination
            self.backward[destination] = source
            return True
        
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
    representations: pairs and neigbours.
    """
    
    def __init__(self, pairs):
        """
        Initialize a Graph object.
        
        Arguments:
        pairs -- [frozenset([node1, node2]), ...]
        
        During initialization, also the neighbour-representation will be created
        neighbourlist = {node1: [neighbour_node1, ...], ...}
        """
        #print "="*50
        #print "="*50
        self.pairs = pairs # a list of 2-frozenset connecting nodes
        #print self.pairs
        self.init_neighbours()

    def init_neighbours(self):
        """Generate a neigbours-representation of the graph."""
        self.neighbours = {}
        
        def add_relation(first, second):    
            if first in self.neighbours:
                self.neighbours[first].append(second)
            else:
                self.neighbours[first] = [second]
        
        
        for pair in self.pairs:
            a, b = pair
            add_relation(a, b)
            add_relation(b, a)
        
        #print self.neighbours

        
class SymmetricGraph(Graph):
    def __init__(self, pairs, initiator=None):
        """
        Initialize the graph and find all symmetries in the graph.
        
        Arguments:
        pairs -- see Graph
        initiator -- a node with as less as possible symmetric equivalents
        
        If possible, make sure that the initiator is a node with as less as
        possible equivalent nodes in the graph. The equivalent nodes of a given
        node are all the nodes on which the given node may me mapped by a
        symmetry transformation.
        """
        Graph.__init__(self, pairs)
        if initiator == None:
            self.initiator = list(self.pairs[0])[0]
        else:
            self.initiator = initiator
        self.init_symmetries()
        self.init_initiator_cycle()
                
    def init_symmetries(self):
        """
        Analyze the symmetries of this graph.
        
        This method uses yield_matches to find all the symmetries in this
        graph. Such symmetry information is essential for further analysis
        of the graph.
        
        self.symmetries is a dictionary (cycle_lists -> permutations). Each
        cycle/permutation describes how the nodes can be interchanged without
        modifying the graph in any aspect.
        
        self.cycles is a list of all cycles that occur in the graph. A cycle
        is defined as undividable subpermutation of a symmetry. It is also a
        list of equivalent nodes.
        """
        self.symmetries = {}
        self.cycles = []
        
        def add_symmetry(symmetry):
            cycles = symmetry.get_closed_cycles()
            if not cycles in self.symmetries:
                self.symmetries[cycles] = symmetry
                for cycle in cycles:
                    if cycle not in self.cycles:
                        self.cycles.append(cycle)
                return True
            else:
                return False
                
        def find_symmetries(source, destination):
            #print "SYMMETRIES %s ~ %s" % (source, destination)
            for symmetry in self.yield_matches(source, destination, self.neighbours, Permutation(), allow_more=False, no_duplicates=False):
                if add_symmetry(symmetry):
                    add_symmetry(symmetry.inverse())
        
        # The unity transformation will certainly be found by this call
        find_symmetries(self.initiator, self.initiator)
        # Non unity transformations will certainly be found in this loop
        for source in self.neighbours:
            for destination in self.neighbours:
                if source >= destination: continue
                find_symmetries(source, destination)
         
                    
        #print "-- number of symmetries: %i" % (len(self.symmetries))
        #for cycles in self.symmetries.iterkeys():
        #    print cycles
        #print "-- number of cycles: %i" % (len(self.cycles))
        #for cycle in self.cycles:
        #    print cycle

    def init_initiator_cycle(self):
        """Find a largest cycle that contains the initiator."""
        self.initiator_cycle = (self.initiator, )
        for cycle in self.cycles:
            if self.initiator in cycle and len(cycle) > len(self.initiator_cycle):
                self.initiator_cycle = cycle
        #print "initiator_cycle", self.initiator_cycle

    def len_cycle_prefix(self, sequence):
        """Returns the largest cycle at the beginning of the given sequence."""
        while (sequence not in self.cycles) and len(sequence) > 0:
            sequence = sequence[:-1]
        return len(sequence)
        
    def yield_matches(self, node, thing, thing_neighbours, match, allow_more=True, no_duplicates=True, spaces=""):
        """
        Generate all matching subgraphs congruent to the graph described by
        thing_neigbours where node and thing are corresponding items between
        this graph and the subgraph.
        
        This code is the core algorithm of the Graph class. It uses extensively
        the yield keyword. If you plan to analyze this algorithm, make sure
        you have a very good understanding of the generator concept in the
        python language. Note that this algorithm is recursive.
        
        Arguments:
        node -- the node that corresponds to thing in the subgraph
        thing -- the item in the subgraph that corresponds to node
        thing_neighbours -- a neighbour representation of the subgraph
        match -- a permutation that describes the partial symmetry at a
                 recursion level. user should give an empty permutation as
                 input
        allow_more -- wether a node in this graph may have more neigbours than
                      the corresponding thing in the subgraph. For symmetry 
                      problems this must be set to False for efficiency reasons.
                      For subgraph problems this may be True, depending on your
                      preferences.
        no_duplicates -- lets the algorithm use symmetry information to a priori
                         avoid duplicate mathces.
        """
        ##print spaces + "BEGIN ONE -- add %s -> %s to %s" % (node, thing, match)
        
        if allow_more and len(self.neighbours[node]) > len(thing_neighbours[thing]):
            ##print spaces + "ONE -- graph node containse more neighbours than thing"
            return
            
        if not allow_more and len(self.neighbours[node]) != len(thing_neighbours[thing]):
            ##print spaces + "ONE -- number of neighbours differs"
            return
            
        if not match.add_relation(node, thing):
            ##print spaces + "ONE -- given relation conflicts with existing match"
            return
            
        if len(match) == len(self.neighbours):
            ##print spaces + "ONE -- match is complete"
            yield match
            return
        
        def worth_trying(relation):
            """
            Wether a candidate mapping between a node and a thing is usefull.
            
            Such a priori exclussions reduce the number of loops below. The 
            following relations are excluded:
                - relations that certainly conflict with the existing match: a&b
                - relations that conflict with the allow_more parameter: c&d
            Note that a relation that is already part of the match, will be
            accepted.
            """
            if match.in_sources(relation[0]): #a
                if match.get_destination(relation[0]) == relation[1]:
                    return True
                else:
                    return False
            elif match.in_destinations(relation[1]): #b
                return False
            elif allow_more and len(self.neighbours[relation[0]]) > len(thing_neighbours[relation[1]]): #c
                return False
            elif not allow_more and len(self.neighbours[relation[0]]) != len(thing_neighbours[relation[1]]): #d
                return False
            else:
                return True

        ##print spaces + "ONE -- possible relations  %s -> %s" % (self.neighbours[node], thing_neighbours[thing])
                
        if no_duplicates:
            def relations_without_symmetric_images(list1, list2, worth_trying=None, num_ordered=0):
                """
                Return all relation sets that do not cause duplicate matches to
                be generated.
                
                Arguments:
                list1 -- list of source items
                list2 -- list of destination items
                worth_trying -- a function that takes a relation (tuple) as
                                argument and returns False if the relation is
                                not usefull.
                num_ordered -- help, I don't know any more what this means!!!
                               Someone please...
                """
                if len(list1) == 0 or len(list2) == 0:
                    yield []
                else:
                    if num_ordered == 0:
                        num_ordered = self.len_cycle_prefix(tuple(list1))
                        
                    if num_ordered > 0:
                        length = len(list2)-len(list1)+1
                        new_num_ordered = num_ordered - 1
                    else:
                        length = len(list2)
                        new_num_ordered = 0
                        
                    for index2 in xrange(length):
                        pair = (list1[0], list2[index2])
                        if worth_trying == None or worth_trying(pair):
                            sublist2 = list2[index2+1:]
                            if num_ordered == 0:
                                sublist2 = list2[:index2] + sublist2
                            for relation_set in relations_without_symmetric_images(list1[1:], sublist2, worth_trying, new_num_ordered):
                                yield [pair] + relation_set                                     
                
            relations = relations_without_symmetric_images
        else:
            relations = all_relations
            
        for relation_set in relations(self.neighbours[node], thing_neighbours[thing], worth_trying):
            # Here we will try to extend the given match with new relations.
            # These extend matches will feed recursively to this algorithm.
            if len(relation_set) == 0:
                continue
            ##print spaces + "BEGIN SET -- add %s to %s" % (relation_set, match)
            former_matches = [copy.copy(match)]
            for node, thing in relation_set:
                new_matches = []
                for former_match in former_matches:
                    if former_match.in_sources(node) and former_match.get_destination(node) == thing:
                        new_matches.append(former_match)
                    # suppose this is more efficient:
                    # elif not (former_match.in_sources(node) or former_match.in_destination(thing)):
                    else:
                        new_matches.extend(list(self.yield_matches(node, thing, thing_neighbours, former_match, allow_more, no_duplicates, spaces=spaces+"  ")))
                former_matches = new_matches
                
            for new_match in new_matches:
                ##print spaces + "SET -- yield %s" % new_match
                yield new_match
  
            ##print spaces + "END SET"
            
        ##print spaces + "END ONE"
     
    def lowest_initiator(self, match):
        """
        Returns true if the thing associated with the initiator is the lowest
        of all things associated with the largest cycle that contains an
        initiator.
        """
        return match.get_destination(self.initiator) == min([match.get_destination(node) for node in self.initiator_cycle])

    def yield_matching_subgraphs(self, graph):
        """
        Returns all macthings subgraphs in terms bijections between nodes and
        things.
        
        This function explains why the 'initiator' should have as less as
        possible equivalent nodes. The cost of the algorithm is proportional
        to this number of equivalent nodes (including the initiator).
        """
        for thing in graph.neighbours:
            for match in self.yield_matches(self.initiator, thing, graph.neighbours, OneToOne()):
                if self.lowest_initiator(match):
                    #print match
                    yield match


class Criterium(object):
    """A base class for all thing_ and relation_criteria."""
    
    def __init__(self, *parameters):
        self.parameters = parameters
        
    def get_tag(self):
        """Return a tag that uniquely identifies the behaviour of a criterium."""
        if "parameters" not in self.__dict__:
            print self, self.__class__
        return (self.__class__, self.parameters)


class MatchFilter(object):
    """
    MatchFilter instances help analyzing and selecting matched subgraphs
    generated by a Graph instance. This is an abstract base class.
    """
    
    def __init__(self, subgraph, calculation_tags):
        """
        Initialize a MatchFilter instance
        
        Arguments:
        subgraph -- the subgraph to be searched for
        calculation_tags -- {node: calculation_tag}. This argument tells the the 
                            MatchFilter what the different types of nodes are in 
                            the subgraph.
        
        Example usage: You want to list all angles 1-0-2 of the sp3 atoms in a
        molecule.
        MatchFilter(
            subgraph = SymmetricGraph([(0, 1), (0, 2), (0, 3), (0, 4)], 0)
            calculation_tags = {
                0: 0,
                1: 1, 2: 1,
                3: 2, 4: 2
            }
        )
        The calculation_tags indicate that the calculation of the angle is
        invariant under the exchange of atom 1 and 2, and under the exchange of
        atom 3 and 4.
        """
        self.subgraph = subgraph
        self.calculation_tags = calculation_tags
        
        # create a symmetry subgroup that only contains symmetry transformations
        # that leave invariant:
        # - the topology of the subgraph
        # - the calculation of the variable based on the matched subgraph
        #   (this is user defined by calculation_tags
        # - the thing and relation criteria (has to be implemented in
        #   subclasses)
        self.lower_symmetries = dict(
            (cycle_representation, symmetry) 
            for cycle_representation, symmetry 
            in self.subgraph.symmetries.iteritems() 
            if self.invariant_tags(symmetry)
        )
        # remove unity
        del self.lower_symmetries[tuple()]
        #print self.lower_symmetries
        
    def invariant_tags(self, symmetry):
        for node, tag in self.calculation_tags.iteritems():
            transformed_tag = self.calculation_tags.get(symmetry.get_destination(node))
            if transformed_tag != None and transformed_tag != tag:
                return False
        return True
            
    def check_thing(self, node, thing):
        raise NotImplementedError

    def check_relation(self, nodes, things):
        raise NotImplementedError
        
    def parse(self, match):
        """
        Return transformed matches that obey the criteria defined in
        check_thing and check_relation and avoid symmetric duplicates.
        """
        def check_match(match):
            """
            For the given (transformed) match return wether it meets the
            criteria.
            """
            for node, neighbours in self.subgraph.neighbours.iteritems():
                if not self.check_thing(node, transformed_match.get_destination(node)):
                    return False
                for neighbour in neighbours:
                    if not self.check_relation(
                        frozenset([node, neighbour]), 
                        frozenset([
                            transformed_match.get_destination(node), 
                            transformed_match.get_destination(neighbour)
                        ])):
                        return False
            return True

        symmetries = self.subgraph.symmetries.copy()
        while len(symmetries) > 0:
            key, symmetry = symmetries.popitem()
            transformed_match = match * symmetry
            # remove all elements of the coset(symmetry, lower_symmetries) in
            # the group symmetries.
            for closed_cycles, lower_symmetry in self.lower_symmetries.iteritems():
                del symmetries[(symmetry * lower_symmetry).get_closed_cycles()]
            if check_match(transformed_match):
                yield transformed_match


class MatchFilterParameterized(MatchFilter):
    """
    This MatchFilter uses dictionaries that map each node and relation of the
    subgraph to a criterion function.
    """
    
    def __init__(self, subgraph, calculation_tags, thing_criteria={}, relation_criteria={}, filter_tags=True):
        """
        Initialize a MatchFilterParameterized
        
        Arguments:
        thing_criteria -- {node: thing_criterion(thing), ...}
        relation_criteria -- {frozenset([node1, node2]: 
                              relation_criterion(thing1, thing2), ...}
        """
        self.thing_criteria = thing_criteria
        self.relation_criteria = relation_criteria
        self.filter_tags = filter_tags
        MatchFilter.__init__(self, subgraph, calculation_tags)
        
    def invariant_tags(self, symmetry):
        if not MatchFilter.invariant_tags(self, symmetry):
            return False
        if not self.filter_tags:
            return True
            
        nodes = set(pair[0] for pair in self.subgraph.pairs) | \
                set(pair[1] for pair in self.subgraph.pairs)
            
        def node_tag(node):
            thing_criterium = self.thing_criteria.get(node)
            if thing_criterium == None:
                return None
            else:
                return thing_criterium.get_tag()
                
        for node in nodes:
            tag = node_tag(node)
            transformed_tag = node_tag(symmetry.get_destination(node))
            if tag != None and transformed_tag != None and transformed_tag != tag:
                return False

        def pair_tag(pair):
            relation_criterium = self.relation_criteria.get(pair)
            if relation_criterium == None:
                return None
            else:
                return relation_criterium.get_tag()
                
        for pair in self.subgraph.pairs:
            tag = node_tag(pair)
            node1, node2 = pair
            transformed_pair = frozenset([symmetry.get_destination(node1), symmetry.get_destination(node2)])
            transformed_tag = node_tag(transformed_pair)
            if tag != None and transformed_tag != None and transformed_tag != tag:
                return False
                
        return True
            
        
    def check_thing(self, node, thing):
        thing_criterium = self.thing_criteria.get(node)
        if thing_criterium == None:
            return True
        else:
            return thing_criterium(thing)
        
    def check_relation(self, nodes, things):
        relation_criterium = self.relation_criteria.get(nodes)
        if relation_criterium == None:
            return True
        else:
            return relation_criterium(things)
