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
"""Basic cluster analysis tool

   Given a mixed set of related and unrelated data pionts, it is often interesting
   to extract clusters of related items. The basic workflow is as follows::

       cf = ClusterFactory()
       while foo:
           cf.add_related(some, related, items)
       for cluster in cf.iter_clusters():
           print cluster
"""


__all__ = ["Cluster", "RuleCluster", "ClusterFactory"]


class Cluster(object):
    """A set of related items

       This is the most elementary implementation of a cluster. In practice
       on is often interested in extending the functionality of a cluster.
    """
    def __init__(self, items):
        """
           Argument:
            | ``items``  --  the items that belong in this cluster
        """
        self.items = set(items)

    def add_item(self, item):
        """Add an item to a cluster"""
        self.items.add(item)

    def update(self, other):
        """Merge another cluster into this cluster"""
        self.items |= other.items


class RuleCluster(Cluster):
    """Clusters based on rules

       This is a typical derived Cluster class where the relation between the
       items is one or more rules, which one would like to know at the end of
       the clustering algorithm.

       An example application is the shake algorithm where it is beneficial
       to group constraints that share certain degrees of freedom into a cluster
       of equations.
    """
    def __init__(self, items, rules=None):
        """
           Argument:
            | ``items``  --  the items that belong in this cluster

           Optional argument:
            | ``rules``  --  a list of rules that binds the items
        """
        Cluster.__init__(self, items)
        if rules is None:
            self.rules = []
        else:
            self.rules = rules

    def update(self, other):
        """Extend the current cluster with data from another cluster"""
        Cluster.update(self, other)
        self.rules.extend(other.rules)


class ClusterFactory(object):
    """A very basic cluster algorithm"""

    def __init__(self, cls=Cluster):
        """
           Optinional argument:
            | ``cls``  --  A class to construct new cluster objects
                           [default=Cluster]
        """
        self.cls = cls
        # mapping: item -> cluster. Each cluster is a tuple of related items.
        self.lookup = {}

    def add_related(self, *objects):
        """Add related items

           The arguments can be individual items or cluster objects containing
           several items.

           When two groups of related items share one or more common members,
           they will be merged into one cluster.
        """
        master = None # this will become the common cluster of all related items
        slaves = set([]) # set of clusters that are going to be merged in the master
        solitaire = set([]) # set of new items that are not yet part of a cluster
        for new in objects:
            if isinstance(new, self.cls):
                if master is None:
                    master = new
                else:
                    slaves.add(new)
                for item in new.items:
                    existing = self.lookup.get(item)
                    if existing is not None:
                        slaves.add(existing)
            else:
                cluster = self.lookup.get(new)
                if cluster is None:
                    #print "solitaire", new
                    solitaire.add(new)
                elif master is None:
                    #print "starting master", new
                    master = cluster
                elif master != cluster:
                    #print "in slave", new
                    slaves.add(cluster)
                #else:
                    ##nothing to do
                    #print "new in master", new

        if master is None:
            master = self.cls([])

        for slave in slaves:
            master.update(slave)
        for item in solitaire:
            master.add_item(item)

        for item in master.items:
            self.lookup[item] = master

    def get_clusters(self):
        """Returns a set with the clusters"""
        return set(self.lookup.values())
