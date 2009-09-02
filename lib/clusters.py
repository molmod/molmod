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
"""Basic cluster analysis tool

   Given a mixed set of related and unrelated data pionts, it is often interesting
   to extract clusters of related items. The basic workflow is as follows:

   >>> cf = ClusterFactory()

   >>> some loop:
   ...     cf.add_related(some, related, items)

   >>> for cluster in cf.iter_clusters():
   ...     print cluster
"""

__all__ = ["ClusterFactory"]


class ClusterFactory(object):
    """A very basic cluster algorithm"""

    def __init__(self):
        """Initialize a ClusterFactory"""
        # mapping: item -> cluster. Each cluster is a tuple of related items.
        self.lookup = {}

    def add_related(self, *group):
        """Add related items

           When two groups of related items share one or more common members,
           they will be merged into one cluster.
        """
        master = None # this will become the common cluster of all related group
        slaves = set([]) # set of clusters that are going to be merged in the master
        solitaire = set([]) # set of new items that are not yet part of a cluster
        for new in group:
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
            master = []
        else:
            master = list(master)

        for slave in slaves:
            master.extend(slave)
        master.extend(solitaire)
        master = tuple(master)

        for item in master:
            self.lookup[item] = master

    def get_clusters(self):
        """Returns a set with the clusters"""
        return set(self.lookup.itervalues())





