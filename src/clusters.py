# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


__all__ = ["Cluster", "ClusterFactoryError", "ClusterFactory"]


class Cluster(object):
    def __init__(self, members=None):
        if members is None:
            self.members = []
        else:
            self.members = members

    def clear(self):
        self.members = []

    def add_member(self, member):
        self.members.append(member)

    def add_cluster(self, cluster):
        self.members.extend(cluster.members)
        cluster.clear()


class ClusterFactoryError(Exception):
    pass


class ClusterFactory(object):
    def __init__(self, ClusterClass=Cluster):
        self.clusters = {}
        self.ClusterClass = ClusterClass

    def add_cluster(self, master):
        slaves = set([]) # set of clusters that is going to be merged in the master
        new = [] # members that are not yet in other clusters
        for member in master.members:
            cluster = self.clusters.get(member)
            if cluster is None:
                #print "solitaire", member
                new.append(member)
            else:
                #print "in slave", member
                slaves.add(cluster)
            #else:
            #    #print "in master", member

        for slave in slaves:
            for member in slave.members:
                self.clusters[member] = master
            master.add_cluster(slave)
        for member in new:
            self.clusters[member] = master

    def add_members(self, members):
        master = None # this will become the common cluster of all related members
        slaves = set([]) # set of clusters that is going to be merged in the master
        solitaire = [] # members that are not yet part of a cluster
        for member in members:
            cluster = self.clusters.get(member)
            if cluster is None:
                #print "solitaire", member
                solitaire.append(member)
            elif master is None:
                #print "starting master", member
                master = cluster
            elif master != cluster:
                #print "in slave", member
                slaves.add(cluster)
            #else:
                #print "in master", member

        if master is None:
            master = self.ClusterClass()
        else:
            for slave in slaves:
                for member in slave.members:
                    self.clusters[member] = master
                master.add_cluster(slave)
        for member in solitaire:
            self.clusters[member] = master
            master.add_member(member)

    def get_clusters(self):
        return set(cluster for cluster in self.clusters.itervalues())


