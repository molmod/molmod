# Zeobuilder is an extensible GUI-toolkit for molecular model construction.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of Zeobuilder.
# 
# Zeobuilder is free software; you can redistribute it and/or
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


__all__ = ["Cluster", "ClusterFactoryError", "ClusterFactory"]


class Cluster(object):
    def __init__(self):
        self.clear()
    
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
        
    def add_group(self, members):
        if len(members) < 2:
            raise ClusterFactoryError("At least groups of two members needed.")
        
        master = None # first encountered existing cluster 
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
            #    #print "in master", member
        
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

