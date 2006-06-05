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
"""
Binning is a usefull technique for efficiently calculating all distances between
a set of points, when you are only interested in the distances below a given
cutoff. The algorithm consists of two major steps:
1) Divide the given set of points into bins on a regular grid in space.
2) Calculate the distances (or other usefull things) between points in 
   neighbouring bins.
"""


from molmod.moldata import periodic

import math, numpy, copy

__all__ = ["PositionedObject", "SparseBinnedObjects",
           "AnalyseNeighbouringObjects", "IntraAnalyseNeighbouringObjects",
           "InterAnalyseNeighbouringObjects"]


class PositionedObject(object):
    """
    PositionedObject instances are used to feed a SparseBinnedObjects instance
    with information about the points that should be binned. Throug the use of
    a PositionedObject, each point is associated with a reference.
    """
    
    def __init__(self, reference, point):
        """
        Initialize a PositionedObject instance.
        
        Arguments:
        reference -- A user defined reference that is associated with the point.
        point -- A numpy array with shape (3,)
        """
        self.reference = reference
        self.point = point


class SparseBinnedObjects(object):
    """
    A SparseBinnedObjects instance divides 3D space into a sparse grid.
    
    Each cell in the grid is called a bin. Each bin can contain a set of
    Positioned objects. This implementation works with sparse bins: A bin is
    only created in memory when an object is encountered that belongs in that
    bin.
    
    All bins are uniquely defined by their indices i,j,k as defined in __init__.
    """
    def __init__(self, yield_positioned_objects, gridsize=1):
        """
        Initialize a BinnedObjects instance.
        
        Arguments:
        yield_positioned_objects -- A generator that iterates over all
                                    PositionedObject instances that have to be
                                    binned.
        gridsize -- Defines the size of the bins.
        """
        self.gridsize = gridsize
        self.reciproke = 1/gridsize

        self.bins = {}
            
        for positioned_object in yield_positioned_objects():
            i, j, k = numpy.floor(positioned_object.point*self.reciproke).astype(int)
            bin_x = self.bins.get(i)
            if bin_x != None:
                bin_y = bin_x.get(j)
                if bin_y != None:
                    bin_z = bin_y.get(k)
                    if bin_z == None:
                        bin_z = set()
                        bin_y[k] = bin_z
                else:
                    bin_z = set()
                    bin_x[j] = {k: bin_z}
            else:
                bin_z = set()
                self.bins[i] = {j: {k: bin_z}}
            bin_z.add(positioned_object)

    def get_bin(self, i, j, k):
        """Return the bin at indices i,j,k."""
        bin_x = self.bins.get(i)
        if bin_x != None:
            bin_y = bin_x.get(j)
            if bin_y != None:
                return bin_y.get(k)


    def yield_positioned_objects_around(self, r):
        """
        Iterate over all objects in the bins that surround the bin that
        contains vector r.
        """
        center_indices = numpy.floor(r*self.reciproke).astype(int)
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    bin = self.get_bin(
                        center_indices[0] + i, 
                        center_indices[1] + j, 
                        center_indices[2] + k
                    )
                    if bin != None:
                        for positioned_object in bin:
                            yield positioned_object


class AnalyseNeighbouringObjects(object):
    """
    AnalyseNeighbouringObjects is the base class for 'comparing' points between
    neigbouring bins.
    """
    corners = (numpy.array([0, 0, 0], int),
               numpy.array([1, 0, 0], int),
               numpy.array([0, 1, 0], int),
               numpy.array([0, 0, 1], int),
               numpy.array([0, 1, 1], int),
               numpy.array([1, 0, 1], int),
               numpy.array([1, 1, 0], int),
               numpy.array([1, 1, 1], int))

    def __init__(self, compare_function):
        """
        Intialize a AnalyseNeighbouringObjects instance.
        
        compare_function -- for each point pair that live in neighbouring bins,
                            this function is called. It should take four
                            parameters:
                            reference1, reference2, point1, point2
        """
        self.compare_function = compare_function
        # All these parameters have to be defined by the base class
        self.gridsize = None
        self.reciproke = None
        self.compare_indices = None
        self.binned_objects1 = None
        self.binned_objects2 = None

    def __call__(self):
        for i, bin1_x in self.binned_objects1.bins.iteritems():
            for j, bin1_y in bin1_x.iteritems():
                for k, center_bin in bin1_y.iteritems():
                    neighbour_bins = []
                    for di,dj,dk in self.compare_indices:
                        bin2_z = self.binned_objects2.get_bin(i+di, j+dj, k+dk)
                        if bin2_z != None:
                            neighbour_bins.append(bin2_z)
                    
                    for neighbour_bin in neighbour_bins:
                        for result in self.compare_bins(center_bin, neighbour_bin):
                            yield result

    def compare_bins(self, bin1, bin2):
        raise NotImplementedError


class IntraAnalyseNeighbouringObjects(AnalyseNeighbouringObjects):
    """
    IntraAnalyseNeighbouringObjects instances compare all points within one
    molecule.
    """
    def __init__(self, binned_objects, compare_function):
        AnalyseNeighbouringObjects.__init__(self, compare_function)
        self.gridsize = binned_objects.gridsize
        self.reciproke = binned_objects.reciproke
        self.compare_indices = [(0, 0, 0), (1, 1, 1), 
                                (1, 0, 0), (0, 1, 0), (0, 0, 1),
                                (0, 1, 1), (1, 0, 1), (1, 1, 0), 
                                (0, 1, -1), (-1, 0, 1), (1, -1, 0),
                                (1, 1, -1), (1, -1, -1), (1, -1, 1)]
        self.binned_objects1 = binned_objects
        self.binned_objects2 = binned_objects

    def compare_bins(self, bin1, bin2):
        for positioned1 in bin1:
            for positioned2 in bin2:
                if bin1 == bin2 and positioned1.reference >= positioned2.reference:
                    continue
                result = self.compare_function(
                    positioned1.reference, 
                    positioned2.reference, 
                    positioned1.point, 
                    positioned2.point
                )
                if result != None:
                    yield (frozenset((positioned1.reference, positioned2.reference)), result)


class InterAnalyseNeighbouringObjects(AnalyseNeighbouringObjects):
    """
    InterAnalyseNeighbouringObjects instances compare 'all' points between two
    molecules.
    """
    def __init__(self, binned_objects1, binned_objects2, compare_function):
        AnalyseNeighbouringObjects.__init__(self, compare_function)
        assert binned_objects1.gridsize==binned_objects2.gridsize
        self.gridsize = binned_objects1.gridsize
        self.reciproke = binned_objects1.reciproke
        self.compare_indices = []
        positions = [-1, 0, 1]
        for x in positions:
            for y in positions:
                for z in positions:
                    self.compare_indices.append((x, y, z))
        self.binned_objects1 = binned_objects1
        self.binned_objects2 = binned_objects2

    def compare_bins(self, bin1, bin2):
        for positioned1 in bin1:
            for positioned2 in bin2:
                result = self.compare_function(
                    positioned1.reference, 
                    positioned2.reference, 
                    positioned1.point, 
                    positioned2.point)
                if result != None:
                    yield (frozenset((positioned1.reference, positioned2.reference)), result)


