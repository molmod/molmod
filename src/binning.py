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
a set of vectors, when you are only interested in the distances below a given
cutoff. The algorithm consists of two major steps:
1) Divide the given set of vectors into bins on a regular grid in space.
2) Calculate the distances (or other usefull things) between vectors in 
   neighboring bins.
"""


import math, numpy, copy

__all__ = ["PositionedObject", "SparseBinnedObjects",
           "AnalyseNeighboringObjects", "IntraAnalyseNeighboringObjects",
           "InterAnalyseNeighboringObjects"]


class PositionedObject(object):
    """
    PositionedObject instances are used to feed a SparseBinnedObjects instance
    with information about the vectors that should be binned. Throug the use of
    a PositionedObject, each vector is associated with an id.
    """
    
    def __init__(self, id, vector):
        """
        Initialize a PositionedObject instance.
        
        Arguments:
        id -- A user defined id that is associated with the vector.
        vector -- A numpy array with shape (3,)
        """
        self.id = id
        self.vector = vector


def yield_combinations(l, n):
    if n == 1:
        for item in l:
            yield [item]
    else:
        for item in l:
            for combinations in yield_combinations(l, n-1):
                yield [item] + combinations
    

class SparseBinnedObjects(object):
    """
    A SparseBinnedObjects instance divides 3D space into a sparse grid.
    
    Each cell in the grid is called a bin. Each bin can contain a set of
    Positioned objects. This implementation works with sparse bins: A bin is
    only created in memory when an object is encountered that belongs in that
    bin.
    
    All bins are uniquely defined by their indices i,j,k as defined in __init__.
    """

    deltas = [
         numpy.array(combination, float)
         for combination
         in yield_combinations([-1, 0, 1], 3)
    ]


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
            indices = tuple(numpy.floor(positioned_object.vector*self.reciproke).astype(int))
            bin = self.bins.get(indices)
            if bin is None:
                bin = set()
                self.bins[indices] = bin
            bin.add(positioned_object)

    def yield_surrounding(self, r, deltas=None):
        """
        Iterate over all objects in the bins that surround the bin that
        contains vector r.
        """
        if deltas is None:
            deltas = self.deltas
        center = numpy.floor(r*self.reciproke).astype(int)
        for delta in deltas:
            bin = self.bins.get(tuple(center + delta))
            if bin is not None:
                for positioned_object in bin:
                    yield bin, positioned_object


class AnalyseNeighboringObjects(object):
    """
    AnalyseNeighboringObjects is the base class for 'comparing' vectors between
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
        Intialize a AnalyseNeighboringObjects instance.
        
        compare_function -- for each positioned_object pair that lives in
                            neighboring bins, this function is called.
        """
        self.compare_function = compare_function
        # All these parameters have to be defined by the base class
        self.compare_indices = None
        self.binned_objects1 = None
        self.binned_objects2 = None

    def __call__(self, unit_cell=None):
        if unit_cell is None:
            for result in self.call_delta(numpy.zeros(3, float)):
                yield result
        else:
            active_a, active_b, active_c = unit_cell.cell_active.astype(int)
            for index_a in xrange(-active_a, active_a+1):
                for index_b in xrange(-active_b, active_b+1):
                    for index_c in xrange(-active_c, active_c+1):
                        delta = numpy.dot(unit_cell.cell, numpy.array([index_a, index_b, index_c]))
                        for result in self.call_delta(delta):
                            yield result        
    
    def call_delta(self, delta):
        for center_bin in self.binned_objects1.bins.itervalues():
            for positioned1 in center_bin:
                for neighbor_bin, positioned2 in self.binned_objects2.yield_surrounding(positioned1.vector + delta, self.compare_indices):
                    if self.allow(center_bin, neighbor_bin, positioned1, positioned2):
                        result = self.compare_function(positioned1, positioned2)
                        if result is not None:
                            yield (positioned1, positioned2), result
                        

    def allow(self, bin1, bin2, positioned1, positioned2):
        return True


class IntraAnalyseNeighboringObjects(AnalyseNeighboringObjects):
    """
    IntraAnalyseNeighboringObjects instances compare all vectors within one
    molecule.
    """
    def __init__(self, binned_objects, compare_function):
        AnalyseNeighboringObjects.__init__(self, compare_function)
        self.compare_indices = [(0, 0, 0), (1, 1, 1), 
                                (1, 0, 0), (0, 1, 0), (0, 0, 1),
                                (0, 1, 1), (1, 0, 1), (1, 1, 0), 
                                (0, 1, -1), (-1, 0, 1), (1, -1, 0),
                                (1, 1, -1), (1, -1, -1), (1, -1, 1)]
        self.binned_objects1 = binned_objects
        self.binned_objects2 = binned_objects

    def allow(self, bin1, bin2, positioned1, positioned2):
        return not (bin1 == bin2 and positioned1 >= positioned2)


class InterAnalyseNeighboringObjects(AnalyseNeighboringObjects):
    """
    InterAnalyseNeighboringObjects instances compare 'all' vectors between two
    molecules.
    """
    def __init__(self, binned_objects1, binned_objects2, compare_function):
        AnalyseNeighboringObjects.__init__(self, compare_function)
        assert binned_objects1.gridsize==binned_objects2.gridsize
        self.compare_indices = None
        self.binned_objects1 = binned_objects1
        self.binned_objects2 = binned_objects2
