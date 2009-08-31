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
"""
Binning is a usefull technique for efficiently calculating all distances between
a set of coordinates, when you are only interested in the distances below a given
cutoff. The algorithm consists of two major steps:
1) Divide the given set of coordinates into bins on a regular grid in space.
2) Calculate the distances (or other usefull things) between coordinates in
   neighboring bins.
"""


import numpy, copy

__all__ = ["PositionedObject", "SparseBinnedObjects",
           "AnalyseNeighboringObjects", "IntraAnalyseNeighboringObjects",
           "InterAnalyseNeighboringObjects"]


class PositionedObject(object):
    """
    PositionedObject instances are used to feed a SparseBinnedObjects instance
    with information about the coordinates that should be binned. Throug the use of
    a PositionedObject, each coordinate is associated with an id.
    """

    def __init__(self, id, coordinate):
        """
        Initialize a PositionedObject instance.

        Arguments:
        id -- A user defined id that is associated with the coordinate.
        coordinate -- A numpy array with shape (3,)
        """
        assert isinstance(coordinate, numpy.ndarray)
        assert coordinate.shape == (3,)
        self.id = id
        self.coordinate = coordinate


class SparseBinnedObjects(object):
    """
    A SparseBinnedObjects instance divides 3D space into a sparse grid.

    Each cell in the grid is called a bin. Each bin can contain a set of
    Positioned objects. This implementation works with sparse bins: A bin is
    only created in memory when an object is encountered that belongs in that
    bin.

    All bins are uniquely defined by their indices i,j,k as defined in __init__.
    """

    def __init__(self, positioned_objects, gridsize=1):
        """
        Initialize a BinnedObjects instance.

        Arguments:
        positioned_objects -- An iterator over the PositionedObject instances
                              that have to be binned.
        gridsize -- Defines the size of the bins.
        """
        self.gridsize = gridsize
        self.reciproke = 1/gridsize

        self.bins = {}

        for positioned_object in positioned_objects:
            indices = tuple(numpy.floor(positioned_object.coordinate*self.reciproke).astype(int))
            bin = self.bins.get(indices)
            if bin is None:
                bin = set()
                self.bins[indices] = bin
            bin.add(positioned_object)

    def yield_surrounding(self, r, deltas):
        """
        Iterate over all objects in the bins that surround the bin that
        contains coordinate r.
        """
        center = numpy.floor(r*self.reciproke).astype(int)
        for delta in deltas:
            bin = self.bins.get(tuple(center + delta))
            if bin is not None:
                for positioned_object in bin:
                    yield bin, positioned_object


class AnalyseNeighboringObjects(object):
    """
    AnalyseNeighboringObjects is the base class for 'comparing' coordinates between
    neigbouring bins.
    """

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
            active_a, active_b, active_c = unit_cell.active.astype(int)
            for index_a in xrange(-active_a, active_a+1):
                for index_b in xrange(-active_b, active_b+1):
                    for index_c in xrange(-active_c, active_c+1):
                        delta = unit_cell.to_cartesian(numpy.array([index_a, index_b, index_c]))
                        for result in self.call_delta(delta):
                            yield result

    def call_delta(self, delta):
        for center_bin in self.binned_objects1.bins.itervalues():
            for positioned1 in center_bin:
                for neighbor_bin, positioned2 in self.binned_objects2.yield_surrounding(positioned1.coordinate + delta, self.compare_indices):
                    if self.allow(center_bin, neighbor_bin, positioned1, positioned2):
                        result = self.compare_function(positioned1, positioned2)
                        if result is not None:
                            yield (positioned1, positioned2), result


    def allow(self, bin1, bin2, positioned1, positioned2):
        return True


class IntraAnalyseNeighboringObjects(AnalyseNeighboringObjects):
    """
    IntraAnalyseNeighboringObjects instances compare all coordinates within one
    molecule.
    """
    def __init__(self, binned_objects, compare_function):
        AnalyseNeighboringObjects.__init__(self, compare_function)
        self.compare_indices = numpy.array([
            (0, 0, 0), (1, 1, 1),
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (0, 1, 1), (1, 0, 1), (1, 1, 0),
            (0, 1, -1), (-1, 0, 1), (1, -1, 0),
            (1, 1, -1), (1, -1, -1), (1, -1, 1)
        ], int)
        self.binned_objects1 = binned_objects
        self.binned_objects2 = binned_objects

    def allow(self, bin1, bin2, positioned1, positioned2):
        return not (bin1 == bin2 and id(positioned1) >= id(positioned2))


class InterAnalyseNeighboringObjects(AnalyseNeighboringObjects):
    """
    InterAnalyseNeighboringObjects instances compare 'all' coordinates between two
    molecules.
    """
    def __init__(self, binned_objects1, binned_objects2, compare_function):
        AnalyseNeighboringObjects.__init__(self, compare_function)
        assert binned_objects1.gridsize==binned_objects2.gridsize
        self.compare_indices = numpy.array([
            (-1, -1, -1), (-1, -1,  0), (-1, -1,  1),
            (-1,  0, -1), (-1,  0,  0), (-1,  0,  1),
            (-1,  1, -1), (-1,  1,  0), (-1,  1,  1),

            ( 0, -1, -1), ( 0, -1,  0), ( 0, -1,  1),
            ( 0,  0, -1), ( 0,  0,  0), ( 0,  0,  1),
            ( 0,  1, -1), ( 0,  1,  0), ( 0,  1,  1),

            ( 1, -1, -1), ( 1, -1,  0), ( 1, -1,  1),
            ( 1,  0, -1), ( 1,  0,  0), ( 1,  0,  1),
            ( 1,  1, -1), ( 1,  1,  0), ( 1,  1,  1),
        ], int)
        self.binned_objects1 = binned_objects1
        self.binned_objects2 = binned_objects2




