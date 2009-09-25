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
"""Efficiently compute all distances below a cutoff

   Binning is a useful technique for efficiently calculating all distances
   between a number of coordinates when you are only interested in the distances
   below a given cutoff. The algorithm consists of two major steps:

     1) Divide the given set of coordinates into bins on a regular grid
     2) Calculate the distances (or other useful things) between coordinates in
        neighboring bins.
"""


from molmod.unit_cells import UnitCell
import numpy


__all__ = ["PairSearch"]


class PairSearch(object):
    """Efficient iterator over all pairs of coordinates with a distance below a cutoff.

       Example usage:
       >>> coordinates = numpy.random.uniform(0,10,(10,3))
       >>> for i, j, distance, delta in  PairSearch(coordinates, 2.5):
       ...     print i, j, distance
    """

    def __init__(self, coordinates, cutoff, unit_cell=None, grid=None):
        """Initialize a PairSearch object

           Argument:
             coordinates  --  A Nx3 numpy array with Cartesian coordinates
             radius  --  The cutoff radius for the pair distances. Distances
                         larger than the cutoff will be neglected in the
                         pair search.

           Optional argument:
             unit_cell  --  Specifies the periodic boundary conditions
             grid  --  Specification of the grid, can be a floating point number
                       which will result in cubic bins with edge length equal to
                       the given number. Otherwise a UnitCell object can be
                       specified to construct non-cubic bins. In the latter case
                       and when a unit_cell is given, the unit cell vectors must
                       be integer linear combinations of the grid cell vectors
                       (for those directions that are active in the unit cell).
                       If this is not the case, a ValueError is raised.

           The default value of grid depends on other parameters:
             1) When no unit cell is given, it is equal to cutoff/2.9.
             2) When a unit cell is given, the grid cell is as close to cubic
                as possible, with spacings below cutoff/2 that are integer
                divisions of the unit cell spacings


        """
        self.cutoff = cutoff
        self.unit_cell = unit_cell

        if grid is None:
            if unit_cell is None:
                grid = cutoff/2.9
            else:
                # The following would be faster, but the code is not reliable
                # enough yet.
                #grid = unit_cell.get_optimal_subcell(cutoff/2.0)
                divisions = numpy.ceil(unit_cell.spacings/cutoff)
                grid = unit_cell/divisions

        if isinstance(grid, float):
            self.grid_cell = UnitCell(numpy.array([[grid, 0, 0], [0, grid, 0], [0, 0, grid]]))
        elif isinstance(grid, UnitCell):
            self.grid_cell = grid
        else:
            raise TypeError("Grid must be None, a float or a UnitCell instance.")

        if unit_cell is not None:
            # The columns of integer_matrix are the unit cell vectors in
            # fractional coordinates of the grid cell.
            integer_matrix = self.grid_cell.to_fractional(unit_cell.matrix.transpose()).transpose()
            if abs((integer_matrix - numpy.round(integer_matrix))*self.unit_cell.active).max() > 1e-6:
                raise ValueError("The unit cell vectors are not an integer linear combination of grid cell vectors.")
            integer_matrix = integer_matrix.round()
            self.integer_cell = UnitCell(integer_matrix, unit_cell.active)

        self.bins = {}

        fractional = self.grid_cell.to_fractional(coordinates)
        for i in xrange(len(coordinates)):
            if unit_cell is None:
                key = fractional[i].astype(int)
            else:
                key = numpy.round(self.integer_cell.shortest_vector(fractional[i].astype(int))).astype(int)
            key = tuple(key)
            bin = self.bins.get(key)
            if bin is None:
                bin = []
                self.bins[key] = bin
            bin.append((i,coordinates[i]))

        neighbor_indexes = self.grid_cell.get_radius_indexes(cutoff)
        if self.unit_cell is None:
            self.neighbor_indexes = neighbor_indexes
        else:
            self.neighbor_indexes = []
            for index in neighbor_indexes:
                fr_index = self.integer_cell.to_fractional(index)
                if fr_index.max() < 0.5 and fr_index.min() >= -0.5:
                    self.neighbor_indexes.append(index)
            self.neighbor_indexes = numpy.array(self.neighbor_indexes)

    def __iter__(self):
        """Iterate over all pairs with a distance below the cutoff"""
        for key0, bin0 in self.bins.iteritems():
            for key1, bin1 in self.iter_surrounding(key0):
                for i0, c0 in bin0:
                    for i1, c1 in bin1:
                        if i1 >= i0:
                            continue
                        delta = c1 - c0
                        if self.unit_cell is not None:
                            delta = self.unit_cell.shortest_vector(delta)
                        distance = numpy.linalg.norm(delta)
                        if distance <= self.cutoff:
                            yield i0, i1, delta, distance


    def iter_surrounding(self, center_key):
        """Iterate over all bins that surround the given bin"""
        for shift in self.neighbor_indexes:
            key = numpy.add(center_key, shift)
            if self.unit_cell is not None:
                key = numpy.round(self.integer_cell.shortest_vector(key))
            key = tuple(key.astype(int))
            bin = self.bins.get(key)
            if bin is not None:
                yield key, bin


