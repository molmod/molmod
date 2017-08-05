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
"""Efficiently compute all distances below a cutoff

   Binning is a useful technique for efficiently calculating all distances
   between a number of coordinates when you are only interested in the distances
   below a given cutoff. The algorithm consists of two major steps:

     1) Divide the given set of coordinates into bins on a regular grid
     2) Calculate the distances (or other useful things) between coordinates in
        neighboring bins.
"""


from __future__ import division

from builtins import range
import numpy as np

from molmod.unit_cells import UnitCell


__all__ = ["PairSearchIntra", "PairSearchInter"]


class Binning(object):
    """Division of coordinates in regular bins"""
    def __init__(self, coordinates, cutoff, grid_cell, integer_cell=None):
        """Initialize a Binning object

           Arguments:
            | ``coordinates``  --  a Nx3 numpy array with coordinates to be
                                   triaged into bins
            | ``cutoff``  --  The maximum distance between coordinates pairs.
                              This affects the variable self.neighbor_bins,
                              which is a set of relative integer bin coordinates
                              that lie within the cuttof of the central bin.
            | ``grid_cell``  --  A unit cell  object specifying the size and
                                 shape of the bins

           Optional argument:
            | ``integer_cell``  --  the periodicity of the system in terms if
                                    integer grid cells.
        """
        self.grid_cell = grid_cell
        self.integer_cell = integer_cell

        # setup the bins
        self._bins = {}

        fractional = grid_cell.to_fractional(coordinates)
        for i in range(len(coordinates)):
            key = tuple(fractional[i].astype(int))
            if integer_cell is not None:
                key = self.wrap_key(key)
            bin = self._bins.get(key)
            if bin is None:
                bin = []
                self._bins[key] = bin
            bin.append((i, coordinates[i]))

        # compute the neigbouring bins within the cutoff
        if self.integer_cell is None:
            self.neighbor_indexes = grid_cell.get_radius_indexes(cutoff)
        else:
            max_ranges = np.diag(self.integer_cell.matrix).astype(int)
            max_ranges[True^self.integer_cell.active] = -1
            self.neighbor_indexes = grid_cell.get_radius_indexes(cutoff, max_ranges)

    def __iter__(self):
        """Iterate over (key,bin) pairs"""
        return iter(self._bins.items())

    def iter_surrounding(self, center_key):
        """Iterate over all bins surrounding the given bin"""
        for shift in self.neighbor_indexes:
            key = tuple(np.add(center_key, shift).astype(int))
            if self.integer_cell is not None:
                key = self.wrap_key(key)
            bin = self._bins.get(key)
            if bin is not None:
                yield key, bin

    def wrap_key(self, key):
        """Translate the key into the central cell

           This method is only applicable in case of a periodic system.
        """
        return tuple(np.round(
            self.integer_cell.shortest_vector(key)
        ).astype(int))


class PairSearchBase(object):
    """Base class for :class:`PairSearchIntra` and :class:`PairSearchInter`"""
    def _setup_grid(self, cutoff, unit_cell, grid):
        """Choose a proper grid for the binning process"""
        if grid is None:
            # automatically choose a decent grid
            if unit_cell is None:
                grid = cutoff/2.9
            else:
                # The following would be faster, but it is not reliable
                # enough yet.
                #grid = unit_cell.get_optimal_subcell(cutoff/2.0)
                divisions = np.ceil(unit_cell.spacings/cutoff)
                divisions[divisions<1] = 1
                grid = unit_cell/divisions

        if isinstance(grid, float):
            grid_cell = UnitCell(np.array([
                [grid, 0, 0],
                [0, grid, 0],
                [0, 0, grid]
            ]))
        elif isinstance(grid, UnitCell):
            grid_cell = grid
        else:
            raise TypeError("Grid must be None, a float or a UnitCell instance.")

        if unit_cell is not None:
            # The columns of integer_matrix are the unit cell vectors in
            # fractional coordinates of the grid cell.
            integer_matrix = grid_cell.to_fractional(unit_cell.matrix.transpose()).transpose()
            if abs((integer_matrix - np.round(integer_matrix))*self.unit_cell.active).max() > 1e-6:
                raise ValueError("The unit cell vectors are not an integer linear combination of grid cell vectors.")
            integer_matrix = integer_matrix.round()
            integer_cell = UnitCell(integer_matrix, unit_cell.active)
        else:
            integer_cell = None

        return grid_cell, integer_cell


class PairSearchIntra(PairSearchBase):
    """Iterator over all pairs of coordinates with a distance below a cutoff.

       Example usage::

           coordinates = np.random.uniform(0,10,(10,3))
           for i, j, delta, distance in PairSearchIntra(coordinates, 2.5):
               print i, j, distance

       Note that for periodic systems the minimum image convention is applied.
    """

    def __init__(self, coordinates, cutoff, unit_cell=None, grid=None):
        """
           Arguments:
            | ``coordinates``  --  A Nx3 numpy array with Cartesian coordinates
            | ``radius``  --  The cutoff radius for the pair distances.
                              Distances larger than the cutoff will be neglected
                              in the pair search.

           Optional arguments:
            | ``unit_cell``  --  Specifies the periodic boundary conditions
            | ``grid``  --  Specification of the grid, can be a floating point
                        number which will result in cubic bins with edge length
                        equal to the given number. Otherwise a UnitCell object
                        can be specified to construct non-cubic bins. In the
                        latter case and when a unit_cell is given, the unit cell
                        vectors must be integer linear combinations of the grid
                        cell vectors (for those directions that are active in
                        the unit cell). If this is not the case, a ValueError is
                        raised.

           The default value of grid depends on other parameters:

             1) When no unit cell is given, it is equal to cutoff/2.9.
             2) When a unit cell is given, the grid cell is as close to cubic
                as possible, with spacings below cutoff/2 that are integer
                divisions of the unit cell spacings
        """
        self.cutoff = cutoff
        self.unit_cell = unit_cell
        grid_cell, integer_cell = self._setup_grid(cutoff, unit_cell, grid)
        self.bins = Binning(coordinates, cutoff, grid_cell, integer_cell)

    def __iter__(self):
        """Iterate over all pairs with a distance below the cutoff"""
        for key0, bin0 in self.bins:
            for key1, bin1 in self.bins.iter_surrounding(key0):
                for i0, c0 in bin0:
                    for i1, c1 in bin1:
                        if i1 >= i0:
                            continue
                        delta = c1 - c0
                        if self.unit_cell is not None:
                            delta = self.unit_cell.shortest_vector(delta)
                        distance = np.linalg.norm(delta)
                        if distance <= self.cutoff:
                            yield i0, i1, delta, distance

class PairSearchInter(PairSearchBase):
    """Iterator over all pairs of coordinates with a distance below a cutoff.

       Example usage::

           coordinates0 = np.random.uniform(0,10,(10,3))
           coordinates1 = np.random.uniform(0,10,(10,3))
           for i, j, delta, distance in PairSearchInter(coordinates0, coordinates1, 2.5):
               print i, j, distance

       Note that for periodic systems the minimum image convention is applied.
    """

    def __init__(self, coordinates0, coordinates1, cutoff, unit_cell=None, grid=None):
        """
           Arguments:
            | ``coordinates0``  --  A Nx3 numpy array with Cartesian coordinates
            | ``coordinates1``  --  A Nx3 numpy array with Cartesian coordinates
            | ``radius``  --  The cutoff radius for the pair distances.
                              Distances larger than the cutoff will be neglected
                              in the pair search.

           Optional arguments:
            | ``unit_cell``  --  Specifies the periodic boundary conditions
            | ``grid``  --  Specification of the grid, can be a floating point
                        number which will result in cubic bins with edge length
                        equal to the given number. Otherwise a UnitCell object
                        can be specified to construct non-cubic bins. In the
                        latter case and when a unit_cell is given, the unit cell
                        vectors must be integer linear combinations of the grid
                        cell vectors (for those directions that are active in
                        the unit cell). If this is not the case, a ValueError is
                        raised.

           The default value of grid depends on other parameters:
             1) When no unit cell is given, it is equal to cutoff/2.9.
             2) When a unit cell is given, the grid cell is as close to cubic
                as possible, with spacings below cutoff/2 that are integer
                divisions of the unit cell spacings
        """
        self.cutoff = cutoff
        self.unit_cell = unit_cell
        grid_cell, integer_cell = self._setup_grid(cutoff, unit_cell, grid)
        self.bins0 = Binning(coordinates0, cutoff, grid_cell, integer_cell)
        self.bins1 = Binning(coordinates1, cutoff, grid_cell, integer_cell)

    def __iter__(self):
        """Iterate over all pairs with a distance below the cutoff"""
        for key0, bin0 in self.bins0:
            for key1, bin1 in self.bins1.iter_surrounding(key0):
                for i0, c0 in bin0:
                    for i1, c1 in bin1:
                        delta = c1 - c0
                        if self.unit_cell is not None:
                            delta = self.unit_cell.shortest_vector(delta)
                        distance = np.linalg.norm(delta)
                        if distance <= self.cutoff:
                            yield i0, i1, delta, distance
