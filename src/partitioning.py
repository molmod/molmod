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


import numpy

from molmod.helpers import partitioning as partitioningf90
from molmod.data import periodic


def set_hirshfeld_weights(grid, atom_grids, label="weights"):
    atom_grid_size = None
    for atom_grid in atom_grids:
        if atom_grid_size is None:
            atom_grid_size = len(atom_grid.points)
        elif atom_grid_size != len(atom_grid.points):
            raise HirshfeldException("All atom densities have to be defined on grids with an equal number of datapoints.")

    atom_indices = numpy.zeros(len(grid.molecule.numbers), float) - 1
    atom_densities = numpy.zeros((len(atom_grids), 2, atom_grid_size), float)
    for index, atom_grid in enumerate(atom_grids):
        for counter, number in enumerate(grid.molecule.numbers):
            if number == atom_grid.molecule.numbers[0]:
                atom_indices[counter] = index
        atom_densities[index,0] = numpy.sqrt(((atom_grid.points - atom_grid.molecule.coordinates[0])**2).sum(axis=1))
        atom_densities[index,1] = atom_grid.density
    if (atom_indices == -1).any():
        missing_numbers = grid.molecule.numbers.compress(atom_indices == -1)
        raise HirshfeldException("The density profiles for the folowing atoms are missing:" % (
            " ".join(periodic[number].symbol for number in missing_numbers)
        ))
    grid.__dict__[label] = partitioningf90.hirshfeld_weights(grid.points, atom_densities, atom_indices, grid.molecule.coordinates)


def set_voronoi_weights(grid, weights_label="weights"):
    grid.__dict__[weights_label] = partitioningf90.voronoi_weights(grid.points, grid.molecule.coordinates)


def charges_and_dipoles(grid, weights_label="weights"):
    result = partitioningf90.charges_and_dipoles(grid.points, grid.volumes, grid.density, grid.__dict__[weights_label], grid.molecule.coordinates)
    result *= -1
    result[:,0] += grid.molecule.numbers
    return result

def exterior_multipoles(grid, weights_label="weights"):
    result = partitioningf90.ext_multipoles(grid.points, grid.volumes, grid.density, grid.__dict__[weights_label], grid.molecule.coordinates)
    result *= -1
    result[:,0] += grid.molecule.numbers
    return result

def interior_multipoles(grid, weights_label="weights"):
    return partitioningf90.int_multipoles(grid.points, grid.volumes, grid.density, grid.__dict__[weights_label], grid.molecule.coordinates, grid.molecule.numbers)

def partitioned_interactions(grid, weights_label="weights"):
    return partitioningf90.partitioned_interactions(grid.points, grid.volumes, grid.density, grid.__dict__[weights_label], grid.molecule.coordinates, grid.molecule.numbers)


def interaction_ext_int(extm, intm):
    return (extm*intm).sum()
    signs = numpy.zeros(len(extm), float)
    lmax = int(numpy.sqrt(len(extm)))-1
    counter = 0
    for l in xrange(lmax+1):
        signs[counter] = 1
        counter += 1
        for m in xrange(l):
            signs[counter] = 1
            counter += 1
            signs[counter] = -1
            counter += 1
    return (extm*intm*signs).sum()
