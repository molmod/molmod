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

from molmod.helpers import hirshfeld as hirshfeldf90


def hirshfeld_cd(density_cube, atom_cubes):
    m = density_cube.molecule
    atom_indices = numpy.zeros(len(m.numbers), float)
    atom_rhos = []
    for index, (number, atom_cube) in enumerate(atom_cubes):
        for atom_counter, atom_number in enumerate(m.numbers):
            if atom_number == number:
                atom_indices[atom_counter] = index
        atom_rho = numpy.zeros((2, len(atom_cube.grid_data)), float)
        atom_rho[0] = numpy.sqrt(((atom_cube.grid_data[:,0:3] - atom_cube.molecule.coordinates[0])**2).sum(axis=1))
        atom_rho[1] = atom_cube.grid_data[:,3]
        atom_rhos.append(atom_rho)
    atom_rho = numpy.concatenate(atom_rhos)
    result = hirshfeldf90.hirshfeld_cd(
        density_cube.grid_data,
        density_cube.volumes,
        atom_rho,
        atom_indices,
        m.coordinates,
    )
    result *= -1
    result[:,0] += m.numbers
    return result
