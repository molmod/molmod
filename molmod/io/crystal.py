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
"""Tools for reading output from the Crystal 06 example API program."""


from __future__ import absolute_import

from builtins import range
import numpy as np

from .common import FileFormatError
from molmod import UnitCell, Molecule


__all__ = ['CrystalAPIOut']


class CrystalAPIOut(object):
    """Reader for output files generate by the Crystal 06 example API program.

       Only those fields related to the the density matrix are read, i.e. also
       the basis set definition etc.
    """
    def __init__(self, filename):
        """
           Arguments:
            | ``filename`` -- The file to load.
        """
        self.filename = filename
        # auxiliary skip function
        def skip_to(f, linestart):
            while True:
                line = f.readline()
                if line.startswith(linestart):
                    return line
                if len(line) == 0:
                    return

        with open(filename) as f:
            # Unit cell parameters
            line = skip_to(f, ' DIRECT LATTICE VECTOR COMPONENTS (BOHR)')
            if line is None:
                raise FileFormatError('Could not find the lattice vectors')
            f.readline()
            f.readline()
            vectors = []
            for i in range(3):
                line = f.readline()
                vectors.append([float(word) for word in line.split()[1:]])
            vectors = np.array(vectors)
            self.unit_cell = UnitCell(vectors.transpose())
            # Atomic numbers and coordinates
            line = skip_to(f, ' ATOM CARTESIAN COORDINATES (BOHR)')
            if line is None:
                raise FileFormatError('Could not find the atomic coordinates')
            f.readline()
            f.readline()
            f.readline()
            numbers = []
            symbols = []
            coordinates = []
            while True:
                line = f.readline()
                if line.startswith(' *****'):
                    break
                words = line.split()
                numbers.append(int(words[1]))
                symbols.append(words[2])
                coordinates.append([float(words[3]), float(words[4]), float(words[5])])
            self.mol = Molecule(np.array(numbers), np.array(coordinates), symbols=symbols, unit_cell=self.unit_cell)
            # Basis set specification
            line = skip_to(f, ' VARIATIONAL BASIS SET')
            if line is None:
                raise FileFormatError('Could not find the basis set')
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            self.basisset = {}
            last_basis = None
            last_contraction = None
            while True:
                line = f.readline()
                if line.startswith(' *****'):
                    break
                if line.startswith('         '):
                    # add info to the last atomic basis set
                    assert last_basis is not None
                    subshell = line[36:40].strip()
                    if len(subshell) == 0:
                        subshell = last_contraction[0]
                        # add a primitive to the contraction
                        exponent = float(line[40:50])
                        if subshell == 'S':
                            values = exponent, float(line[50:60])
                        elif subshell == 'SP':
                            values = exponent, float(line[50:60]), float(line[60:70])
                        else:
                            values = exponent, float(line[70:80])
                        last_contraction[1].append(values)
                    else:
                        # define a new contraction
                        last_contraction = (subshell, [])
                        last_basis.append(last_contraction)
                else:
                    # add new atoms
                    symbol = line.split()[1]
                    if symbol not in self.basisset:
                        last_basis = []
                        self.basisset[symbol] = last_basis
            # Compute the total number of basis functions (and orbitals).
            self.num_basis = 0
            subshell_counts = {'S': 1, 'P': 3, 'SP': 4, 'D': 5, 'F': 7, 'G': 9}
            for symbol, basis in self.basisset.items():
                symbol_count = symbols.count(symbol)
                for subshell, contraction in basis:
                    self.num_basis += symbol_count*subshell_counts[subshell]
            # Density matrix.
            line = skip_to(f, ' DENSITY MATRIX DIRECT LATTICE')
            if line is None:
                raise FileFormatError('Could not find the density matrix')
            f.readline()
            f.readline()
            f.readline()
            self.density_matrix = np.zeros((self.num_basis, self.num_basis), float)
            j0 = 0
            while j0 < self.num_basis:
                f.readline()
                f.readline()
                f.readline()
                for i in range(self.num_basis):
                    line = f.readline()
                    words = line.split()[1:]
                    for j1, word in enumerate(words):
                        self.density_matrix[i,j0 + j1] = float(word)
                j0 += 10
            # Ignore the remainder.
