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
"""Tools for reading SDF (Structure Data Format) files

   For more information about the file format, see
   http://www.epa.gov/ncct/dsstox/MoreonSDF.html
"""


from builtins import range
from builtins import object
import numpy as np

from molmod.units import angstrom
from molmod.periodic import periodic
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.io.common import FileFormatError


__all__ = ["SDFReader"]


class SDFReader(object):
    """A basic reader for SDF files.

       Use this reader as an iterator:

         >>> sr = SDFReader("somefile.sdf")
         >>> for mol in sr:
         ...     print mol.title
    """
    def __init__(self, f):
        """
           Argument:
            | ``f``  --  a filename or a file-like object
        """
        if isinstance(f, str):
            self.filename = f
            self.f = open(f)
            self._auto_close = True
        else:
            # try to treat f as a file-like object and hope for the best.
            self.f = f
            self._auto_close = False

    def __del__(self):
        if self._auto_close:
            self.f.close()

    def __iter__(self):
        return self

    def __next__(self):
        """Load the next molecule from the SDF file

           This method is part of the iterator protocol.
        """
        while True:
            title = next(self.f)
            if len(title) == 0:
                raise StopIteration
            else:
                title = title.strip()
            next(self.f) # skip line
            next(self.f) # skip empty line
            words = next(self.f).split()
            if len(words) < 2:
                raise FileFormatError("Expecting at least two numbers at fourth line.")
            try:
                num_atoms = int(words[0])
                num_bonds = int(words[1])
            except ValueError:
                raise FileFormatError("Expecting at least two numbers at fourth line.")

            numbers = np.zeros(num_atoms, int)
            coordinates = np.zeros((num_atoms, 3), float)
            for i in range(num_atoms):
                words = next(self.f).split()
                if len(words) < 4:
                    raise FileFormatError("Expecting at least four words on an atom line.")
                try:
                    coordinates[i, 0] = float(words[0])
                    coordinates[i, 1] = float(words[1])
                    coordinates[i, 2] = float(words[2])
                except ValueError:
                    raise FileFormatError("Coordinates must be floating point numbers.")
                atom = periodic[words[3]]
                if atom is None:
                    raise FileFormatError("Unrecognized atom symbol: %s" % words[3])
                numbers[i] = atom.number
            coordinates *= angstrom

            edges = []
            orders = np.zeros(num_bonds, int)
            for i in range(num_bonds):
                words = next(self.f).split()
                if len(words) < 3:
                    raise FileFormatError("Expecting at least three numbers on a bond line.")
                try:
                    edges.append((int(words[0])-1, int(words[1])-1))
                    orders[i] = int(words[2])
                except ValueError:
                    raise FileFormatError("Expecting at least three numbers on a bond line.")

            formal_charges = np.zeros(len(numbers), int)

            line = next(self.f)
            while line != "M  END\n":
                if line.startswith("M  CHG"):
                    words = line[6:].split()[1:] # drop the first number which is the number of charges
                    i = 0
                    while i < len(words)-1:
                        try:
                            formal_charges[int(words[i])-1] = int(words[i+1])
                        except ValueError:
                            raise FileFormatError("Expecting only integer formal charges.")
                        i += 2
                line = next(self.f)

            # Read on to the next molecule
            for line in self.f:
                if line == "$$$$\n":
                    break

            molecule = Molecule(numbers, coordinates, title)
            molecule.formal_charges = formal_charges
            molecule.formal_charges.setflags(write=False)
            molecule.graph = MolecularGraph(edges, numbers, orders)
            return molecule
