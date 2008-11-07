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


from molmod.units import angstrom
from molmod.data.periodic import periodic
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph

import numpy


__all__ = ["SDFError", "SDFReader"]


class SDFError(Exception):
    pass


class SDFReader(object):
    def __init__(self, f):
        if isinstance(f, basestring):
            self.filename = f
            self.f = file(f)
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

    def next(self):
        while True:
            title = self.f.next()
            if len(title) == 0:
                raise StopIteration
            else:
                title = title.strip()
            self.f.next() # skip line
            self.f.next() # skip empty line
            words = self.f.next().split()
            if len(words) < 2:
                raise SDFError("Expecting at least two numbers at fourth line.")
            try:
                num_atoms = int(words[0])
                num_bonds = int(words[1])
            except ValueError:
                raise SDFError("Expecting at least two numbers at fourth line.")

            numbers = numpy.zeros(num_atoms, int)
            coordinates = numpy.zeros((num_atoms,3), float)
            for i in xrange(num_atoms):
                words = self.f.next().split()
                if len(words) < 4:
                    raise SDFError("Expecting at least four words on an atom line.")
                try:
                    coordinates[i,0] = float(words[0])
                    coordinates[i,1] = float(words[1])
                    coordinates[i,2] = float(words[2])
                except ValueError:
                    raise SDFError("Coordinates must be floating point numbers.")
                atom = periodic[words[3]]
                if atom is None:
                    raise SDFError("Unrecognized atom symbol: %s" % words[3])
                numbers[i] = atom.number
            coordinates *= angstrom

            pairs = []
            orders = numpy.zeros(num_bonds, int)
            for i in xrange(num_bonds):
                words = self.f.next().split()
                if len(words) < 3:
                    raise SDFError("Expecting at least three numbers on a bond line.")
                try:
                    pairs.append((int(words[0])-1,int(words[1])-1))
                    orders[i] = int(words[2])
                except ValueError:
                    raise SDFError("Expecting at least three numbers on a bond line.")

            formal_charges = numpy.zeros(len(numbers), int)

            line = self.f.next()
            while line != "M  END\n":
                if line.startswith("M  CHG"):
                    words = line[6:].split()[1:] # drop the first number which is the number of charges
                    i = 0
                    while i < len(words)-1:
                        try:
                            formal_charges[int(words[i])-1] = int(words[i+1])
                        except ValueError:
                            raise SDFError("Expecting only integer formal charges.")
                        i += 2
                line = self.f.next()

            # Read on to the next molecule
            for line in self.f:
                if line == "$$$$\n":
                    break

            molecule = Molecule(numbers, coordinates, title)
            molecule.formal_charges = formal_charges
            molecule.charge = formal_charges.sum()
            molecule.graph = MolecularGraph(pairs, numbers, orders)
            return molecule



