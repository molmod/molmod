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

from molmod.units import from_unit

import string, numpy, copy

__all__ = ["PeriodicData"]


class AtomInfo(object):
    def __init__(self):
        self.radius = None

    def add_attribute(self, name, value):
        if name.endswith("radius") and self.radius is None:
            self.radius = value
        self.__dict__[name] = value


class PeriodicData(object):
    """
    Objects of the PeriodicData class centralize information about the
    periodic system. The data is loaded during initialization.
    """

    def __init__(self, filename):
        # Initialize empty lists
        self.atoms_by_number = {}
        self.atoms_by_symbol = {}

        self.max_radius = 0.0

        convertors = []
        names = []

        def append_convertor(word):
            if word == "str":
                convertors.append(str)
            elif word == "int":
                convertors.append(int)
            elif word == "float":
                convertors.append(float)
            elif word == "bool":
                convertors.append(eval)
            else:
                convertor = from_unit.get(word)
                if convertor is not None:
                    convertors.append(lambda s: convertor(float(s)))
                else:
                    convertors.append(float)

        f = file(filename)
        lines_read = 0
        for line in f:
            words = string.split(line)
            if (len(words) > 0) and (words[0][0] != "#"):
                if lines_read == 0:
                    # load all the attribute names
                    names = words
                elif lines_read == 1:
                    # load all the conversion factors
                    for word in words:
                        append_convertor(word)
                else:
                    atom_info = AtomInfo()
                    for name, convertor, word in zip(names, convertors, words):
                        if word == "NA":
                            atom_info.add_attribute(name, None)
                        else:
                            atom_info.add_attribute(name, convertor(word))
                    self.add_atom_info(atom_info)
                    if self.max_radius < atom_info.radius:
                        self.max_radius = atom_info.radius
                lines_read += 1
        f.close()
        
    def add_atom_info(self, atom_info):
        self.atoms_by_number[atom_info.number] = atom_info
        self.atoms_by_symbol[atom_info.symbol.lower()] = atom_info

    def __getitem__(self, index):
        result = self.atoms_by_number.get(index)
        if (result is None) and isinstance(index, str):
            return self.atoms_by_symbol.get(index.lower())
        else:
            return result

    def yield_numbers(self):
        for number in self.atoms_by_number:
            yield number
