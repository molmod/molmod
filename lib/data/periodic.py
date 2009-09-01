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


import molmod.units as units
from molmod import context


__all__ = ["AtomInfo", "PeriodicData", "periodic"]


class AtomInfo(object):
    """Data structure for info about an atom"""
    pass


class PeriodicData(object):
    """The entire periodic system

       Objects of the PeriodicData class centralize information about the
       periodic system. The data is loaded during initialization.
    """

    def __init__(self, filename):
        """Initlialize a periodic system structure

           This object is created when importing this module. There is no need
           to do it a second time externally.
        """
        # Initialize empty lists
        self.atoms_by_number = {}
        self.atoms_by_symbol = {}

        self.max_radius = 0.0

        convertors = []
        names = []

        from_unit = {
            "u": (lambda s: float(s)*units.unified),
            "g/cm**3": (lambda s: float(s)*(1e-3*units.gram)/(units.centimeter**3)),
            "a.u.": (lambda s: float(s)),
            "A": (lambda s: float(s)*units.angstrom),
            "pm": (lambda s: float(s)*units.picometer),
        }


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
                    raise TypeError("Can not interpret unit %s." % word)

        f = file(filename)
        lines_read = 0
        for line in f:
            words = line.split()
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
                            setattr(atom_info, name, None)
                        else:
                            value = convertor(word)
                            setattr(atom_info, name, value)
                            if name.endswith("radius") and self.max_radius < value:
                                self.max_radius = value
                    self._add_atom_info(atom_info)
                lines_read += 1
        f.close()

    def _add_atom_info(self, atom_info):
        self.atoms_by_number[atom_info.number] = atom_info
        self.atoms_by_symbol[atom_info.symbol.lower()] = atom_info

    def __len__(self):
        return len(self.atoms_by_symbol)

    def __getitem__(self, index):
        result = self.atoms_by_number.get(index)
        if (result is None) and isinstance(index, basestring):
            return self.atoms_by_symbol.get(index.lower())
        else:
            return result

    def iter_numbers(self):
        """Iterate over all atom numbers in the periodic system"""
        for number in sorted(self.atoms_by_number):
            yield number


periodic = PeriodicData(context.get_share_filename("periodic.csv"))




