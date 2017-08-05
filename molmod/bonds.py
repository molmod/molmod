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
"""Bond length database and bond type definitions

   An object ``bonds`` of the class ``BondData`` is created upon importing this
   module. It loads information about average bond lengths from a csv file
   when it is initialized. Missing data points in the csv file are estimated
   by adding Van der Waals radii of the two atoms of the given bond.

   This module also defines a few constants for different bond types:

   =============  =====
   Name           Value
   =============  =====
   BOND_SINGLE    0
   BOND_DOUBLE    1
   BOND_TRIPLE    2
   BOND_HYBRID    3
   BOND_HYDROGEN  4
   =============  =====
"""


from __future__ import division

from builtins import range
import pkg_resources

from molmod.periodic import periodic
import molmod.units as units


__all__ = [
    "BOND_SINGLE", "BOND_DOUBLE", "BOND_TRIPLE", "BOND_HYBRID",
    "BOND_HYDROGEN", "bond_types", "BondData", "bonds"
]


BOND_SINGLE = 1
BOND_DOUBLE = 2
BOND_TRIPLE = 3
BOND_HYBRID = 4
BOND_HYDROGEN = 5

bond_types = [
    BOND_SINGLE, BOND_DOUBLE, BOND_TRIPLE, BOND_HYBRID, BOND_HYDROGEN
]


class BondData(object):
    """Database with bond lengths"""
    bond_tolerance = 1.2

    def __init__(self):
        """
           This object is created when importing this module. There is no need
           to do it a second time manually.
        """
        self.lengths = dict([bond_type, {}] for bond_type in bond_types)
        self._load_bond_data()
        self._approximate_unkown_bond_lengths()
        self.max_length = max(
            max(lengths.values())
            for lengths
            in self.lengths.values()
            if len(lengths) > 0
        )

    def _load_bond_data(self):
        """Load the bond data from the given file

           It's assumed that the uncommented lines in the data file have the
           following format:
           symbol1 symbol2 number1 number2 bond_length_single_a bond_length_double_a bond_length_triple_a bond_length_single_b bond_length_double_b bond_length_triple_b ..."
           where a, b, ... stand for different sources.
        """

        def read_units(unit_names):
            """convert unit_names into conversion factors"""
            tmp = {
                "A": units.angstrom,
                "pm": units.picometer,
                "nm": units.nanometer,
            }
            return [tmp[unit_name] for unit_name in unit_names]

        def read_length(BOND_TYPE, words, col):
            """Read the bondlengths from a single line in the data file"""
            nlow = int(words[2])
            nhigh = int(words[3])
            for i, conversion in zip(range((len(words) - 4) // 3), conversions):
                word = words[col + 3 + i*3]
                if word != 'NA':
                    self.lengths[BOND_TYPE][frozenset([nlow, nhigh])] = float(word)*conversion
                    return

        with pkg_resources.resource_stream(__name__, 'data/bonds.csv') as f:
            for line in f:
                words = line.decode('utf-8').split()
                if (len(words) > 0) and (words[0][0] != "#"):
                    if words[0] == "unit":
                        conversions = read_units(words[1:])
                    else:
                        read_length(BOND_SINGLE, words, 1)
                        read_length(BOND_DOUBLE, words, 2)
                        read_length(BOND_TRIPLE, words, 3)

    def _approximate_unkown_bond_lengths(self):
        """Completes the bond length database with approximations based on VDW radii"""
        dataset = self.lengths[BOND_SINGLE]
        for n1 in periodic.iter_numbers():
            for n2 in periodic.iter_numbers():
                if n1 <= n2:
                    pair = frozenset([n1, n2])
                    atom1 = periodic[n1]
                    atom2 = periodic[n2]
                    #if (pair not in dataset) and hasattr(atom1, "covalent_radius") and hasattr(atom2, "covalent_radius"):
                    if (pair not in dataset) and (atom1.covalent_radius is not None) and (atom2.covalent_radius is not None):
                        dataset[pair] = (atom1.covalent_radius + atom2.covalent_radius)
                    #print "%3i  %3i  %s %30s %30s" % (n1, n2, dataset.get(pair), atom1, atom2)

    def bonded(self, n1, n2, distance):
        """Return the estimated bond type

           Arguments:
            | ``n1``  --  the atom number of the first atom in the bond
            | ``n2``  --  the atom number of the second atom the bond
            | ``distance``  --  the distance between the two atoms

           This method checks whether for the given pair of atom numbers, the
           given distance corresponds to a certain bond length. The best
           matching bond type will be returned. If the distance is a factor
           ``self.bond_tolerance`` larger than a tabulated distance, the
           algorithm will not relate them.
        """
        if distance > self.max_length * self.bond_tolerance:
            return None

        deviation = 0.0
        pair = frozenset([n1, n2])
        result = None
        for bond_type in bond_types:
            bond_length = self.lengths[bond_type].get(pair)
            if (bond_length is not None) and \
               (distance < bond_length * self.bond_tolerance):
                if result is None:
                    result = bond_type
                    deviation = abs(bond_length - distance)
                else:
                    new_deviation = abs(bond_length - distance)
                    if deviation > new_deviation:
                        result = bond_type
                        deviation = new_deviation
        return result

    def get_length(self, n1, n2, bond_type=BOND_SINGLE):
        """Return the length of a bond between n1 and n2 of type bond_type

           Arguments:
            | ``n1``  --  the atom number of the first atom in the bond
            | ``n2``  --  the atom number of the second atom the bond

           Optional argument:
            | ``bond_type``  --  the type of bond [default=BOND_SINGLE]

           This is a safe method for querying a bond_length. If no answer can be
           found, this get_length returns None.
        """
        dataset = self.lengths.get(bond_type)
        if dataset == None:
            return None
        return dataset.get(frozenset([n1, n2]))


bonds = BondData()
