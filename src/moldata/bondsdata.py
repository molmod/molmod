# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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

import string, Numeric

__all__ = [
    "BondData", "BOND_SINGLE", "BOND_DOUBLE", "BOND_TRIPLE", "bond_types"
]


BOND_SINGLE = 1
BOND_DOUBLE = 2
BOND_TRIPLE = 3

bond_types = [BOND_SINGLE, BOND_DOUBLE, BOND_TRIPLE]

class BondData(object):
    """
    This object loads information about average bond lengths from a csv file
    when it is initialized. Missing data points in the csv file are estimated
    by adding Van der Waals radii of the two atoms of the given bond.
    """
    bond_tolerance = 1.2

    def __init__(self, filename, periodic_data):
        self.bond_lengths = dict([bond_type, {}] for bond_type in bond_types)
        self.periodic_data = periodic_data
        self.load_bond_data(filename)
        self.approximate_unkown_bond_lengths()
        self.max_length = max(self.bond_lengths.itervalues())
        
    def read_length(self, BOND_TYPE, cells, col):
        """This is a helper method for load_bond_data."""
        nlow = int(cells[2])
        nhigh = int(cells[3])
        for i in range((len(cells) - 4) / 3):
            cell = cells[col + 3 + i*3]
            if cell != 'NA':
                self.bond_lengths[BOND_TYPE][frozenset([nlow, nhigh])] = float(cell)
                return
                
    def load_bond_data(self, filename):
        """
        Load the bond data from the given file.
        
        It's assumed that the uncommented lines in the data file have the
        following format:
        "symbol1 symbol2 number1 number2 bond_length_single_a\
        bond_length_double_a bond_length_triple_a bond_length_single_b\
        bond_length_double_b bond_length_triple_b ..."
        where a, b, ... stand for different sources.
        """
        bond_file = file(filename)
        for line in bond_file:
            cells = string.split(line)
            if (len(cells) > 0) and (cells[0][0] != "#"):
                self.read_length(BOND_SINGLE, cells, 1)
                self.read_length(BOND_DOUBLE, cells, 2)
                self.read_length(BOND_TRIPLE, cells, 3)
        bond_file.close()

    def approximate_unkown_bond_lengths(self):
        dataset = self.bond_lengths[BOND_SINGLE]
        for index1, n1 in enumerate(self.periodic_data.numbers):
            for n2 in self.periodic_data.numbers[:index1]:
                pair = frozenset([n1, n2])
                if (pair not in dataset) and (n1 in self.periodic_data.radius)\
                   and (n2 in self.periodic_data.radius):
                    dataset[pair] = self.periodic_data.radius[n1] + self.periodic_data.radius[n2]
                
    def bonded(self, n1, n2, distance):
        """
        Return the estimated bond type.
        
        This method checks wether for the given pair of atom numbers, the
        given distance corresponds to a certain bond_length. The best
        matching bond type will be returned. If the distance is a factor
        self.bond_tolerance larger than a tabulated distance, the algorithm
        will not relate them.
        """
        if distance > self.max_length * self.bond_tolerance: return None
        deviation = 0.0
        pair = frozenset(n1,n2)
        bond_type = None
        for bt in bond_types:
            bond_length =  self.bond_lengths[bt][pair]
            if distance < bond_length * self.bond_tolerance:
                if bond_type == None:
                    bond_type = bt
                    deviation = abs(bond_length - distance)
                else:
                    new_deviation = abs(bond_length - distance)
                    if deviation > new_deviation:
                        bond_type = bt
                        deviation = new_deviation
        return bond_type
        
    def get_length(self, n1, n2, bond_type):
        """
        Return the length of a bond between n1 and n2 of type bond_type.
        
        This is a safe method for querying a bond_length. If no answer can be
        found, this get_length returns None.
        """
        dataset = self.bond_lengths.get(bond_type)
        if dataset == None:
            return None
        return dataset.get(frozenset([n1, n2]))
    
