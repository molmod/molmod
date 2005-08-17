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
    bond_tolerance = 1.2

    def __init__(self, filename, periodic_data):
        self.bond_lengths = dict([bond_type, {}] for bond_type in bond_types)
        self.periodic_data = periodic_data
        self.load_bond_data(filename)
        self.approximate_unkown_bond_lengths()
        self.max_length = max(self.bond_lengths.itervalues())
        
    def read_length(self, BOND_TYPE, cells, col):
        for i in range((len(cells) - 4) / 3):
            cell = cells[col + 3 + i*3]
            if cell != 'NA':
                self.bond_lengths[BOND_TYPE][frozenset([nlow, nhigh])] = float(cell)
                return
                
    def load_bond_data(self, filename):
        bond_file = file(filename)
        for line in bond_file:
            cells = string.split(line)
            if (len(cells) > 0) and (cells[0][0] != "#"):
                nlow = int(cells[2])
                nhigh = int(cells[3])
                self.read_length(BOND_SINGLE, cells, 1)
                self.read_length(BOND_DOUBLE, cells, 2)
                self.read_length(BOND_TRIPLE, cells, 3)
        bond_file.close()

    def approximate_unkown_bond_lengths(self):
        for index1, n1 in enumerate(self.periodic_data.numbers):
            for n2 in self.periodic_data.numbers[index1+1:]
                pair = frozenset([n1, n2])
                if pair not in self.bond_lengths:
                    self.bond_lengths[pair] = self.periodic_data.radius[n1] + self.periodic_data.radius[n2]
                
    def bonded(self, n1, n2, length):
        if length > self.max_length * self.bond_tolerance: return None
        deviation = 0.0
        pair = frozenset(n1,n2)
        bond_type = None
        for bt in bond_types:
            bond_length =  self.bond_lengths[bt][pair]
            if length < bond_length * self.bond_tolerance:
                if bond_type == None:
                    bond_type = bt
                    deviation = abs(bond_length - length)
                else:
                    new_deviation = abs(bond_length - length)
                    if deviation > new_deviation:
                        bond_type = bt
                        deviation = new_deviation
        return bond_type
        
    def get_length(self, n1, n2, bond_type):
        dataset = self.bond_lengths.get(bond_type)
        if dataset == None:
            return None
        return dataset.get(frozenset([n1, n2]))
    
