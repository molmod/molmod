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

from pychem.units import from_unified

import string, numpy, copy

__all__ = ["PeriodicData"]


class PeriodicData(object):
    """
    Objects of the PeriodicData class centralize information about the 
    periodic system. The data is loaded during initialization.
    """
    def __init__(self, filename):
        # Initialize empty lists
        self.name = {}
        self.symbol = {}
        self.symbol_reverse = {}
        self.numbers = []
        self.row = {}
        self.col = {}    
        self.mass = {}
        self.radius = {}
        self.valence = {}
        self.lonepairs = {}
        self.artificial = {}
        self.color = {}

        self.max_radius = 0.0
        
        periodic_file = file(filename)
        for line in periodic_file:
            cells = string.split(line)
            if (len(cells) > 0) and (cells[0][0] != "#"):
                n = int(cells[2])
                self.name[n] = cells[0]
                self.symbol[n] = cells[1]
                self.symbol_reverse[cells[1].lower()] = n
                self.numbers.append(n)
                self.row[n] = int(cells[3])
                self.col[n] = int(cells[4])
                if cells[5] != "NA": self.mass[n] = from_unified(float(cells[5]))
                if cells[7] != "NA": 
                    radius = float(cells[7])
                    if radius > self.max_radius: self.max_radius = radius
                    self.radius[n] = radius
                if cells[9] != "NA": self.valence[n] = int(cells[9])
                if cells[10] != "NA": self.lonepairs[n] = int(cells[10])
                self.artificial[n] = int(cells[11])
                if cells[12] != "NA": self.color[n] = numpy.array([float(cells[12]), float(cells[13]), float(cells[14]), 1.0])
        periodic_file.close()
        
    def symbol_lookup(self, symbol):
        """Return the atom number of the given symbol (case insensitive)."""
        return self.symbol_reverse.get(symbol.lower())   
            
