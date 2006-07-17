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
""" 
periodic is an instance of the PeriodicData class in the module 
molmod.data.periodicdata. It represents the data in periodic.csv.

bonds is an instance of the BondsData class in the module 
molmod.data.bondsdata. It represents the data in bonds.csv.
"""

from periodicdata import PeriodicData
from bondsdata import BondData, BOND_SINGLE, BOND_DOUBLE, BOND_TRIPLE, bond_types

from molmod import context

__all__ = [
    "periodic", 
    "bonds", "BOND_SINGLE", "BOND_DOUBLE", "BOND_TRIPLE", "bond_types"
]

periodic = PeriodicData(context.share_path + "periodic.csv")
bonds = BondData(context.share_path + "bonds.csv", periodic)
