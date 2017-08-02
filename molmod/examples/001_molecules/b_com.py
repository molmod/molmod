#!/usr/bin/env python

from molmod import *

# 0) Load the molecule from file
mol1 = Molecule.from_file("ibuprofen.sdf")

# 1) The atomic masses are not included in the sdf file. We need them below to
# compute the center of mass, so we assign standard masses.
mol1.set_default_masses()
print(mol1.masses/amu)

# 2) Create modified coordinates. The attribute coordinates is an N by 3 array,
# while the com attribute is an array with three components. The following
# subtraction is performed elementwise at the level of the rows. Each row
# of the coordinates array corresponds to one atom. From each atomic coordinate
# the center of mass is subtracted to construct a new coordinates array.
new_coordinates = mol1.coordinates - mol1.com

# 3) Make a copy of mol1 with updated coordinates and write it to a file.
mol2 = mol1.copy_with(coordinates=new_coordinates)
mol2.write_to_file("ibuprofen_com.xyz")
print("Written file ibuprofen_com.xyz")
