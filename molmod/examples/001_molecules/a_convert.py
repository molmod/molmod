#!/usr/bin/env python

from molmod import *

mol = Molecule.from_file("ibuprofen.sdf")
mol.write_to_file("ibuprofen.xyz")
