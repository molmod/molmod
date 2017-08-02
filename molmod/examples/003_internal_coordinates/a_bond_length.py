#!/usr/bin/env python

from __future__ import print_function

from molmod import *

# 0) Load the molecule
mol = Molecule.from_file("dopamine.xyz")
# 1) Compute the length of the hydrogen bond.
# (Atoms 1 and 20 form a hydrogen bond.)
d = bond_length(mol.coordinates[[1, 20]])[0]
# The return value of bond_length is a singleton by default. The final
# part `[0]` takes the first value of this singleton.
print("Hydrogen bond length [Angstrom] =", d/angstrom)
