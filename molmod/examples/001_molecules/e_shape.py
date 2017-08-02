#!/usr/bin/env python

from __future__ import print_function

from molmod import *

from numpy import dot, sqrt
from numpy.linalg import eigvalsh


# 0) Load the molecule.
mol = Molecule.from_file("ibuprofen.sdf")

# 1) Compute the arithmetic center of the atomic coordinates.
center = mol.coordinates.mean(axis=0)
# Without the axis=0 argument, the average over all X, Y and Z components would
# be computed. Now it just computes the averages for each column. axis=1 would
# refer to averaging over rows.

# 2) Move the arithmetic center of the coordinates to the origin. The same
# comments from the b_com.py apply.
centered = mol.coordinates - mol.coordinates.mean(axis=0)

# 3) Compute the covariance matrix of the centered coordinates.
covar = dot(centered.transpose(), centered)/mol.size

# 4) Compute the eigenvalues of the symmetric covariance matrix.
evals = eigvalsh(covar)

# 5) The spread along the three eigenvectors is computed as the standard deviation
# of the atomic coordinates along that direction:
c, b, a = sqrt(evals)
print("Spread along the long axis [A]:", a/angstrom)
print("Spread along the intermediate axis [A]:", b/angstrom)
print("Spread along the short axis [A]:", c/angstrom)

# 6) Test in which category this shape belongs. The factor R is set to 1.5.
R = 1.5
if b < R*c:
    if a < R*b:
        shape = 'equant' # sphere-like
    else:
        shape = 'prolate' # sigar-like
else:
    if a < R*b:
        shape = 'oblate' # pancake-like
    else:
        shape = 'bladed' # keyboard-like
print("The shape category:", shape)
