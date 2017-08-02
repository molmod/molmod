#!/usr/bin/env python

from __future__ import print_function

from molmod import *
from numpy import sqrt, pi

m = 1.0*amu
k = 3200*kjmol/angstrom**2

freq = sqrt(k/m)/(2*pi)

print("Force constant [kcal mol^-1 angstrom^-2]:", k/(kcalmol/angstrom**2))
print("Wavenumber [cm^-1]:", (freq/lightspeed*centimeter))
