# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
#
# References for the conversion values:
#    * B. J. Mohr and B. N. Taylor,
#      CODATA recommended values of the fundamental physical
#      constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)
#    * The NIST Reference on Constants, Units, and Uncertainty
#      (http://physics.nist.gov/cuu/Constants/)
#    * 1 calorie = 4.184 Joules
#
# Naming conventions in this module: unit is the value of one external unit
# in internal - i.e. atomic - units. e.g. If you want to have a distance of
# five angstrom in internal units: 5*angstrom. If you want to convert a length
# of 5 internal units to angstrom: 5/angstrom. It is recommended to perform
# this kind of conversions, only when data is read from the input and data is
# written to the output.


from constants import avogadro, lightspeed


def parse_unit(s):
    return eval(s)


# *** Charge ***

C = 1/1.602176462e-19
# for compatibility with previous verions
coulomb = C

# Mol

mol = avogadro

# *** Mass ***

kg = 1/9.10938188e-31

g = 1e-3*kg
mg = 1e-3*g
u = 1e-3*kg/mol
# for compatibility with previous verions
unified = u

# *** Length ***

m = 1/0.5291772083e-10

cm = 1e-2*m
mm = 1e-3*m
um = 1e-6*m
nm = 1e-9*m
A = 1e-10*m
pm = 1e-12*m
# for compatibility with previous verions
meter = m
angstrom = A
nanometer = nm


# *** Energy ***

J = 1/4.35974381e-18

cal = 4.184*J
kJmol = 1e3*J/mol
kcalmol = 1e3*cal/mol
eV = (1.0/C)*J
# for compatibility with previous verions
kjmol = kJmol
joule = J
calorie = cal
ev = eV

# *** Angles ***

degree = 0.017453292519943295

# *** Time ***

s = 1/2.418884326500e-17

ns = 1e-9*s
fs = 1e-15*s
ps = 1e-12*s
# for compatibility with previous verions
second = s
nanosecond = ns
picosecond = ps
femtosecond = fs

# *** Frequency ***

invcm = lightspeed/(0.01*meter)

# *** Temperature ***

K = 1
kelvin = K



