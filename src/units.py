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
#
# Source of numbers: http://physics.nist.gov/cuu/Constants/
# Except for calorie: 1 calorie = 4.184 Joules (thank google) :-)

ATOMARY = 0

LENGTH = 0
ANGSTROM = 1
NANOMETER = 2

ENERGY = 1
JOULE = 3
CALORIE = 4
KJMOL = 5
KCALMOL = 6
EV = 7

MASS = 2
UNIFIED = 8

CHARGE = 3
COULOMB = 9

ANGLE = 4
RADIALS = 10
DEGREES = 11

suffices = {
    ATOMARY: "a.u.",
    ANGSTROM: "A",
    NANOMETER: "nm",
    JOULE: "J",
    CALORIE: "cal",
    KJMOL: "kJ/mol",
    KCALMOL: "kcal/mol",
    EV: "eV",
    UNIFIED: "u",
    COULOMB: "C",
    RADIALS: "rad",
    DEGREES: "°"
}

measures = {
    LENGTH: [ATOMARY, ANGSTROM, NANOMETER],
    ENERGY: [ATOMARY, JOULE, CALORIE, KJMOL, KCALMOL, EV],
    MASS: [ATOMARY, UNIFIED],
    CHARGE: [ATOMARY, COULOMB],
    ANGLE: [RADIALS, DEGREES]
}

measure_names = {
    LENGTH: "Length",
    ENERGY: "Energy",
    MASS: "Mass",
    CHARGE: "Charge",
    ANGLE: "Angle"
}

## Conversion to internal units (a.u)

# Length

angstrom = lambda x: x * 1.8897261249
nanometer = lambda x: x * 18.897261249

# Energy

joule = lambda x: x * 2.2937125689189235e+17
calorie = lambda x: x * 9.596893388356777e+17
kjmol = lambda x: x * 0.00038087988615327677
kcalmol = lambda x: x * 0.0015859838459422444
ev = lambda x: x * 0.036749324444879071

# Mass

unified = lambda x: x * 1822.8884798405547

# Charge

coulomb = lambda x: x * 6.2415094796077179e+18

# Angles

degrees = lambda x: x * 0.017453292519943295

# In a dictionary

unit = {
    ATOMARY: lambda x: x,
    ANGSTROM: angstrom,
    NANOMETER: nanometer,
    JOULE: joule,
    CALORIE: calorie,
    KJMOL: kjmol,
    KCALMOL: kcalmol,
    EV: ev,
    UNIFIED: unified,
    COULOMB: coulomb,
    RADIALS: lambda x: x,
    DEGREES: degrees
}

## Conversion from internal units (a.u.)

# Length

to_angstrom = lambda x: x * 0.5291772108
to_nanometer = lambda x: x * 0.05291772108

# Energy

to_joule = lambda x: x * 4.35974417e-18
to_calorie = lambda x: x * 1.0420038647227532e-18
to_kjmol = lambda x: x * 2625.4996295540054
to_kcalmol = lambda x: x * 630.52344609846432
to_ev = lambda x: x * 27.211384565719481

# Mass

to_unified = lambda x: x * 0.0005485799109814269

# Charge

to_coulomb = lambda x: x * 1.60217653e-19

# Degrees

to_degrees = lambda x: x * 57.295779513082323

# In a dictionary

to_unit = {
    ATOMARY: lambda x: x,
    ANGSTROM: to_angstrom,
    NANOMETER: to_nanometer,
    JOULE: to_joule,
    CALORIE: to_calorie,
    KJMOL: to_kjmol,
    KCALMOL: to_kcalmol,
    EV: to_ev,
    UNIFIED: to_unified,
    COULOMB: to_coulomb,
    RADIALS: lambda x: x,
    DEGREES: to_degrees
}
