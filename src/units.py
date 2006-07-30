# -*- coding: utf-8 -*-
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
#
# Source of the conversion values: The NIST Reference on Constants, Units, 
# and Uncertainty (http://physics.nist.gov/cuu/Constants/)
# Except for calorie: 1 calorie = 4.184 Joules
#
#
# Naming conventions in this module:
#    * to_unit is a function that converts a value in internal units to a
#      value in external units
#    * unit is the value of one external unit in internal units
#    * from_unit is a function that does the inverse of to_unit
#  The internal units are the atomic units. External units can be anything else.


# Constants that depict measures and units.
        
ATOMARY = 0

LENGTH = 0
ANGSTROM = 1
NANOMETER = 2
PICOMETER = 13

ENERGY = 1
KJMOL = 3
KCALMOL = 4
EV = 5

MASS = 2
UNIFIED = 6

CHARGE = 3

ANGLE = 4
RADIAL = 7
DEGREE = 8

TIME = 5
NANOSECOND = 9
PICOSECOND = 10
FEMTOSECOND = 11

DIPOLE = 6
DEBYE = 12

measures = [LENGTH, ENERGY, MASS, CHARGE, ANGLE, TIME, DIPOLE]
assert len(measures) == len(set(measures)), "Some measures have the same id."

units = [
    ATOMARY, ANGSTROM, NANOMETER, PICOMETER, KJMOL, KCALMOL, EV, UNIFIED, 
    RADIAL, DEGREE, NANOSECOND, PICOSECOND, FEMTOSECOND, DEBYE
]
assert len(units) == len(set(units)), "Some units have the same id."


suffices = {
    ATOMARY: "a.u.",
    ANGSTROM: "A",
    NANOMETER: "nm",
    PICOMETER: "pm",
    KJMOL: "kJ/mol",
    KCALMOL: "kcal/mol",
    EV: "eV",
    UNIFIED: "u",
    RADIAL: "rad",
    DEGREE: "°",
    NANOSECOND: "ns",
    PICOSECOND: "ps",
    FEMTOSECOND: "fs",
    DEBYE: "D",
}
assert len(units) == len(suffices), "Some units don't have a suffix."

by_suffices = {
    "a.u.": ATOMARY,
    "A": ANGSTROM,
    "nm": NANOMETER,
    "pm": PICOMETER,
    "kJ/mol": KJMOL,
    "kcal/mol": KCALMOL,
    "eV": EV,
    "u": UNIFIED,
    "rad": RADIAL,
    "°": DEGREE,
    "ns": NANOSECOND,
    "ps": PICOSECOND,
    "fs": FEMTOSECOND,
    "D": DEBYE,
}
assert len(units) == len(by_suffices), "Some units don't have a suffix."

tex_suffices = {
    ATOMARY: r"a.u.",
    ANGSTROM: r"\AA",
    NANOMETER: r"nm",
    PICOMETER: r"pm",
    KJMOL: r"\frac{kJ}{mol}",
    KCALMOL: r"\frac{kcal}{mol}",
    EV: r"eV",
    UNIFIED: r"u",
    RADIAL: r"rad",
    DEGREE: r"^{\circ}",
    NANOSECOND: r"ns",
    PICOSECOND: r"ps",
    FEMTOSECOND: r"fs",
    DEBYE: r"D"
}
assert len(units) == len(suffices), "Some units don't have a tex suffix."

units_by_measure = {
    LENGTH: [ATOMARY, ANGSTROM, NANOMETER, PICOMETER],
    ENERGY: [ATOMARY, KJMOL, KCALMOL, EV],
    MASS: [ATOMARY, UNIFIED],
    CHARGE: [ATOMARY],
    ANGLE: [RADIAL, DEGREE],
    TIME: [ATOMARY, NANOSECOND, PICOSECOND, FEMTOSECOND],
    DIPOLE: [ATOMARY, DEBYE],
}
assert len(measures) == len(units_by_measure), "Some measures don't have units."

measure_names = {
    LENGTH: "Length",
    ENERGY: "Energy",
    MASS: "Mass",
    CHARGE: "Charge",
    ANGLE: "Angle",
    TIME: "Time",
    DIPOLE: "Dipole",
}
assert len(measures) == len(measure_names), "Some measures don't have a name."

# some sanity checks:


# Length

angstrom = 1.8897261249
nanometer = 18.897261249
picometer = 0.018897261249
meter = 1.8897261249e10

def to_angstrom(x): return x * 0.5291772108
def to_nanometer(x): return x * 0.05291772108
def to_picometer(x): return x * 52.91772108
def to_meter(x): return x * 0.5291772108e-10

def from_angstrom(x): return x * 1.8897261249
def from_nanometer(x): return x * 18.897261249
def from_picometer(x): return x * 0.018897261249
def from_meter(x): return x * 1.8897261249e10


# Energy

joule = 2.2937125689189235e17
calorie = 9.596893388356777e17
kjmol = 0.00038087988615327677
kcalmol = 0.0015859838459422444
ev = 0.036749324444879071

def to_joule(x): return x * 4.35974417e-18
def to_calorie(x): return x * 1.0420038647227532e-18
def to_kjmol(x): return x * 2625.4996295540054
def to_kcalmol(x): return x * 630.52344609846432
def to_ev(x): return x * 27.211384565719481

def from_joule(x): return x * 2.2937125689189235e17
def from_calorie(x): return x * 9.596893388356777e17
def from_kjmol(x): return x * 0.00038087988615327677
def from_kcalmol(x): return x * 0.0015859838459422444
def from_ev(x): return x * 0.036749324444879071


# Mass

unified = 1822.8884798405547
kg = 1.0977692384992151e30

def to_unified(x): return x * 0.0005485799109814269
def to_kg(x): return x * 9.1093826e-31

def from_unified(x): return x * 1822.8884798405547
def from_kg(x): return x * 1.0977692384992151e30

# Charge

coulomb = 6.2415094796077179e18

def to_coulomb(x): return x * 1.60217653e-19

def from_coulomb(x): return x * 6.2415094796077179e18

# Angles

degree = 0.017453292519943295

def to_degree(x): return x * 57.295779513082323

def from_degree(x): return x * 0.017453292519943295

# Time

second = 41341373336561368.0
nanosecond = 41341373.336561368
picosecond = 41341.373336561368
femtosecond = 41.341373336561368

def to_second(x): return x * 2.418884326505e-17
def to_nanosecond(x): return x * 2.418884326505e-8
def to_picosecond(x): return x * 2.418884326505e-5
def to_femtosecond(x): return x * 2.418884326505e-2

def from_second(x): return x * 41341373336561368.0
def from_nanosecond(x): return x * 41341373.336561368
def from_picosecond(x): return x * 41341.373336561368
def from_femtosecond(x): return x * 41.341373336561368

# dipole

debye = 0.39343018283144093 # = 3.33564e-30*coulomb*meter

def to_debye(x): return x * 2.541747033242832

def from_debye(x): return x * 0.39343018283144093

# In a dictionary

def unity(x): return x

unit = {
    ATOMARY: 1,
    ANGSTROM: angstrom,
    NANOMETER: nanometer,
    PICOMETER: picometer,
    KJMOL: kjmol,
    KCALMOL: kcalmol,
    EV: ev,
    UNIFIED: unified,
    RADIAL: 1,
    DEGREE: degree,
    NANOSECOND: nanosecond,
    PICOSECOND: picosecond,
    FEMTOSECOND: femtosecond,
    DEBYE: debye,
}
assert len(units) == len(unit), "Some units don't have a conversion."

to_unit = {
    ATOMARY: unity,
    ANGSTROM: to_angstrom,
    NANOMETER: to_nanometer,
    PICOMETER: to_picometer,
    KJMOL: to_kjmol,
    KCALMOL: to_kcalmol,
    EV: to_ev,
    UNIFIED: to_unified,
    RADIAL: unity,
    DEGREE: to_degree,
    NANOSECOND: to_nanosecond,
    PICOSECOND: to_picosecond,
    FEMTOSECOND: to_femtosecond,
    DEBYE: to_debye,
}
assert len(units) == len(to_unit), "Some units don't have a to_ function."

from_unit = {
    ATOMARY: unity,
    ANGSTROM: from_angstrom,
    NANOMETER: from_nanometer,
    PICOMETER: from_picometer,
    KJMOL: from_kjmol,
    KCALMOL: from_kcalmol,
    EV: from_ev,
    UNIFIED: from_unified,
    RADIAL: unity,
    DEGREE: from_degree,
    NANOSECOND: from_nanosecond,
    PICOSECOND: from_picosecond,
    FEMTOSECOND: from_femtosecond,
    DEBYE: from_debye,
}
assert len(units) == len(from_unit), "Some units don't have a from_ function."

