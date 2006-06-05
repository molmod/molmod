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

# Constants that depict measures and units

class Counter:
    def __init__(self):
        self.counter = 0
        
    def __call__(self):
        self.counter += 1
        return self.counter
        
measure_counter = Counter()
unit_counter = Counter()

ATOMARY = measure_counter()

LENGTH = measure_counter()
ANGSTROM = unit_counter()
NANOMETER = unit_counter()
METER = unit_counter()

ENERGY = measure_counter()
JOULE = unit_counter()
CALORIE = unit_counter()
KJMOL = unit_counter()
KCALMOL = unit_counter()
EV = unit_counter()

MASS = measure_counter()
UNIFIED = unit_counter()
KG = unit_counter()

CHARGE = measure_counter()
COULOMB = unit_counter()

ANGLE = measure_counter()
RADIALS = unit_counter()
DEGREES = unit_counter()

TIME = measure_counter()
SECONDS = unit_counter()


suffices = {
    ATOMARY: "a.u.",
    ANGSTROM: "A",
    NANOMETER: "nm",
    METER: "m",
    JOULE: "J",
    CALORIE: "cal",
    KJMOL: "kJ/mol",
    KCALMOL: "kcal/mol",
    EV: "eV",
    UNIFIED: "u",
    KG: "kg", 
    COULOMB: "C",
    RADIALS: "rad",
    DEGREES: "°",
    SECONDS: "s"
}

tex_suffices = {
    ATOMARY: r"a.u.",
    ANGSTROM: r"\AA",
    NANOMETER: r"nm",
    METER: r"m",
    JOULE: r"J",
    CALORIE: r"cal",
    KJMOL: r"\frac{kJ}{mol}",
    KCALMOL: r"\frac{kcal}{mol}",
    EV: r"eV",
    UNIFIED: r"u",
    KG: r"kg", 
    COULOMB: r"C",
    RADIALS: r"rad",
    DEGREES: r"^{\circ}",
    SECONDS: r"s"
}

measures = {
    LENGTH: [ATOMARY, ANGSTROM, NANOMETER, METER],
    ENERGY: [ATOMARY, JOULE, CALORIE, KJMOL, KCALMOL, EV],
    MASS: [ATOMARY, UNIFIED, KG],
    CHARGE: [ATOMARY, COULOMB],
    ANGLE: [RADIALS, DEGREES],
    TIME: [SECONDS]
}

measure_names = {
    LENGTH: "Length",
    ENERGY: "Energy",
    MASS: "Mass",
    CHARGE: "Charge",
    ANGLE: "Angle",
    TIME: "Time"
}


# Usefull constants

avogadro = 6.0221415e23

# Length

angstrom = 1.8897261249
nanometer = 18.897261249
meter = 1.8897261249e10

def to_angstrom(x): return x * 0.5291772108
def to_nanometer(x): return x * 0.05291772108
def to_meter(x): return x * 0.5291772108e-10

def from_angstrom(x): return x * 1.8897261249
def from_nanometer(x): return x * 18.897261249
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

degrees = 0.017453292519943295

def to_degrees(x): return x * 57.295779513082323

def from_degrees(x): return x * 0.017453292519943295

# Time

second = 41341373336561368.0

def to_seconds(x): return x * 2.418884326505e-17

def from_seconds(x): return x * 41341373336561368.0

# In a dictionary

def unity(x): return x

unit = {
    ATOMARY: 1,
    ANGSTROM: angstrom,
    NANOMETER: nanometer,
    METER: meter,
    JOULE: joule,
    CALORIE: calorie,
    KJMOL: kjmol,
    KCALMOL: kcalmol,
    EV: ev,
    UNIFIED: unified,
    COULOMB: coulomb,
    RADIALS: 1,
    DEGREES: degrees
}

to_unit = {
    ATOMARY: unity,
    ANGSTROM: to_angstrom,
    NANOMETER: to_nanometer,
    METER: to_meter,
    JOULE: to_joule,
    CALORIE: to_calorie,
    KJMOL: to_kjmol,
    KCALMOL: to_kcalmol,
    EV: to_ev,
    UNIFIED: to_unified,
    COULOMB: to_coulomb,
    RADIALS: unity,
    DEGREES: to_degrees
}

from_unit = {
    ATOMARY: unity,
    ANGSTROM: from_angstrom,
    NANOMETER: from_nanometer,
    METER: from_meter,
    JOULE: from_joule,
    CALORIE: from_calorie,
    KJMOL: from_kjmol,
    KCALMOL: from_kcalmol,
    EV: from_ev,
    UNIFIED: from_unified,
    COULOMB: from_coulomb,
    RADIALS: unity,
    DEGREES: from_degrees
}
