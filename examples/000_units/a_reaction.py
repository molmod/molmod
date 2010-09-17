#!/usr/bin/env python

from molmod import *

energy_react = -157.31456213
energy_prod = -157.31397873

reaction_energy = energy_prod - energy_react

print "Reaction energy [kJ mol^-1]:", (reaction_energy/kjmol)
