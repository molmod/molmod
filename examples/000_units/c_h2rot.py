#!/usr/bin/env python
# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
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
#--
#!/usr/bin/env python

from molmod import *
from numpy import exp, pi

# The parameters
m = 1.00794*amu
d = 74*picometer
temp = 300

# The moment of inertia
imom = m*d**2/2
# The rotational temperature:
rot_temp = 1.0/(2*imom*boltzmann)
# A function that computes a term in the rotational partition function:
def q_rot_term(j):
    # j is the index of the term
    return (2*j+1)*exp(-rot_temp*j*(j+1)/temp)

# Approximate the partition function with 10 terms. This is good enough for H2.
q_rot = 0.0
for i in xrange(10):
    q_rot += q_rot_term(i)

# Compute the probability of the state j=0
prob = q_rot_term(0)/q_rot

print "The rotational temperature [Kelvin]:", rot_temp
print "Probability [%]:", prob*100
