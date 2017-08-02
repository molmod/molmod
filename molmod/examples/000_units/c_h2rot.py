#!/usr/bin/env python

from __future__ import print_function

from molmod import *
from numpy import exp, pi

# The parameters
m = 1.00794*amu
d = 74*picometer
temp = 300

# The moment of inertia
imom = m*d**2/2
# The rotational temperature:
rot_temp = 1/(2*imom*boltzmann)
# A function that computes a term in the rotational partition function:
def q_rot_term(j):
    # j is the index of the term
    return (2*j+1)*exp(-rot_temp*j*(j+1)/temp)

# Approximate the partition function with 10 terms. This is good enough for H2.
q_rot = 0.0
for i in range(10):
    q_rot += q_rot_term(i)

# Compute the probability of the state j=0
prob = q_rot_term(0)/q_rot

print("The rotational temperature [Kelvin]:", rot_temp)
print("Probability [%]:", prob*100)
