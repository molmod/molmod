# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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
# --
"""Computation of molecular volumes with monte carlo methods

This approach is slow, but robust and simple to implement.
"""

from molmod.ext import monte_carlo_volumes
from molmod.periodic import periodic
from molmod.units import angstrom

import numpy


__all__ = ["estimate_volumes", "estimate_radius"]


def estimate_volumes(molecule, probe_radius=2.65, max_error=6.0, bigbox=True):
    """Returns a Monte Carlo estimates of molecular volumes.

       Arguments:
         molecule  --  A molecule object.
         probe_radius  --  A representative radius of the solvent
         max_error  --  An upper limit for the error on the estimates

       Returns:
         vdw_volume  --  Estimate of the Van Der Waals volume
         sas_volume  --  Estimate of the volume within the Surface Accessible
                         surface
         ses_volume  --  Estimate of the volume within the Solvent Excluded
                         surface

       It is most efficient to rotate the molecule so that its bounding box
       (including the Van Der Waals spheres) in Cartesian coordinates has a
       minimal volume.
    """

    counts = numpy.zeros(5, numpy.int64)
    radii = numpy.array([periodic[number].vdw_radius for number in molecule.numbers])

    def mean_error(box_volume, count, total):
        """compute the mean and the error on the mean of a binomial distro"""
        fraction = float(count)/total
        mean = box_volume*fraction
        error = box_volume*numpy.sqrt(fraction*(1-fraction)/total)
        return mean, error

    error = 2*max_error
    while error > max_error:
        box_volume = monte_carlo_volumes(probe_radius, molecule.coordinates, radii, 10000*molecule.size, 500, bigbox, counts)
        vdw_volume, vdw_error = mean_error(box_volume, counts[0], counts[-1])
        sas_volume, sas_error = mean_error(box_volume, counts[1], counts[-1])
        ses_volume, ses_error = mean_error(box_volume, counts[2], counts[-1])
        error = max(vdw_error, sas_error, ses_error)

    return vdw_volume, sas_volume, ses_volume, error


def estimate_radius(mol, max_error=0.02):
    """Estimates the effective radius of a solvent molecule.

       Arguments:
         mol  --  The solvent molecule object
         max_error  --  A convergence threshold for the self-consistent estimate
                        of the radius

       This function estimates the effective solvent radius in a self-consistent
       way. The initial guess for the radius is 1.0*angstrom. Based on this
       radius, the solvent excluded volume is computed. Then a new radius is
       computed as the radius of a sphere with a volume equal to the solvent
       excluded volume. This is repeated until the change in radii between two
       iterations and the sampling error is smaller than the given threshold.
    """
    radius = 1.0*angstrom
    while True:
        max_vol_error = min(4*numpy.pi*radius*radius*max_error, 2.65)
        ses_volume, vol_error = estimate_volumes(mol, radius, max_vol_error, False)[2:]
        new_radius = (3*ses_volume/4/numpy.pi)**(1.0/3)
        error = max(vol_error/(4*numpy.pi*new_radius*new_radius), abs(new_radius-radius))
        radius = new_radius
        #print radius/angstrom, error/angstrom
        if error < max_error:
            return radius, error
