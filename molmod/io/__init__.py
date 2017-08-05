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
# --
"""Collection of routines to handle computational chemistry file formats

   It is not the intention to duplicate the work done in the OpenBabel project.
   Several file formats below are rather specific for molecular dynamics
   simulations and not all of them deal with representing molecular systems. The
   scope is more generic. The selection of supported formats is purely driven
   by the convenience of the MolMod developers.
"""

from molmod.io.atrj import *
from molmod.io.chk import *
from molmod.io.cml import *
from molmod.io.common import *
from molmod.io.cp2k import *
from molmod.io.cpmd import *
from molmod.io.crystal import *
from molmod.io.cube import *
from molmod.io.dlpoly import *
from molmod.io.fchk import *
from molmod.io.gamess import *
from molmod.io.gromacs import *
from molmod.io.lammps import *
from molmod.io.number_state import *
from molmod.io.pdb import *
from molmod.io.psf import *
from molmod.io.sdf import *
from molmod.io.xyz import *
