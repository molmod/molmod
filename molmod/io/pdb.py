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
"""Basic support for the PDB format"""


from __future__ import print_function, division

from builtins import range
import numpy as np

from molmod.periodic import periodic
from molmod.units import angstrom
from molmod.molecules import Molecule
from molmod.io.common import FileFormatError


__all__ = ["load_pdb", "dump_pdb"]


def dump_pdb(filename, molecule, atomnames=None, resnames=None, chain_ids=None, occupancies=None, betas=None):
    """Writes a single molecule to a pdb file.

       This function is based on the pdb file specification:
       http://www.wwpdb.org/documentation/format32/sect9.html
       For convenience, the relevant table is copied and the character indexes are
       transformed to C-style (starting from zero)

       =======        ============  ==========   ==========================================
       COLUMNS        DATA  TYPE    FIELD        DEFINITION
       =======        ============  ==========   ==========================================
        0 -  5        Record name   "ATOM  "
        6 - 10        Integer       serial       Atom  serial number.
       12 - 15        Atom          name         Atom name.
       16             Character     altLoc       Alternate location indicator.
       17 - 19        Residue name  resName      Residue name.
       21             Character     chainID      Chain identifier.
       22 - 25        Integer       resSeq       Residue sequence number.
       26             AChar         iCode        Code for insertion of residues.
       30 - 37        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
       38 - 45        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
       46 - 53        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
       54 - 59        Real(6.2)     occupancy    Occupancy.
       60 - 65        Real(6.2)     tempFactor   Temperature  factor.
       76 - 77        LString(2)    element      Element symbol, right-justified.
       78 - 79        LString(2)    charge       Charge  on the atom.
       =======        ============  ==========   ==========================================
    """

    with open(filename, "w") as f:
        res_id = 1
        old_resname = None

        for i in range(molecule.size):
            symbol = periodic[molecule.numbers[i]].symbol
            if atomnames is None:
                atomname = symbol
            else:
                atomname = atomnames[i]
            if resnames is None:
                resname = "OXO"
            else:
                resname = resnames[i]
            if resname != old_resname:
                res_id += 1
            if chain_ids is None:
                chain_id = "A"
            else:
                chain_id = chain_ids[i]
            if occupancies is None:
                occupancy = 1.0
            else:
                occupancy = occupancies[i]
            if betas is None:
                beta = 1.0
            else:
                beta = betas[i]

            print("ATOM   %4i  %3s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  " % (
                i+1, atomname.ljust(3), resname.ljust(3), chain_id, res_id,
                molecule.coordinates[i, 0]/angstrom,
                molecule.coordinates[i, 1]/angstrom,
                molecule.coordinates[i, 2]/angstrom,
                occupancy, beta, symbol.ljust(2)
            ), file=f)
            old_resname = resname


def load_pdb(filename):
    """Loads a single molecule from a pdb file.

       This function does support only a small fragment from the pdb specification.
       It assumes that there is only one molecular geometry in the pdb file.
    """
    with open(filename) as f:
        numbers = []
        coordinates = []
        occupancies = []
        betas = []
        for line in f:
            if line.startswith("ATOM"):
                symbol = line[76:78].strip()
                numbers.append(periodic[symbol].number)
                coordinates.append([float(line[30:38])*angstrom, float(line[38:46])*angstrom, float(line[46:54])*angstrom])
                occupancies.append(float(line[54:60]))
                betas.append(float(line[60:66]))
    if len(numbers) > 0:
        molecule = Molecule(numbers, coordinates)
        molecule.occupancies = np.array(occupancies)
        molecule.betas = np.array(betas)
        return molecule
    else:
        raise FileFormatError("No molecule found in pdb file %s" % filename)
