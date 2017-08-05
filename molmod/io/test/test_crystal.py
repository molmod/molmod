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


from __future__ import print_function

import pkg_resources

from molmod.io import *
from molmod import *


def test_crystal_quartz():
    cryout = CrystalAPIOut(pkg_resources.resource_filename(__name__, "../../data/test/crystal_api.out"))
    # unit cell
    assert abs(cryout.unit_cell.matrix[0,0] - 8.041028888) < 1e-8
    assert cryout.unit_cell.matrix[2,0] == 0.0
    assert abs(cryout.unit_cell.matrix[1,0] - -4.642490193) < 1e-8
    assert cryout.unit_cell.active.sum() == 3
    # molecule
    assert cryout.mol.size == 9
    assert cryout.mol.numbers[4] == 8
    assert abs(cryout.mol.coordinates[5,1] - -3.167106809804E+00) < 1e-8
    # basis
    assert len(cryout.basisset) == 2
    assert len(cryout.basisset['SI']) == 6
    assert len(cryout.basisset['O']) == 4
    assert len(cryout.basisset['SI'][0]) == 2
    assert len(cryout.basisset['SI'][5]) == 2
    assert len(cryout.basisset['O'][0]) == 2
    assert len(cryout.basisset['O'][3]) == 2
    assert cryout.basisset['SI'][0][0] == 'S'
    assert cryout.basisset['SI'][1][0] == 'SP'
    assert cryout.basisset['SI'][5][0] == 'D'
    assert cryout.basisset['O'][0][0] == 'S'
    assert cryout.basisset['O'][1][0] == 'SP'
    assert cryout.basisset['O'][3][0] == 'D'
    assert len(cryout.basisset['SI'][0][1]) == 8
    assert len(cryout.basisset['SI'][1][1]) == 6
    assert len(cryout.basisset['SI'][5][1]) == 1
    assert len(cryout.basisset['O'][0][1]) == 6
    assert len(cryout.basisset['O'][1][1]) == 3
    assert len(cryout.basisset['O'][3][1]) == 1
    assert len(cryout.basisset['SI'][0][1][0]) == 2
    assert len(cryout.basisset['SI'][1][1][0]) == 3
    assert len(cryout.basisset['SI'][5][1][0]) == 2
    assert len(cryout.basisset['O'][0][1][0]) == 2
    assert len(cryout.basisset['O'][1][1][0]) == 3
    assert len(cryout.basisset['O'][3][1][0]) == 2
    assert len(cryout.basisset['SI'][0][1][7]) == 2
    assert len(cryout.basisset['SI'][1][1][5]) == 3
    assert len(cryout.basisset['O'][0][1][5]) == 2
    assert len(cryout.basisset['O'][1][1][2]) == 3
    assert abs(cryout.basisset['SI'][0][1][0][0] - 8.765E+04) < 1e+0
    assert abs(cryout.basisset['SI'][1][1][0][2] - 2.360E+01) < 1e-3
    assert abs(cryout.basisset['SI'][5][1][0][1] - 1.637E+00) < 1e-4
    assert abs(cryout.basisset['O'][0][1][0][1] - 2.023E+00) < 1e-4
    assert abs(cryout.basisset['O'][1][1][0][1] - -1.503E+00) < 1e-4
    assert abs(cryout.basisset['O'][3][1][0][0] - 8.000E-01) < 1e-5
    assert abs(cryout.basisset['SI'][0][1][7][1] - 1.683E+00) < 1e-4
    assert abs(cryout.basisset['SI'][1][1][5][0] - 7.365E-01) < 1e-5
    assert abs(cryout.basisset['O'][0][1][5][0] - 5.800E+00) < 1e-1
    assert abs(cryout.basisset['O'][1][1][2][1] - 1.980E+00) < 1e-4
    # density matrix
    print(cryout.density_matrix[:5,:5])
    assert cryout.density_matrix[0,0] == 0.0
    assert abs(cryout.density_matrix[2,2] - 2.1646E-02) < 1e-5
    print(cryout.density_matrix[1,4])
    assert abs(cryout.density_matrix[4,1] - -1.9742E-01) < 1e-4
