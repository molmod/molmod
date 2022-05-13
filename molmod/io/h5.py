# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2019 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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
"""Tools for reading and writing h5 files"""


from __future__ import print_function


from __future__ import division

from builtins import range
import numpy as np
import h5py as h5

from molmod.io.common import SlicedReader, FileFormatError
from molmod.periodic import periodic
from molmod.molecules import Molecule
from molmod.units import angstrom


__all__ = ["h5Reader"]


class h5Reader(SlicedReader):
    """A reader for h5 trajectory files
    """

    def __init__(self, f, subdirectory=None, sub=slice(None)):
        """Initialize an h5 reader

           Arguments:
            | ``f``  --  a filename or a file-like object
            | ``directory`` -- a subdirectory in h5 where the data is stored

           Optional arguments:
            | ``sub``  --  a slice indicating which frames to read/skip
            | ``file_unit``  --  the conversion constant to convert data into atomic
                                  units [default=angstrom]

           After initialization, the following attributes are defined:
            | ``symbols``  --  The atom symbols
            | ``numbers``  --  The atom numbers
        """
        SlicedReader.__init__(self, f, sub)
        self.f = f
        self.subdirectory = subdirectory


    def read_subdirectory(self):
        f_h5_1 = h5.File(self.f, mode = 'r')
        data = np.array(f_h5_1['%s'%self.subdirectory])
        return data

    def read_size(self):
        raise NotImplementedError

    def _read_frame(self):
        raise NotImplementedError

    def _skip_frame(self):
        raise NotImplementedError

    def get_first_molecule(self):
        raise NotImplementedError
