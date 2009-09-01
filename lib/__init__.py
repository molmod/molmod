# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import sys, os


__all__ = ["context"]


class Error(Exception):
    pass


class Context(object):
    def __init__(self):
        # find the data files
        fn_datadir = os.path.join(os.path.dirname(__file__), "datadir.txt")
        if os.path.isfile(fn_datadir):
            f = file(fn_datadir)
            datadir = f.readline().strip()
            f.close()
            self.share_dir = os.path.join(datadir, "share", "molmod")
        else:
            self.share_dir = "../share" # When running from the build directory for the tests.
        if not os.path.isdir(self.share_dir):
            raise RuntimeError("Share dir '%s' does not exist." % self.share_dir)

    def get_share_filename(self, filename):
        result = os.path.join(self.share_dir, filename)
        if not os.path.isfile(result):
            raise ValueError("Data file '%s' not found." % result)
        return result


context = Context()


from clusters import *
from graphs import *
from ic import *
from io import *
from minimizer import *
from molecules import *
from molecular_graphs import *
from quaternions import *
from randomize import *
from similarity import *
from toyff import *
from transformations import *
from unit_cells import *
from units import *
from vectors import *
from volume import *
from zmatrix import *


