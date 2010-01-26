#!/usr/bin/python
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


import sys, os

if '-i' in sys.argv:
    # use the installed library for testing
    sys.argv.remove('-i')
else:
    import glob
    if '-c' in sys.argv:
        sys.argv.remove('-c')
        os.system("cd ../; rm -rf build; cd ext; rm -rf build")
    retcode = os.system("(cd ..; python setup.py build)")
    if retcode != 0: sys.exit(retcode)
    retcode = os.system("(cd ../ext; python setup.py build; cp -avr build/lib*/* ../build/lib*/)")
    if retcode != 0: sys.exit(retcode)
    sys.path.insert(0, glob.glob("../build/lib*")[0])

if not os.path.exists("output"):
    os.mkdir("output")


import unittest

from binning import *
from clusters import *
from data import *
from graphs import *
from ic import *
from io import *
from minimizer import *
from molecules import *
from molecular_graphs import *
from pairff import *
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

unittest.main()


