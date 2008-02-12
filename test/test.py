# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


import sys, glob, os
retcode = os.system("(cd ..; python setup.py build)")
if retcode != 0: sys.exit(retcode)
sys.path.insert(0, glob.glob("../build/lib*")[0])

import os
if not os.path.exists("output"):
    os.mkdir("output")

from molmod import context
context.share_path = "../share/"

import unittest

from binning import *
from graphs import *
from molecules import *
from molecular_graphs import *
from internal_coordinates import *
from pairff import *
from flexff import *
from data import *
from unit_cell import *
from clusters import *
from transformations import *
from environments import *
from descriptors import *
#from similarity import *
from tunneling import *
from quaternions import *
from io import *

unittest.main()



