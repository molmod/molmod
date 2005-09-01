# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 
# --


import init_files

from pychem import context
context.share_path = "pychem/"

import unittest
suite = unittest.TestSuite()


#import binning
#suite.addTest(binning.suite)

#import graphs
#suite.addTest(graphs.suite)

#import interfaces
#suite.addTest(interfaces.suite)

import molecular_graphs
suite.addTest(molecular_graphs.suite)

test_result = unittest.TextTestRunner(verbosity=2).run(suite)

if len(test_result.failures) > 0 or len(test_result.errors) > 0:
    import sys
    sys.exit(1)
