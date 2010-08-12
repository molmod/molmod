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
"""Runtime context information for shared files

   This module defines a ``context`` object of the class :class:`Context`. The
   instance can be used to retrieve shared files::

   >>> from molmod import context
   >>> context.get_share_filename("periodic.csv")
"""


import os


__all__ = ["Context", "context"]


class Context(object):
    """Global variable to find and use share directory"""
    def __init__(self):
        """Initialize the Context object

           This is done once when importing a molmod module. There is no need
           to do this manually a second time.
        """
        # find the data files
        fn_datadir = os.path.join(os.path.dirname(__file__), "datadir.txt")
        if os.path.isfile(fn_datadir):
            f = file(fn_datadir)
            datadir = f.readline().strip()
            f.close()
            self.share_dir = os.path.join(datadir, "share", "molmod")
        else:
            # When running from the build directory for the tests.
            self.share_dir = "../share"
        if not os.path.isdir(self.share_dir):
            raise RuntimeError("Share dir '%s' does not exist." % self.share_dir)

    def get_share_filename(self, filename):
        """Retrieve the full path for a given filename in the share folder

           Argument:
            | filename  --  A file from the share directory (without the
                            ``share/`` prefix)
        """
        result = os.path.join(self.share_dir, filename)
        if not os.path.isfile(result):
            raise ValueError("Data file '%s' not found." % result)
        return result


context = Context()
