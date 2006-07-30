# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of MolMod.
# 
# MolMod is free software; you can redistribute it and/or
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

from molmod.interfaces.base import IOJob
from molmod.interfaces.mpqc.keyval import KeyValWriter

import os, glob


class OOMpqcJob(IOJob):
    binary = "mpqc"

    def __init__(self, prefix, title, keyval):
        self.keyval = keyval
        IOJob.__init__(self, prefix, title)
        
    def write_input(self, f):
        print >> f, "%% %s" % self.title
        KeyValWriter(f, self.keyval)

    def external_command(self):
        return "%s -o %s%s %s%s" % (self.binary, self.prefix, self.output_extension, self.prefix, self.input_extension)
        
    def remove_temporary_files(self):
        for temp_filename in glob.glob("%s*.tmp" % self.prefix):
            os.remove(temp_filename)

    def determine_completed(self):
        self.completed = (os.system("grep \"End Time\" %s%s &> /dev/null" % (self.prefix, self.output_extension)) == 0)
        
