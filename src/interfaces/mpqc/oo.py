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

from pychem.interfaces.base import Job
from pychem.interfaces.mpqc.keyval import KeyValObject
from pychem.moldata import periodic

import os, glob


class OOMpqcJob(Job):
    def __init__(self, prefix, title, keyval, output_parser):
        self.keyval = keyval
        self.output_parser = output_parser
        Job.__init__(self, prefix, title)
        
    def write_input(self, f):
        print >> f, "%% %s" % self.title
        self.keyval.write_stream(f)

    def external_command(self):
        return "mpqc -o %s.out %s.in" % (self.filename, self.filename)
        
    def remove_temporary_files(self):
        for temp_filename in glob.glob("%s*.tmp" % self.filename):
            os.remove(temp_filename)

    def determine_completed(self):
        self.completed = (os.system("grep \"End Time\" %s.out &> /dev/null" % self.filename) == 0)
        
    def read_output(self):
        self.__dict__.update(self.output_parser.parse(self.filename))
