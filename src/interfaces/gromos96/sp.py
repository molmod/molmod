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

from pychem.interfaces.base import TemplateJob
from pychem.units import to_nanometer

import os, copy, glob
import Numeric


__all__ = ["Gromos96SP"]


class Gromos96SP(TemplateJob):
    output_filenames = ["FORCES", "sxmd2.dat", "rxmd2.dat", "omd.out", "ogmt.out", "mta1.dat"]
    template_filenames = ["igmt.dat", "sxmd1.dat", "imd.dat"]
    template = "gromos96_sp"
    scripts = ["jgmt.sh", "jmd.sh"]
    
    def __init__(self, prefix, title, input_molecule, parameters, topology, box_size, gromos_root="/usr/local/gromos96", output_parser=None):
        mapping = dict((key, str(value)) for key, value in parameters.iteritems())
        mapping['topology'] = topology
        mapping['box_size'] = str(to_nanometer(box_size)) # cubic boxes only
        mapping['coordinates'] = "\n".join([" "*24 + " %f %f %f" % tuple(coordinate) for coordinate in to_nanometer(input_molecule.coordinates)])
        mapping['gromos_root'] = gromos_root
        TemplateJob.__init__(self, prefix, title, mapping)
        self.output_parser = output_parser

    def determine_completed(self):
        self.completed = (
            os.system("grep \"PROMD WRITING FINAL CONFIGURATION\" %somd.out &> /dev/null" % self.directory) == 0 and
            os.system("grep \"WRTOPO: OK!\" %sogmt.out &> /dev/null" % self.directory) == 0
        )
    
    def read_output(self):
        if self.output_parser != None:
            slash_pos = self.prefix.rfind("/")
            directory = self.prefix[:slash_pos]
            prefix = self.prefix[slash_pos+1:]
            self.__dict__.update(self.output_parser.parse(directory, prefix))        
