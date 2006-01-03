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


from pychem.interfaces.output_parsers import FileParser, MultiLineParser

import re, Numeric

class ForcesParser(MultiLineParser):
    filename = "FORCES"
    extension = False

    def reset(self):
        MultiLineParser.reset(self)
        self.forces = None
        
    def start_collecting(self):
        if self.forces != None:
            self.active = False
        else:
            self.forces = []

    def collect(self, line):
        self.forces.append(float(line))

    def stop_collecting(self):
        self.forces = Numeric.array(self.forces)
        self.forces.shape = (-1, 3)
         
    def result(self):
        return self.forces


class Forces1Parser(ForcesParser):
    def __init__(self, label="forces", condition=None):
        MultiLineParser.__init__(self, label, 
            activator=re.compile(r"^ ---$"), 
            deactivator=re.compile(r"^    $"), 
            condition=condition
        )


class Forces2Parser(ForcesParser):
    def __init__(self, label="forces", condition=None):
        MultiLineParser.__init__(self, label, 
            activator=re.compile(r"^    $"), 
            deactivator=re.compile(r"^ ---$"), 
            condition=condition
        )
