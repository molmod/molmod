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
from pychem.moldata import periodic
from pychem.molecules import Molecule

import re, Numeric

class NumStepsParser(FileParser):
    extension = ".out"

    def __init__(self, label='num_steps', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"MAXIMUM NUMBER OF STEPS:\s+(?P<num_steps>\S+)\s+STEPS")
        
    def reset(self):
        self.num_steps = None
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.num_steps = float(match.group("num_steps"))
        
    def result(self):
        return self.num_steps


class NumEveryParser(FileParser):
    extension = ".out"

    def __init__(self, label='num_every', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"PRINT INTERMEDIATE RESULTS EVERY\s+(?P<num_every>\S+)\s+STEPS")
        
    def reset(self):
        self.num_every = None
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.num_every = float(match.group("num_every"))
        
    def result(self):
        return self.num_every


class ElementsParser(MultiLineParser):
    extension = ".out"

    def __init__(self, label='elements', condition=None):
        MultiLineParser.__init__(
            self,
            label,
            re.compile(r"\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* ATOMS \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*"),
            re.compile(r"\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*"),
            condition
        )
        self.re = re.compile(r"\s+\d+\s+(?P<symbol>\S+)\s+\S+\s+\S+\s+\S+")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.done = False
        
    def start_collecting(self):
        self.elements = []

    def collect(self, line):
        if not self.done:
            match = self.re.search(line)
            if match != None:
                self.elements.append(periodic.symbol_lookup(match.group("symbol")))
        
    def stop_collecting(self):
        self.done = True

    def result(self):
        return self.elements


class TrajectoryParserMixin(object):
    extension = ".out"
    
    def __init__(self, num_steps_parser, num_every_parser):
        self.num_steps_parser = num_steps_parser
        self.num_every_parser = num_every_parser
        
    def allocate_if_necessary(self):
        return self.num_steps_parser.result() / self.num_every_parser.result()
        
        
class CoordinatesGradientsParser(MultiLineParser, TrajectoryParserMixin):
    extension = ".out"
    
    def __init__(self, num_steps_parser, num_every_parser, elements_parser, label='coor_grad', condition=None):
        TrajectoryParserMixin.__init__(self, num_steps_parser, num_every_parser)
        MultiLineParser.__init__(
            self,
            label,
            re.compile(r"ATOM\s+COORDINATES\s+GRADIENTS"),
            re.compile(r"^$"),
            condition
        )
        self.elements_parser = elements_parser
        self.re = re.compile(r"\d+\s+\S+\s+(?P<x>\S+)\s+(?P<y>\S+)\s+(?P<z>\S+)\s+(?P<gx>\S+)\s+(?P<gy>\S+)\s+(?P<gz>\S+)")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.coordinates = None
        self.gradients = None
        self.step_counter = None
        
    def allocate_if_necessary(self):
        if self.step_counter == None:
            self.num_steps = TrajectoryParserMixin.allocate_if_necessary(self) + 1
            num_atoms = len(self.elements_parser.result())
            self.coordinates = Numeric.zeros((self.num_steps, num_atoms, 3), Numeric.Float)
            self.gradients = Numeric.zeros((self.num_steps, num_atoms, 3), Numeric.Float)
            self.step_counter = 0
        
    def start_collecting(self):
        self.current_gradient = []
        self.current_coordinates = []
        self.allocate_if_necessary()
        
    def collect(self, line):
        match = self.re.search(line)
        if match != None:
            self.current_coordinates.append([
                float(match.group("x")),
                float(match.group("y")),
                float(match.group("z"))
            ])
            self.current_gradient.append([
                float(match.group("gx")),
                float(match.group("gy")),
                float(match.group("gz"))
            ])

    def stop_collecting(self):
        self.coordinates[self.step_counter,:,:] = Numeric.array(self.current_coordinates)
        self.gradients[self.step_counter,:,:] = Numeric.array(self.current_gradient)
        self.step_counter += 1

    def result(self):
        #assert self.step_counter == self.num_steps
        return self.coordinates, self.gradients
        

class EnergiesParser(MultiLineParser, TrajectoryParserMixin):
    extension = ".out"

    def __init__(self, num_steps_parser, num_every_parser, energy_name='TOTAL ENERGY', label='energies', condition=None):
        FileParser.__init__(self, label, condition)
        TrajectoryParserMixin.__init__(self, num_steps_parser, num_every_parser)
        self.re = re.compile(r"%s =\s+(?P<energy>\S+)\s+A.U." % energy_name)
        
    def reset(self):
        MultiLineParser.reset(self)
        self.energies = None
        self.step_counter = None
        
    def allocate_if_necessary(self):
        if self.step_counter == None:
            self.num_steps = TrajectoryParserMixin.allocate_if_necessary(self)
            self.energies = Numeric.zeros(self.num_steps, Numeric.Float)
            self.step_counter = 0
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.allocate_if_necessary()
            self.energies[self.step_counter] = float(match.group("energy"))
        
    def result(self):
        #assert self.step_counter == self.num_steps
        return self.energies
