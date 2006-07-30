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


from molmod.interfaces.output_parsers import FileParser, MultiLineParser
from molmod.data import periodic
from molmod.molecules import Molecule
from molmod.units import from_unified

import re, numpy

class NumStepsParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='num_steps', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"MAXIMUM NUMBER OF STEPS:\s+(?P<num_steps>\S+)\s+STEPS")
        
    def reset(self):
        self.num_steps = None
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.num_steps = int(match.group("num_steps"))
        
    def result(self):
        return self.num_steps


class NumEveryParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='num_every', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"PRINT INTERMEDIATE RESULTS EVERY\s+(?P<num_every>\S+)\s+STEPS")
        
    def reset(self):
        self.num_every = None
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.num_every = int(match.group("num_every"))
        
    def result(self):
        return self.num_every


class TimeStepsParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='time_steps', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"TIME STEP FOR (?P<component>\S+):\s+(?P<time_step>\S+)")
        
    def reset(self):
        self.time_step_electrons = None
        self.time_step_ions = None
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            component = match.group("component")
            time_step = float(match.group("time_step"))
            if component == "ELECTRONS":
                self.time_step_electrons = time_step
            elif component == "IONS":
                self.time_step_ions = time_step
        
    def result(self):
        return self.time_step_electrons, self.time_step_ions


class ElementsParser(MultiLineParser):
    filename = ".out"
    extension = True

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
                self.elements.append(periodic[match.group("symbol")].number)
        
    def stop_collecting(self):
        self.done = True

    def result(self):
        return self.elements


class MassesParser(MultiLineParser):
    filename = ".out"
    extension = True

    def __init__(self, label='masses', condition=None):
        MultiLineParser.__init__(
            self,
            label,
            re.compile(r"ATOM\s+MASS\s+RAGGIO\s+NLCC\s+PSEUDOPOTENTIAL"),
            re.compile(r"^$"),
            condition
        )
        self.re = re.compile(r"\*\s+(?P<symbol>\S+)\s+(?P<mass>\S+)\s+\S+\s+\S+\s+\S+\s+\*")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.done = False
        
    def start_collecting(self):
        self.masses = {}

    def collect(self, line):
        if not self.done:
            match = self.re.search(line)
            if match != None:
                number = periodic[match.group("symbol")].number
                mass = from_unified(float(match.group("mass")))
                self.masses[number] = mass
        
    def stop_collecting(self):
        self.done = True

    def result(self):
        return self.masses


class EveryParserMixin(object):
    filename = ".out"
    extension = True
    
    def __init__(self, num_steps_parser, num_every_parser):
        self.num_steps_parser = num_steps_parser
        self.num_every_parser = num_every_parser
        
    def allocate_if_necessary(self):
        return self.num_steps_parser.result() / self.num_every_parser.result()


class CoordinatesGradientsParser(MultiLineParser, EveryParserMixin):
    filename = ".out"
    extension = True
    
    def __init__(self, num_steps_parser, num_every_parser, elements_parser, label='coor_grad', condition=None):
        EveryParserMixin.__init__(self, num_steps_parser, num_every_parser)
        MultiLineParser.__init__(
            self,
            label,
            re.compile(r"ATOM\s+COORDINATES\s+GRADIENTS"),
            re.compile(r"^$"),
            condition,
            [num_steps_parser, num_every_parser]
        )
        self.elements_parser = elements_parser
        self.re = re.compile(r"\d+\s+\S+(?P<x>.{8})(?P<y>.{8})(?P<z>.{8}) (?P<gx>.{11})(?P<gy>.{11})(?P<gz>.{11})")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.coordinates = None
        self.gradients = None
        self.step_counter = None
        
    def allocate_if_necessary(self):
        if self.step_counter == None:
            self.num_steps = EveryParserMixin.allocate_if_necessary(self) + 1
            num_atoms = len(self.elements_parser.result())
            self.coordinates = numpy.zeros((self.num_steps, num_atoms, 3), float)
            self.gradients = numpy.zeros((self.num_steps, num_atoms, 3), float)
            self.step_counter = 0
        
    def start_collecting(self):
        self.allocate_if_necessary()
        if self.step_counter < self.num_steps:
            self.current_gradient = []
            self.current_coordinates = []
        
    def collect(self, line):
        if self.step_counter < self.num_steps:
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
        if self.step_counter < self.num_steps:
            self.coordinates[self.step_counter,:,:] = numpy.array(self.current_coordinates)
            self.gradients[self.step_counter,:,:] = numpy.array(self.current_gradient)
            self.step_counter += 1

    def result(self):
        #assert self.step_counter == self.num_steps
        return self.coordinates, self.gradients
        

class EnergiesParser(MultiLineParser, EveryParserMixin):
    filename = ".out"
    extension = True

    def __init__(self, num_steps_parser, num_every_parser, energy_name='TOTAL ENERGY', label='energies', condition=None):
        FileParser.__init__(self, label, condition, [num_steps_parser, num_every_parser])
        EveryParserMixin.__init__(self, num_steps_parser, num_every_parser)
        self.re = re.compile(r"%s =\s+(?P<energy>\S+)\s+A.U." % energy_name)
        
    def reset(self):
        MultiLineParser.reset(self)
        self.energies = None
        self.step_counter = None
        
    def allocate_if_necessary(self):
        if self.step_counter == None:
            self.num_steps = EveryParserMixin.allocate_if_necessary(self)
            self.energies = numpy.zeros(self.num_steps, float)
            self.step_counter = 0
        
    def parse(self, line):
        if (self.step_counter == None) or (self.step_counter < self.num_steps):
            match = self.re.search(line)
            if match != None:
                self.allocate_if_necessary()
                self.energies[self.step_counter] = float(match.group("energy"))
                self.step_counter += 1
        
    def result(self):
        #assert self.step_counter == self.num_steps
        return self.energies


class EnergiesFileParser(FileParser):
    filename = "ENERGIES"
    extionsion = False
    
    def __init__(self, num_steps_parser, label='energies_table', condition=None):
        FileParser.__init__(self, label, condition, [num_steps_parser])
        self.num_steps_parser = num_steps_parser
        
    def reset(self):
        self.num_steps = self.num_steps_parser.result()
        self.energies = numpy.zeros((self.num_steps, 5), float)
        self.counter = 0
        
    def parse(self, line):
        if self.counter < self.num_steps:
            self.energies[self.counter] = [float(word) for word in line.split()[1:6]]
            self.counter += 1
        
    def result(self):
        return self.energies


class TrajectoryFileParser(FileParser):
    filename = "TRAJECTORY"
    extionsion = False
    
    def __init__(self, num_steps_parser, elements_parser, label='coor_velo', condition=None):
        FileParser.__init__(self, label, condition, [num_steps_parser, elements_parser])
        self.num_steps_parser = num_steps_parser
        self.elements_parser = elements_parser
        
    def reset(self):
        self.num_atoms = len(self.elements_parser.result())
        self.num_steps = self.num_steps_parser.result()
        self.coordinates = numpy.zeros((self.num_steps, self.num_atoms, 3), float)
        self.velocities = numpy.zeros((self.num_steps, self.num_atoms, 3), float)
        self.counter = 0
        self.atom_counter = 0
        
    def parse(self, line):
        if self.counter < self.num_steps:
            data = [float(word) for word in line.split()[1:7]]
            self.coordinates[self.counter, self.atom_counter] = data[0:3]
            self.velocities[self.counter, self.atom_counter] = data[3:6]
            self.atom_counter += 1
            if self.atom_counter >= self.num_atoms:
                self.atom_counter = 0
                self.counter += 1
        
    def result(self):
        return self.coordinates, self.velocities


class CellDimensionParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='cell_dimension', condition=None):
        FileParser.__init__(self, label, condition)
        
    def reset(self):
        self.cell_dimension = None
        
    def parse(self, line):
        if line[0:16] == " CELL DIMENSION:":
            words = line[16:].split()
            self.cell_dimension = [float(word) for word in words]
        
    def result(self):
        return self.cell_dimension


