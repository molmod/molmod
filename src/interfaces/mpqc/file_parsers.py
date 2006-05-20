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

import re, numpy

class ScfEnergiesParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='scf_energies', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"total scf energy\s+=\s+(?P<energy>\S+)")
        
    def reset(self):
        self.energies = []
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.energies.append(float(match.group("energy")))
        
    def result(self):
        return self.energies


class MolecularEnergiesParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='molecular_energies', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"Value of the MolecularEnergy:\s+(?P<energy>\S+)")
        
    def reset(self):
        self.energies = []
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.energies.append(float(match.group("energy")))
        
    def result(self):
        return self.energies


class EnergyAccuracyParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='energy_accuracy', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"value_accuracy\s+=\s+(?P<energy_accuracy>\S+)")
        
    def reset(self):
        self.energy_accuracy = None
        
    def parse(self, line):
        if self.energy_accuracy == None:
            match = self.re.search(line)
            if match != None:
                self.energy_accuracy = float(match.group("energy_accuracy"))
        
    def result(self):
        return self.energy_accuracy


class WarningParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='warnings', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"WARNING:")
        
    def reset(self):
        self.warnings = False
        
    def parse(self, line):
        if not self.warnings:
            match = self.re.search(line)
            if match != None:
                self.warnings = True
        
    def result(self):
        return self.warnings


class OptimizationConvergedParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='optimization_converged', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"The optimization has converged.")
        
    def reset(self):
        self.converged = False
        
    def parse(self, line):
        if not self.converged:
            match = self.re.search(line)
            if match != None:
                self.converged = True
        
    def result(self):
        return self.converged


class OutputMoleculesParser(MultiLineParser):
    filename = ".out"
    extension = True

    def __init__(self, label='output_molecules', condition=None):
        activator = re.compile(r"n\s+atoms\s+geometry")
        deactivator = re.compile(r"}$")
        MultiLineParser.__init__(self, label, activator, deactivator, condition)
        self.re = re.compile(r"(?P<symbol>\S+)\s*\[\s*(?P<x>\S+)\s*(?P<y>\S+)\s*(?P<z>\S+)\s*\]")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.molecules = []
        
    def start_collecting(self):
        self.current_atoms = []

    def collect(self, line):
        match = self.re.search(line)
        self.current_atoms.append([
            periodic.symbol_lookup(match.group("symbol")),
            float(match.group("x")),
            float(match.group("y")),
            float(match.group("z"))
        ])

    def stop_collecting(self):
        self.molecules.append(Molecule(self.current_atoms))
        del self.current_atoms
        
    def result(self):
        return self.molecules


class GradientsParser(MultiLineParser):
    filename = ".out"
    extension = True

    def __init__(self, label='gradients', condition=None):
        activator = re.compile(r"Total Gradient")
        deactivator = re.compile(r"^$")
        MultiLineParser.__init__(self, label, activator, deactivator, condition)
        self.re = re.compile(r"\d+\s+\S+\s+(?P<gradient_x>\S+)\s+(?P<gradient_y>\S+)\s+(?P<gradient_z>\S+)")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.gradients = []
        
    def start_collecting(self):
        self.current_gradient = []

    def collect(self, line):
        match = self.re.search(line)
        if match != None:
            self.current_gradient.append([
                float(match.group("gradient_x")),
                float(match.group("gradient_y")),
                float(match.group("gradient_z"))
            ])
        
    def stop_collecting(self):
        gradient = numpy.array(self.current_gradient, float)
        self.gradients.append(gradient)
        del self.current_gradient

    def result(self):
        return self.gradients


class GradientAccuracyParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='gradient_accuracy', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"gradient_accuracy\s+=\s+(?P<gradient_accuracy>\S+)")
        
    def reset(self):
        self.gradient_accuracy = None
        
    def parse(self, line):
        if self.gradient_accuracy == None:
            match = self.re.search(line)
            if match != None:
                self.gradient_accuracy = float(match.group("gradient_accuracy"))
        
    def result(self):
        return self.gradient_accuracy


class HessianParser(FileParser):
    filename = ".hess"
    extension = True

    def __init__(self, label='hessian', condition=None):
        FileParser.__init__(self, label, condition)
        self.re_num_atoms = re.compile(r"(?P<num_atoms>\d+)\s+atoms")
        
    def reset(self):
        self.num_atoms = None
        self.begin_line = None
        self.end_line = None
        self.current_line = 0
        self.hessian_elements = []
        
    def parse(self, line):
        if self.num_atoms == None:
            match = self.re_num_atoms.search(line)
            if match != None:
                self.num_atoms = int(match.group("num_atoms"))
                num_elements = self.num_atoms*3 * (self.num_atoms*3 + 1) / 2
                num_lines = num_elements / 5
                if num_elements % 5 > 0: num_lines += 1
                self.begin_line = self.num_atoms + 2
                self.end_line = self.begin_line + num_lines
        elif (self.current_line >= self.begin_line) and (self.current_line < self.end_line):
            self.hessian_elements.extend(float(word) for word in line.split())
            #print line
            #print [float(word) for word in line.split()]
        self.current_line += 1

    def result(self):
        if self.num_atoms != None:
            result = numpy.zeros((self.num_atoms*3, self.num_atoms*3), float)
            counter = 0
            for i in xrange(self.num_atoms*3):
                for j in xrange(0, i):
                    result[i,j] = self.hessian_elements[counter]
                    result[j,i] = self.hessian_elements[counter]
                    counter += 1
                result[i,i] = self.hessian_elements[counter]
                counter += 1
            return result


class RawGridParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='grid', trigger=None, condition=None):
        self.trigger = trigger
        FileParser.__init__(self, label, condition)

    def reset(self):
        self.grid = None
        self.active = (self.trigger == None)
        self.read_nrecords = (self.trigger == None)
        
    def parse(self, line):
        if self.active:
            words = line.split()
            if len(words) == 4:
                try:
                    self.grid[self.counter, 0] = float(words[0])
                    self.grid[self.counter, 1] = float(words[1])
                    self.grid[self.counter, 2] = float(words[2])
                    self.grid[self.counter, 3] = float(words[3])
                except ValueError:
                    self.active = False
            else:
                self.active = False
            self.counter += 1
            if self.counter >= self.nrecords:
                self.active = False
        elif self.read_nrecords:
            assert line.startswith("# Number of records:")
            words = line.split()
            self.nrecords = int(words[-1])
            self.grid = numpy.zeros((self.nrecords, 4), float)
            self.counter = 0
            self.active = True
            self.read_nrecords = False
        elif self.grid == None:
            if line == self.trigger: self.read_nrecords = True

    def result(self):
        return self.grid


class TotalChargeParser(FileParser):
    filename = ".out"
    extension = True

    def __init__(self, label='total_charge', condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"total charge =\s+(?P<total_charge>\S+)")
        
    def reset(self):
        self.total_charge = None
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.total_charge = float(match.group("total_charge"))
        
    def result(self):
        return self.total_charge


class NPAParser(MultiLineParser):
    filename = ".out"
    extension = True

    def __init__(self, label='npa_charges', condition=None):
        activator = re.compile(r"Natural Population Analysis:")
        deactivator = re.compile(r"^$")
        MultiLineParser.__init__(self, label, activator, deactivator, condition)
        self.re = re.compile(r"^\s+\d+\s+\S+\s+(?P<npa_charge>\S+)")
        
    def start_collecting(self):
        self.npa_charges = []

    def collect(self, line):
        match = self.re.search(line)
        if match != None:
            self.npa_charges.append(float(match.group("npa_charge")))

    def stop_collecting(self):
        pass

    def result(self):
        return numpy.array(self.npa_charges)
