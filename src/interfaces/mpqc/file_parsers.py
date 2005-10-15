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

from pychem.interfaces.output_parsers import FileParser
from pychem.moldata import periodic
from pychem.molecules import Molecule

import re, Numeric

class ScfEnergiesParser(FileParser):
    extension="out"

    def __init__(self, label, condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"total scf energy =\s(?P<energy>\S+)")
        
    def reset(self):
        self.energies = []
        
    def parse(self, line):
        match = self.re.search(line)
        if match != None:
            self.energies.append(float(match.group("energy")))
        
    def result(self):
        return self.energies


class MolecularEnergiesParser(FileParser):
    extension="out"

    def __init__(self, label, condition=None):
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


class WarningParser(FileParser):
    extension="out"

    def __init__(self, label, condition=None):
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
    extension="out"

    def __init__(self, label, condition=None):
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


class MultiLineParser(FileParser):
    def __init__(self, label, activator, deactivator, condition=None):
        FileParser.__init__(self, label, condition)
        self.activator = activator
        self.deactivator = deactivator

    def reset(self):
        self.active = False

    def parse(self, line):
        if self.active:
            if self.deactivator.search(line) != None:
                self.active = False
                self.stop_collecting()
            else:
                self.collect(line)
        elif self.activator.search(line) != None:
            self.active = True
            self.start_collecting()

    def start_collecting(self):
        raise NotImplementedError

    def collect(self, line):
        raise NotImplementedError

    def stop_collecting(self):
        raise NotImplementedError


class OutputMoleculesParser(MultiLineParser):
    extension="out"

    def __init__(self, label, condition=None):
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
            periodic.reverse_symbol_lookup(match.group("symbol")),
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
    extension="out"

    def __init__(self, label, condition=None):
        activator = re.compile(r"Gradient of the MolecularEnergy:")
        deactivator = re.compile(r"^$")
        MultiLineParser.__init__(self, label, activator, deactivator, condition)
        self.re = re.compile(r"\d+\s+(?P<gradient>\S+)")
        
    def reset(self):
        MultiLineParser.reset(self)
        self.gradients = []
        
    def start_collecting(self):
        self.current_gradient = []

    def collect(self, line):
        match = self.re.search(line)
        self.current_gradient.append(float(match.group("gradient")))
        
    def stop_collecting(self):
        gradient = Numeric.array(self.current_gradient, Numeric.Float)
        gradient.shape = (-1,3)
        self.gradients.append(gradient)
        del self.current_gradient

    def result(self):
        return self.gradients


class GradientAccuracyParser(MultiLineParser):
    extension="out"

    def __init__(self, label, condition=None):
        FileParser.__init__(self, label, condition)
        self.re = re.compile(r"gradient accuracy = (?P<gradient_accuracy>\S+)")
        
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
    extension="hess"

    def __init__(self, label, condition=None):
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
            result = Numeric.zeros((self.num_atoms*3, self.num_atoms*3), Numeric.Float)
            counter = 0
            for i in xrange(self.num_atoms*3):
                for j in xrange(0, i):
                    result[i,j] = self.hessian_elements[counter]
                    result[j,i] = self.hessian_elements[counter]
                    counter += 1
                result[i,i] = self.hessian_elements[counter]
                counter += 1
            return result
