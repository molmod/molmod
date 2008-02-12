# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


from molmod.io.output_parsers import FileParser, MultiLineParser

from molmod.molecules import Molecule
from molmod.units import angstrom, unified

import numpy

import re


__all__ = [
    "LinkParser", "ThermoChemParser", "HessianParser", "FrequenciesParser",
    "LowFrequenciesParser", "SelectedFrequenciesParser", "MassParser",
    "GradientParser", "InputOrientationGradientParser",
    "StandardOrientationGradientParser", "ConfigurationParser",
    "CoordinatesParser", "StandardOrientationCoordinatesParser",
    "InputOrientationCoordinatesParser", "OptimizedParser", "IsOptimizedParser",
    "OptimizedCoordinatesParser", "SCFParser", "EnergyParser"
]


class LinkParser(MultiLineParser):
    filename = ".log"
    extension = True

    def __init__(self, link, label, activator=None, deactivator=None, condition=None, depends_on=[]):
        MultiLineParser.__init__(self, label, activator, deactivator, condition, depends_on)
        self.link = str(link)

    def reset(self):
        MultiLineParser.reset(self)
        self.in_link = False

    def parse(self, line):
        if line[:11] == " Leave Link" and line[13:13+len(self.link)] == self.link:
            self.in_link = False
        if self.in_link:
            MultiLineParser.parse(self, line)
        if line[:8] == " (Enter " and line[-6-len(self.link):-6] == self.link:
            self.in_link = True


class ThermoChemParser(LinkParser):
    def __init__(self, label, activator=None, deactivator=None, condition=None, depends_on=[]):
        LinkParser.__init__(self, "716", label, activator, deactivator, condition, depends_on)


class HessianParser(ThermoChemParser):
    def __init__(self, label="hessian", condition=None):
        ThermoChemParser.__init__(self, label,
            activator=re.compile("Force constants in Cartesian coordinates:"),
            deactivator=re.compile(r"^\s*\b[^0-9+\-]"),
            condition=condition
        )

    def reset(self):
        ThermoChemParser.reset(self)
        self.hessian = []

    def start_collecting(self):
        self.hessian = []

    def collect(self, line):
        words = line.split()
        if (len(words) > 1) and (words[1].find("D") >= 0):
            row_number = int(words[0])
            if row_number > len(self.hessian):
                row = []
                self.hessian.append(row)
            else:
                row = self.hessian[row_number-1]
            for word in words[1:]:
                row.append(float(word.replace("D", "e")))

    def stop_collecting(self):
        hessian = numpy.zeros((len(self.hessian), len(self.hessian)), float)
        for row_index, row in enumerate(self.hessian):
            for col_index, value in enumerate(row):
                hessian[row_index, col_index] = value
                if row_index != col_index:
                    hessian[col_index, row_index] = value
        self.hessian = hessian

    def result(self):
        return self.hessian


class FrequenciesParser(ThermoChemParser):
    def __init__(self, label, pattern, condition):
        # returns the frequencies in cm-1
        ThermoChemParser.__init__(self, label, None, None, condition)
        self.pattern = pattern

    def reset(self):
        ThermoChemParser.reset(self)
        self.frequencies = []

    def collect(self, line):
        if line[:len(self.pattern)] == self.pattern:
            words = line[len(self.pattern):].split()
            self.frequencies.extend(float(word) for word in words)

    def result(self):
        return numpy.array(self.frequencies)


class LowFrequenciesParser(FrequenciesParser):
    def __init__(self, label="low_frequencies", condition=None):
        FrequenciesParser.__init__(self, label, " Low frequencies ---", condition)


class SelectedFrequenciesParser(FrequenciesParser):
    def __init__(self, label="selected_frequencies", condition=None):
        FrequenciesParser.__init__(self, label, " Frequencies --", condition)


class MassParser(ThermoChemParser):
    def __init__(self, label="masses", condition=None):
        ThermoChemParser.__init__(self, label,
            activator=re.compile("Temperature\s+\S+\s+Kelvin.\s+Pressure\s+\S+\s+Atm."),
            deactivator=re.compile("Molecular mass:\s+\S+\s+amu."),
            condition=condition
        )
        self.re = re.compile("Atom\s*\d+\s+has atomic number\s+\d+\s+and mass\s+(?P<mass>\S+)")

    def reset(self):
        ThermoChemParser.reset(self)
        self.masses = []

    def start_collecting(self):
        self.masses = []

    def collect(self, line):
        match = self.re.search(line)
        if match is not None:
            self.masses.append(float(match.group("mass"))*unified)

    def stop_collecting(self):
        self.masses = numpy.array(self.masses, float)

    def result(self):
        return self.masses


class GradientParser(ThermoChemParser):
    def __init__(self, label, activator, deactivator, condition=None):
        ThermoChemParser.__init__(self, label, activator, deactivator, condition)
        self.re = re.compile("\d+\s+\d+\s+(?P<fx>\S+)\s+(?P<fy>\S+)\s+(?P<fz>\S+)")

    def reset(self):
        ThermoChemParser.reset(self)
        self.gradient_list = []

    def start_collecting(self):
        self.gradient = []

    def collect(self, line):
        match = self.re.search(line)
        if match is not None:
            self.gradient.append([
                -float(match.group("fx")),
                -float(match.group("fy")),
                -float(match.group("fz"))
            ])

    def stop_collecting(self):
        self.gradient_list.append(numpy.array(self.gradient, float))

    def result(self):
        return self.gradient_list


class InputOrientationGradientParser(GradientParser):
    def __init__(self, label="io_gradient_list", condition=None):
        GradientParser.__init__(self, label,
            activator=re.compile("\*\*\*\*\* Axes restored to original set \*\*\*\*\*"),
            deactivator=re.compile("Cartesian Forces:"),
            condition=condition
        )


class StandardOrientationGradientParser(GradientParser):
    def __init__(self, label="so_gradient_list", condition=None):
        GradientParser.__init__(self, label,
            activator=re.compile("Forces in standard orientation"),
            deactivator=re.compile("\*\*\*\*\* Axes restored to original set \*\*\*\*\*"),
            condition=condition
        )


class ConfigurationParser(LinkParser):
    def __init__(self, label, activator=None, deactivator=None, condition=None, depends_on=[]):
        LinkParser.__init__(self, "202", label, activator, deactivator, condition, depends_on)


class CoordinatesParser(ConfigurationParser):
    def __init__(self, label, activator, deactivator, condition=None):
        ConfigurationParser.__init__(self, label, activator, deactivator, condition)
        self.re = re.compile("\d+\s+\d+\s+\d+\s+(?P<x>\S+)\s+(?P<y>\S+)\s+(?P<z>\S+)")

    def reset(self):
        ConfigurationParser.reset(self)
        self.coordinates = []

    def start_collecting(self):
        self.current_coordinates = []

    def collect(self, line):
        match = self.re.search(line)
        if match is not None:
            self.current_coordinates.append([
                float(match.group("x")) * angstrom,
                float(match.group("y")) * angstrom,
                float(match.group("z")) * angstrom,
            ])

    def stop_collecting(self):
        self.coordinates.append(numpy.array(self.current_coordinates, float))

    def result(self):
        return self.coordinates


class StandardOrientationCoordinatesParser(CoordinatesParser):
    def __init__(self, label="so_coordinates_list", condition=None):
        CoordinatesParser.__init__(self, label,
            re.compile("Standard orientation"),
            re.compile("Rotational constants"),
            condition
        )


class InputOrientationCoordinatesParser(CoordinatesParser):
    def __init__(self, label="io_coordinates_list", condition=None):
        CoordinatesParser.__init__(self, label,
            re.compile("Input orientation"),
            re.compile("Standard orientation"),
            condition
        )


class OptimizedParser(LinkParser):
    def __init__(self, label, activator=None, deactivator=None, condition=None, depends_on=[]):
        LinkParser.__init__(self, "103", label, activator, deactivator, condition, depends_on)


class IsOptimizedParser(OptimizedParser):
    def __init__(self, label="optimized", condition=None):
        OptimizedParser.__init__(self, label, None, None, condition)
        self.re = re.compile("-- Stationary point found\.")

    def reset(self):
        OptimizedParser.reset(self)
        self.optimized = False

    def collect(self, line):
        if not self.optimized and self.re.search(line) is not None:
            self.optimized = True

    def result(self):
        return self.optimized


class OptimizedCoordinatesParser(OptimizedParser):
    def __init__(self, label="optimized_coordinates", condition=None):
        OptimizedParser.__init__(self, label,
            re.compile("Optimized Parameters"),
            re.compile("GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad"),
            condition
        )
        self.re = re.compile("\S+\s+R\(\d+,-\d\)\s+(?P<coordinate>\S+)\s+-DE/DX")

    def reset(self):
        OptimizedParser.reset(self)
        self.completed = False

    def start_collecting(self):
        if not self.completed:
            self.optimized_coordinates = []

    def collect(self, line):
        if not self.completed:
            match = self.re.search(line)
            if match is not None:
                self.optimized_coordinates.append(float(match.group("coordinate"))*angstrom)

    def stop_collecting(self):
        if not self.completed:
            self.optimized_coordinates = numpy.array(self.optimized_coordinates, float)
            self.optimized_coordinates.shape = (-1, 3)
            self.completed = True

    def result(self):
        return self.optimized_coordinates


class SCFParser(LinkParser):
    def __init__(self, label, activator=None, deactivator=None, condition=None, depends_on=[]):
        LinkParser.__init__(self, "502", label, activator, deactivator, condition, depends_on)


class EnergyParser(SCFParser):
    def __init__(self, label="energies", condition=None):
        SCFParser.__init__(self, label, None, None, condition)
        self.re = re.compile("SCF Done:\s+E\S+\s+=\s+(?P<energy>\S+)\s+A.U.")

    def reset(self):
        SCFParser.reset(self)
        self.energies = []

    def collect(self, line):
        match = self.re.search(line)
        if match is not None:
            self.energies.append(float(match.group("energy")))

    def result(self):
        return numpy.array(self.energies)


class PopulationParser(LinkParser):
    def __init__(self, label, activator=None, deactivator=None, condition=None, depends_on=[]):
        LinkParser.__init__(self, "602", label, activator, deactivator, condition, depends_on)


class ESPFitParser(PopulationParser):
    def __init__(self, label="espfit", condition=None):
        PopulationParser.__init__(self, label,
            re.compile("Charges from ESP fit"),
            re.compile("-----------------------------------------------------------------"),
            condition
        )

    def start_collecting(self):
        self.espfit = []
        self.started = False

    def collect(self, line):
        if self.started:
            self.espfit.append([float(word) for word in line.split()[2:]])
        if line.startswith("              1"):
            self.started = True

    def stop_collecting(self):
        pass

    def result(self):
        return numpy.array(self.espfit)





