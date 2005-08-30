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

from pychem.moldata import periodic
import os, copy, glob
import Numeric
from pychem import context

from pychem.units import to_angstrom, from_angstrom


__all__ = [
    "ExternalError", "MpqcJob", "SimpleMpqcJob", 
    "SimpleMpqcJobSinglePoint", "SimpleMpqcJobOptimize"
]


class ExternalError(Exception):
    """
    This error is raised when an external executable doesn't appear to do it's
    job. This is judged on the completeness of the generated output files.
    """
    pass


class MpqcJob(object):
    """
    This is the base class for all MPQC calculations.
    
    Arguments:
    filename -- This is the base name for input (.in), ouput (.out), 
                summary (.smr), ... filenames.
    title -- The description of the calculation in the input file.
    input_molecule -- 
    """
    
    def __init__(self, filename, title, input_molecule):
        self.filename = filename
        self.title = title
        self.input_molecule = input_molecule

    def write_input(self, f):
        raise NotImplementedError

    def run_external(self, overwrite=False):
        """
        Call the external program.

        Note that the calculation is only executed if no output file is found,
        unless overwrite == True
        """
        if (not os.path.isfile("%s.out" % self.filename)) or overwrite:
            os.system("mpqc -o %s.out %s.in" % (self.filename, self.filename))
            for temp_filename in glob.glob("%s.wfn.*.tmp" % self.filename):
                os.remove(temp_filename)
            return False
        else:
            return True

    def awk_scriptname(self):
        """Translate the class name into a base name for the awk script filename."""
        result = ""
        class_name = str(self.__class__)
        class_name = class_name[class_name.rfind(".")+1:-2]
        for char in str(class_name):
            if char == char.lower():
                result += char
            elif len(result)==0:
                result += char.lower()
            else:
                result += "_"+char.lower()
        return result

    def summarize_output(self):
        """
        Generate a summary file based on the output of the calculation.
        Return the sumary file interpreted as a python expression.
        """
        os.system(
            "gawk -f %sinterfaces/awk/%s.awk < %s.out > %s.smr" % (
                context.share_path, self.awk_scriptname(),
                self.filename, self.filename
            )
        )
        smr = file("%s.smr" % self.filename)
        result = eval(''.join(smr))
        smr.close()
        return result

    def process_summary(self, f):
        """Process the attributes taken from the summary file and assigned to self."""
        raise NotImplementedError        
        
    def run(self):
        """Perform the complete calculation and analysis."""
        print "running job: %s" % self.filename
        self.write_input(file(self.filename + ".in", 'w'))
        recycled = self.run_external()
        print "output recycled: %s" % recycled
        self.__dict__.update(self.summarize_output())
        print "job completed: %s" % self.job_completed
        if recycled and not self.job_completed:
            print "trying again: %s"
            recycled = self.run_external(overwrite=True)
        print "job completed: %s" % self.job_completed
        if not self.job_completed:
            raise ExternalError("Output file of external job is not complete (%s)" % self.filename)
        self.process_summary()
        return recycled


class SimpleMpqcJob(MpqcJob):
    """MPQC jobs that use the simple input format."""
    
    def __init__(self, filename, title, input_molecule, method, basis):
        """
        Initialize a SimpleMpqcJob instance.
        
        New arguments:
        method -- The type of approximation for the many body electron
                  Hamiltonian
        basis -- The basis set used to describe the wave-function.
        """
        MpqcJob.__init__(self, filename, title, input_molecule)
        self.method = method
        self.basis = basis
        
    def write_input(self, f):
        print >> f, "% " + self.title
        print >> f, "method: " + self.method
        print >> f, "basis: " + self.basis
        print >> f, "charge: " + str(self.input_molecule.charge)
        print >> f, "multiplicity: " + str(self.input_molecule.spin_multiplicity)
        print >> f, "molecule: "
        for number, (x, y, z) in zip(self.input_molecule.numbers, to_angstrom(self.input_molecule.coordinates)):
            print >> f, "   %2s  % 10.7f  % 10.7f  % 10.7f" % (periodic.symbol[number], x, y, z)


def yesno(value):
    if value:
        return "yes"
    else:
        return "no"

class SimpleMpqcJobSinglePoint(SimpleMpqcJob):
    """
    Simple MPQC jobs that doesn't change the geometry of the molecule.
    Only one SCF calculation is performed, with eventual post analysis.
    """
    def __init__(self, filename, title, input_molecule, method, basis, do_gradient=False):
        """
        Initialize a SimpleMpqcJobSinglePoint instance.
        
        Extra arguments:
        do_gradient -- wether calculate the gradient of the energy
        """
        SimpleMpqcJob.__init__(self, filename, title, input_molecule, method, basis)
        self.do_gradient = do_gradient
        self.energy = None
        self.gradient = None
        
    def write_input(self, f):
        SimpleMpqcJob.write_input(self, f)
        print >> f, "optimize: no"
        print >> f, "gradient: " + yesno(self.do_gradient)
        
    def process_summary(self):
        if self.gradient != None:
            self.gradient = Numeric.array(self.gradient, Numeric.Float)
        self.input_molecule.coordinates = from_angstrom(Numeric.array(self.input_coordinates, Numeric.Float))
        del self.input_coordinates


class SimpleMpqcJobOptimize(SimpleMpqcJob):
    """
    A Simple MPQC job that optimizes the geometry of the molecule towards
    lower energies. The default MPQC optimization scheme is used.
    """
    def __init__(self, filename, title, input_molecule, method, basis):
        """Initialize a SimpleMpqcJobOptimize instance."""
        SimpleMpqcJob.__init__(self, filename, title, input_molecule, method, basis)
        self.energies = []
        self.output_molecule = None
        self.gradient = None
        
    def write_input(self, f):
        SimpleMpqcJob.write_input(self, f)
        print >> f, "optimize: yes"
        
    def process_summary(self):
        self.output_molecule = copy.deepcopy(self.input_molecule)
        self.output_molecule.coordinates = from_angstrom(Numeric.array(self.output_coordinates, Numeric.Float))
        del self.output_coordinates
