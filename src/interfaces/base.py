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


from pychem import context


from StringIO import StringIO
from pickle import load, dump
import sha, os, glob


class ExternalError(Exception):
    """
    This error is raised when an external executable doesn't appear to do it's
    job. This is judged on the completeness of the generated output files.
    """
    pass


class Job(object):
    def __init__(self, prefix, title):
        self.prefix = prefix
        self.title = title
        self.ran = False
        input_string_io = StringIO()
        self.write_input(input_string_io)
        self.input_string = input_string_io.getvalue()
        input_string_io.close()
        self.filename = "%s_%s" % (self.prefix, sha.new(self.input_string).hexdigest())
        
    def write_input(self, f):
        raise NotImplementedError
    
    def create_job_file(self):
        f = file(self.filename + ".job", 'w')
        dump(self, f)
        f.close()
    
    def create_input_file(self):
        input_file = file(self.filename + ".in", 'w')
        input_file.write(self.input_string)
        input_file.close()
        
    def cleanup(self):
        for temp_filename in glob.glob("%s*.*" % self.filename):
            os.remove(temp_filename)        

    def output_file_exists(self):
        return os.path.isfile(self.filename + ".out")

    def external_command(self):
        raise NotImplementedError
        
    def remove_temporary_files(self):
        raise NotImplementedError
        
    def determine_completed(self):
        raise NotImplementedError

    def run_external(self, overwrite=False):
        """
        Call the external program.

        Note that the calculation is only executed if no output file is found,
        unless overwrite == True
        """
        recycled = False
        if (not self.output_file_exists()) or overwrite:
            os.system(self.external_command())
            self.remove_temporary_files()
        else:
            recycled = True
        self.determine_completed()
        return recycled
            
    def read_output(self):
        raise NotImplementedError
            
    def run(self, user_overwrite=False):
        """Perform the complete calculation and analysis."""
        if user_overwrite:
            self.cleanup()
        self.create_job_file()
        self.create_input_file()
        recycled = self.run_external(overwrite=user_overwrite)
        if recycled and not self.completed:
            recycled = self.run_external(overwrite=True)
        if not self.completed:
            raise ExternalError("External job could not be completed (%s)" % self.filename)
        self.read_output()
        self.ran = True
        return recycled


class SimpleJob(Job):
    def __init__(self, prefix, title):
        Job.__init__(self, prefix, title)
        self.summary = {}

    def determine_completed(self):
        self.summarize_output()

    def read_output(self):
        self.process_output_summary()

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
            "gawk -f %sinterfaces/awk/%s.awk < %s.out > %s.out.smr" % (
                context.share_path, self.awk_scriptname(),
                self.filename, self.filename
            )
        )
        smr = file("%s.out.smr" % self.filename)
        self.summary.update(eval(''.join(smr)))
        smr.close()
        self.assign_fields(["completed"])
        
    def assign_fields(self, fields):
        for field in fields:
            self.__dict__[field] = self.summary[field]

    def process_output_summary(self):
        """Process the attributes taken from the summary file and assigned to self."""
        raise NotImplementedError


def reload_job(job_filename):
    f = file(job_filename, 'r')
    job = load(f)
    f.close()
    return job
