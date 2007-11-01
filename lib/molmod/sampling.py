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
# --


class Job(object):
    # A empty job class. This is usefull when custom jobs are constructed,
    # that are going to be pickled.
    pass


class Ensemble(object):
    def __init__(self, name):
        self.name = name
        self.jobs = []

    def add_job(self, job):
        # 1) assert that the jobs have the required attributes
        # 2) store the job in a sensible place
        self.jobs.append(job)

    def finalize(self):
        # post process the jobs. (sorting, ...)
        pass


class Sampler(object):
    EnsembleClass = None

    def __init__(self, name):
        self.name = name
        self.ensemble = self.EnsembleClass(name)

    def run(self):
        # should return an ensemble object
        raise NotImplementedError



