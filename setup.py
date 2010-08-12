#!/usr/bin/env python
# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
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


import glob, os
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
from numpy.distutils.command.install_data import install_data


class MyInstallData(install_data):
    """Add a datadir.txt file that points to the root for the data files. It is
       otherwise impossible to figure out the location of these data files at
       runtime.
    """
    def run(self):
        # Do the normal install_data
        install_data.run(self)
        # Create the file datadir.txt. It's exact content is only known
        # at installation time.
        dist = self.distribution
        libdir = dist.command_obj["install_lib"].install_dir
        for name in dist.packages:
            if '.' not in name:
                destination = os.path.join(libdir, name, "datadir.txt")
                print "Creating %s" % destination
                if not self.dry_run:
                    f = file(destination, "w")
                    print >> f, self.install_dir
                    f.close()


setup(
    name='MolMod',
    version='0.004',
    description='MolMod is a collection of molecular modelling tools for python.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    cmdclass={'install_data': MyInstallData},
    package_dir = {'molmod': 'molmod'},
    packages=[
        'molmod',
        'molmod.io',
    ],
    data_files=[
        ('share/molmod', [
            "share/periodic.csv", "share/bonds.csv",
            "share/mass.mas03", "share/nubtab03.asc",
            "share/toyff_angles.txt"
        ]),
    ],
    ext_modules=[
        Extension("molmod.ext", ["molmod/ext.pyf", "molmod/common.c",
            "molmod/ff.c", "molmod/graphs.c", "molmod/similarity.c", "molmod/molecules.c",
            "molmod/unit_cells.c", "molmod/volume.c",
        ]),
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)


