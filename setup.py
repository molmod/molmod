#!/usr/bin/env python
# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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
#--


import os
from glob import glob
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
from distutils.command.install_data import install_data


class MyInstallData(install_data):
    """Ensure data is store at root-level
    """
    def run(self):
        # Do the normal install_data
        install_data.run(self)

setup(
    name='molmod',
    version='1.0',
    description='MolMod is a collection of molecular modelling tools for python.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    cmdclass={'install_data': MyInstallData},
    package_dir = {'molmod': 'molmod'},
    packages=[
        'molmod', 'molmod.test',
        'molmod.io', 'molmod.io.test',
    ],
    data_files=[
        ('share/molmod', [
            "data/periodic.csv", "data/bonds.csv",
            "data/mass.mas03", "data/nubtab03.asc",
            "data/toyff_angles.txt"
        ]),
        ('share/molmod/test', glob('data/test/*')),
    ] + [
        ('share/molmod/examples/%s' % os.path.basename(dn), glob('%s/*.*' % dn))
        for dn in glob('data/examples/???_*')
    ],
    ext_modules=[
        Extension("molmod.ext", ["molmod/ext.pyf", "molmod/common.c",
            "molmod/ff.c", "molmod/graphs.c", "molmod/similarity.c",
            "molmod/molecules.c", "molmod/unit_cells.c",
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


