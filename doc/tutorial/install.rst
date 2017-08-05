..
    : MolMod is a collection of molecular modelling tools for python.
    : Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
    : for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
    : reserved unless otherwise stated.
    :
    : This file is part of MolMod.
    :
    : MolMod is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : MolMod is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Installation instructions
#########################


Disclaimer
==========

MolMod is developed and tested on modern Linux environments. The installation and usage
will therefore be relatively easy on Linux. If you want to use MolMod on other operating
systems such as Windows or OSX, you should have a minimal computer geek status to get it
working. We are always interested in hearing from your installation adventures.


Dependencies
============

The following software is used by MolMod:

* Python >=2.7 (including the development files): http://www.python.org/doc/
* A C compiler e.g. gcc: http://gcc.gnu.org/
* Numpy >=1.0 or later: http://numpy.scipy.org/
* Cython >=0.24.1: http://cython.org/
* Nosetests >=0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/


Installation
============

You can install MolMod with pip, using either of the following two commands:

.. code:: bash

    # system wide (requires root permission) or in virtual env
    pip install numpy Cython
    pip install molmod

    # installs in ~/.local
    pip install numpy Cython --user
    pip install molmod --user

Alternatively, you can use conda. (See https://www.continuum.io/downloads)

.. code:: bash

    conda install -c tovrstra molmod


Testing
=======

The installation can be tested as follows:

.. code:: bash

    nosetests molmod

This will run a series of tests to check the validity of the outcomes generated
by MolMod. If some tests fail, post an issue on https://github.com/molmod/molmod/issues
