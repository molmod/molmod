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

* Python >=2.7, <3.0 (including the development files): http://www.python.org/doc/
* A C compiler e.g. gcc: http://gcc.gnu.org/
* Numpy >=1.0 or later: http://numpy.scipy.org/
* Cython >=0.24.1: http://cython.org/
* Nosetests >=0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/


Installation
============

You can install MolMod with pip, using either of the following two commands:

.. code:: bash

    pip install molmod   # system wide (requires root permission) or in virtual env
    pip install molmod --user   # installs in ~/.local

Alternatively, you can use conda. (See https://www.continuum.io/downloads)

.. code:: bash

    conda install -c tovrstra molmod


Testing
=======

The installation can be tested as follows:

.. code:: bash

    (cd; nosetests -v molmod)

This will run a series of tests to check the validity of the outcomes generated
by MolMod. If some tests fail, post an issue on https://github.com/molmod/molmod/issues
