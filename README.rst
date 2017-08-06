.. image:: https://travis-ci.org/molmod/molmod.svg?branch=master
    :target: https://travis-ci.org/molmod/molmod
.. image:: https://ci.appveyor.com/api/projects/status/g3kpxx1n80vek3hq/branch/master?svg=true
    :target: https://ci.appveyor.com/project/tovrstra/molmod
.. image:: https://anaconda.org/tovrstra/molmod/badges/version.svg
    :target: https://anaconda.org/tovrstra/molmod

MolMod is a collection of molecular modelling tools for python. It is used by other
software developed at the CMM, including Yaff, TAMkin and Zeobuilder.

More information about MolMod can be found on the CMM Code website:
http://molmod.ugent.be/software

MolMod is distributed as open source software under the conditions of the GPL
license version 3. Read the file COPYING for more details, or visit
http://www.gnu.org/licenses/


Installation
============

MolMod can be installed with pip (system wide or in a virtual environment):

.. code:: bash

    pip install Cython numpy
    pip install molmod

Alternatively, you can install MolMod in your home directory:

.. code:: bash

    pip install Cython numpy --user
    pip install molmod --user

Lastly, you can also install MolMod with conda. (See
https://www.continuum.io/downloads)

.. code:: bash

    conda install -c tovrstra molmod


Testing
=======

The tests can be executed as follows:

.. code:: bash

    nosetests molmod
