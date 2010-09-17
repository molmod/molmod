Unit conversion in MolMod
=========================

.. highlight:: python
   :linenothreshold: 5

The MolMod package expresses all quantities internally in a well-defined unit
system, i.e. using `internal units`. Whenever numerical data goes into the
MolMod code, it is assumed to be in atomic units, and results will also be returned
in atomic units. When data is read from a file (through one of the modules
in :mod:`molmod.io`) the numerical values are immediately converted to `internal
units`. In case one needs to print a number on screen or to a file, it can be
converted back to some `non-internal unit`.

The `internal units` in MolMod are the atomic units. Most atomic units are
listed at the `NIST website for physical constants
<http://physics.nist.gov/cuu/Constants/index.html>`_. Note that there is no
atomic unit for temperature, so we simply use Kelvin as the `internal unit`
of temperature. The atomic units system has the advantage, just
like the SI system, to be a consistent unit system. This means that, given any
mathematical expression, one can just plug in all values in atomic units and the
result will come out in atomic units without the need for some conversion
factor. The atomic units are used instead of SI units as most molecular
quantities have a reasonable order of magnitude in atomic units. A second
advantage is that several important constant (which are basically plain
conversion factors) become unity.

The conversion factors to translate values from and to `internal units` are
defined in the module :mod:`molmod.units`. The convention for unit conversion is
that multiplication with a constant converts a value into `internal units` and
that division with the same constant converts a value back to the `non-internal
unit`::

    from molmod import *

    distance = 5*angstrom
    print "Distances in internal units:", distance
    print "Distances in Angstrom:", distance/angstrom

Some physical constants that are relevant for molecular modeling, are defined in
:mod:`molmod.constants`.

A few practical examples of computations with unit conversion are given below.

Energies
~~~~~~~~

This example shows, given a reactant and product energy, how one computes the
reaction energy and prints out the result in kJ/mol.

File: ``examples/000_units/a_reaction.py``

.. literalinclude:: ../examples/000_units/a_reaction.py

This example is very basic, but it demonstrates some general requirements for
any script that uses the MolMod package:

* The first line always reads ``#!/usr/bin/env python``. This only matters when
  the script is executed on a Unix machine. It is the so-called `shebang line
  <http://en.wikipedia.org/wiki/Shebang_%28Unix%29>`_, which is used to
  determine the interpreter when the script is executed on the command line::

      toon@poony ~> ./a_reaction.py

  is equivalent to::

      toon@poony ~> /usr/bin/env python a_reaction.py

  The short form only works when the executable flag of the script file is set.
  If this is not the case yet, it can be changed as follows::

      toon@poony ~> chmod +x a_reaction.py

* The beginning of the file must contain a line ``from molmod import *`` to
  import the entire molmod package. Some alternative methods to import
  libraries in Python are discussed in the `official Python documentation
  <http://docs.python.org/>`_.

Wavenumber
~~~~~~~~~~

Given the mass and force constant of a harmonic spring the following script
computes the spectroscopic wavenumber of the oscillation. The mass and the force
constant are an approximate model for a C-H bond.

File: ``examples/000_units/b_chbond.py``

.. literalinclude:: ../examples/000_units/b_chbond.py


Rotational partition function of Hydrogen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The last example solves a typical statistical physics exam question: "Compute
the probability that a Hydrogen molecule in a dilute gas does not rotate at a
temperature of 300 K. Approximate the Hydrogen molecule as a rigid rotor". The
parameters are included in the source code below. Note that the constants
``planck`` and ``boltzmann`` are defined in :mod:`molmod.constants`.

File: ``examples/000_units/c_h2rot.py``

.. literalinclude:: ../examples/000_units/c_h2rot.py
