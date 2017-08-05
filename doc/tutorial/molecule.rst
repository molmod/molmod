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

Working with molecules
======================

.. highlight:: python
   :linenothreshold: 5


Introduction
~~~~~~~~~~~~


Molecule objects
----------------

The class :class:`molmod.molecules.Molecule` is a fundamental part of the MolMod
package. Instances of this class hold all information to describe the geometry
and the topology of a molecular structure. In practice the following attributes
can be stored in a Molecule object:

* ``numbers``  --  An array of atomic numbers with N elements. This
  is the only mandatory attribute. All other attributes are optional.
* ``coordinates``  --  An array with atomic coordinates (in Bohr) with
  shape (N,3).
* ``title``  --  A short description of the molecule.
* ``masses``  --  An array with atomic masses (in :math:`m_e`)
* ``graph``  --  A ``MolecularGraph`` instance describing the connectivity of
  the atoms.
* ``symbols``  --  A tuple with atomic symbols. These may be the names of the
  elements or some force-field atom types.
* ``unit_cell``  --  A ``UnitCell`` instance describing the periodic boundary
  conditions.

One can create Molecule instances based on several conventional file formats as
follows::

    mol = Molecule.from_file("foo.xyz")

The extension of the file is used to detect the file format. The following
formats are currently supported: ``*.cml``, ``*.fchk``, ``*.pdb``, ``*.sdf``,
``*.xyz``. It is relatively trivial to implement new file formats, but there is
no real urgent need for that because the babel program can convert most files
into one of these formats.

A Molecule object can be written to file as follows::

    mol.write_to_file("bar.xyz")

The currently supported formats for writing files are: ``*.xyz``, ``*.cml``

Derived quantities
------------------

Some derived quantities are accessible as attributes. The contents of these
attributes are only computed when they are accessed for the first time. Once a
derived quantity is computed it is stored internally (in the molecule object) in
case it would be requested a second time. This mechanism is called `caching`. It
imposes good scripting practices when using the MolMod package and also improves
the efficiency of the scripts. These are a few convenient derived quantities:

* ``chemical_formula``  --  A string containing the chemical formula.
* ``com``  --  The center of mass of the molecule.
* ``distance_matrix``  --  A square matrix with all inter-atomic distances in
  a molecule (ignoring periodic boundary conditions for the moment).
* ``inertia_tensor``  --  The inertia tensor of the molecule computed with the
  center of mass as the origin.
* ``mass``  --  The total mass of the system.
* ``size``  --  The number of atoms.

An up-to-date list can be found in the reference documentation:
:class:`molmod.molecules.Molecule`. Most derived quantities depend on optional
attributes. If these optional attributes are not present, an error is raised
when a derived quantity is requested.


Modifying molecule objects
--------------------------

The `caching` mechanism has one important side-effect. When the underlying data
of the derived quantities change, the `cached` results are no longer correct.
This problem is called `cache inconsistency`. There are typically two ways to
impose `cache consistency`: (i) all underlying quantities are kept constant, or
(ii) cached results are deleted when the underlying data change. The first
mechanism is used throughout the MolMod package. When you need to modify a
``Molecule`` object, just create a new one with modified attributes. This is
facilitated by the ``copy_with`` method::

    mol1 = Molecule.from_file("foo.xyz")
    mol2 = mol1.copy_with(title="bar")

Any of the above attributes can be modified through ``copy_with``. Cached
derived quantities are never discarded with this approach.

The concept of derived quantities and read-only attributes is extensively used
in the following classes:

* ``Molecule``, see :class:`molmod.molecules.Molecule`.
* ``Graph``, see :class:`molmod.graphs.Graph`.
* ``MolecularGraph``, see :class:`molmod.molecular_graphs.MolecularGraph`.
* ``UnitCell``, see :class:`molmod.unit_cells.UnitCell`.
* ``Translation``, see :class:`molmod.transformations.Translation`.
* ``Rotation``, see :class:`molmod.transformations.Rotation`.
* ``Complete``, see :class:`molmod.transformations.Complete`.


Examples
~~~~~~~~


Conversion of files
-------------------

Although the conversion of files with molecular systems from one format to the
other is easily done with the babel program, the example is still useful for
didactic purposes.

File: ``molmod/examples/001_molecules/a_convert.py``

.. literalinclude:: ../../molmod/examples/001_molecules/a_convert.py


Center of mass to origin
------------------------

This example shows how a modified molecule object is made.

File: ``molmod/examples/001_molecules/b_com.py``

.. literalinclude:: ../../molmod/examples/001_molecules/b_com.py


All carbons
-----------

The following program prints the X-, Y- and Z-coordinates of all Carbon atoms in
the ibuprofen molecule.

File: ``molmod/examples/001_molecules/c_carbon.py``

.. literalinclude:: ../../molmod/examples/001_molecules/c_carbon.py


Size of a molecule
------------------

One can define the size of a molecule in several ways. The following program
prints out the largest inter-atomic distance in Ångström.

File: ``molmod/examples/001_molecules/d_size.py``

.. literalinclude:: ../../molmod/examples/001_molecules/d_size.py


Shape of a molecule
-------------------

One can define the shape of a molecule based on the covariance of the atomic
coordinates with respect to their geometric center. Some shape categories were
introduced by `Zingg <http://en.wikipedia.org/wiki/Equidimensional>`_. The
following program computes the shape category for the ibuprofen molecule

File: ``molmod/examples/001_molecules/e_shape.py``

.. literalinclude:: ../../molmod/examples/001_molecules/e_shape.py


Problems
~~~~~~~~

* Write a program that makes a list of all (unique) pairs of Carbon atoms that
  have an inter-atomic distance below 2.0 Ångström.

* Try to find for each shape category (equant, prolate, oblate and bladed) some
  molecules from `PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_. Just download a
  bunch of SDF files and use the glob module from Python to loop over all these
  files::

      from glob import glob
      for fn in glob("*.sdf"):
          mol = Molecule.from_file(fn)
          # do something here

  If you feel `lucky`, extend to program so that it downloads the first 100
  molecules from the PubChem database.
