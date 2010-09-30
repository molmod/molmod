Internal coordinates
====================

.. highlight:: python
   :linenothreshold: 5


Introduction
~~~~~~~~~~~~

Internal coordinates are often used to characterize molecular geometries and for
the definition of energy terms in valence force-field models. In general an
internal coordinate can be defined as any function of the Cartesian atomic
coordinates that does not depend on global rotation and translation.

The module :mod:`molmod.ic` contains functions for the
internal coordinates listed in the table below.

==================== =================== =========================================================
Name                 Number of arguments Description
==================== =================== =========================================================
``bond_length``      2                   Distance between two atoms
``pair_distance``    2                   (idem)
``bend_cos``         3                   The cosine of a bending angle
``bend_angle``       4                   The bending angle
``dihed_cos``        4                   The cosine of a dihedral angle
``dihed_angle``      4                   The dihedral angle, following IUPAC standard for the sign
``opbend_angle``     4                   The out of plane bending angle
==================== =================== =========================================================

In addition to the value of
the internal coordinate, it can also compute the first and second order
derivatives of the internal coordinate towards the Cartesian coordinates.
Each internal coordinate function follows the same API style, i.e. ::

    def some_ic(v1, v2, ..., deriv=0):
        ...

All mandatory arguments are 3-element numpy arrays with the atomic coordinates
that define the internal coordinate. The last and optional argument determines
which derivatives are computed. By default no derivatives are computed, and the
internal coordinate is returned in a singleton tuple. When ``deriv==1`` the
internal coordinate and the gradient are returned in a tuple. When ``deriv==2``
the internal coordinate, the gradient and the Hessian are returned in a tuple.

The implementation is based on a generic schemes that takes care of following
the chain rule for the derivatives, which makes it rather easy to extend this
module with new types of internal coordinates.

The examples and problems below use a propane and dopamine geometry. The latter
is optimized at the B3LYP/6-31G(d) level, followed by a frequency computation.
The computations are carried out with Gaussian03. The formatted checkpoint file
of the frequency job is stripped to include only the section used in this
chapter.


Examples
~~~~~~~~


A simple bond length
--------------------

This is just a simple example...

File: ``examples/003_internal_coordinates/a_bond_length.py``

.. literalinclude:: ../../examples/003_internal_coordinates/a_bond_length.py


All bending angles in dopamine
------------------------------

We use the graph features to detect all bending angles in the dopamine molecule.
An overview is printed including the atomic indexes and elements involved in
each angle.

File: ``examples/003_internal_coordinates/b_bending_angles.py``

.. literalinclude:: ../../examples/003_internal_coordinates/b_bending_angles.py



Force-field Hessian -- Optimized molecule
-----------------------------------------


Hartree-Fock Hessian in internal coordinates
--------------------------------------------


Problems
~~~~~~~~


All bond CH-bond lengths in propane
-----------------------------------


Force-field Hessian -- non-optimized molecule
---------------------------------------------


Fitting force-constants
-----------------------
