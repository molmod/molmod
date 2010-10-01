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



Simple Force-field Hessian
--------------------------

In this example we do a simple computation of a force-field Hessian for the
propane molecule. The force-field in this example model has only bond-stretch
and bending angle terms. We further assume that all internal coordinates are at
the rest-value of the corresponding energy term in the force field. This
facilitates the program as we only have to take into account the force constants
of the energy terms.

The mathematical form of the FF model is:

.. math:: E_{\text{FF}} = \sum_{i=1}^{N_\text{bonds}} K_i (b - b_0)^2 +
                          \sum_{i=1}^{N_\text{bends}} K_i (\theta - \theta_0)^2.

The Cartesian gradient of the force-field energy becomes:

.. math:: \frac{\partial E_{\text{FF}}}{\partial x_{j}} =
            \sum_{i=1}^{N_\text{bonds}} 2 K_i (b - b_0) \frac{\partial b}{\partial x_{j}} +
            \sum_{i=1}^{N_\text{bends}} 2 K_i (\theta - \theta_0) \frac{\partial \theta}{\partial x_{j}}

The Cartesian Hessian of the force-field becomes:

.. math:: \frac{\partial^2 E_{\text{FF}}}{\partial x_{j1} \partial x_{j2}} =
            \sum_{i=1}^{N_\text{bonds}} 2 K_i \left[
                \frac{\partial b}{\partial x_{j1}} \frac{\partial b}{\partial x_{j2}} +
                (b - b_0) \frac{\partial^2 b}{\partial x_{j1} \partial x_{j2}}
            \right] +
            \sum_{i=1}^{N_\text{bends}} 2 K_i \left[
                \frac{\partial \theta}{\partial x_{j1}} \frac{\partial \theta}{\partial x_{j2}} +
                (\theta - \theta_0) \frac{\partial^2 \theta}{\partial x_{j1} \partial x_{j2}} +
            \right]

For the sake of a simple example, we assume that all bond lengths and valence
angles are at their optimum such that the second derivatives of the internal
coordinates towards the Cartesian coordinates drop out of the expression for the
Hessian:

.. math:: \frac{\partial^2 E_{\text{FF}}}{\partial x_{j1} \partial x_{j2}} \approx
            \sum_{i=1}^{N_\text{bonds}} 2 K_i
                \frac{\partial b}{\partial x_{j1}} \frac{\partial b}{\partial x_{j2}} +
            \sum_{i=1}^{N_\text{bends}} 2 K_i
                \frac{\partial \theta}{\partial x_{j1}} \frac{\partial \theta}{\partial x_{j2}} +

We use the following force constants.

============== ================= ======================================
Energy term    Force constant    Unit
============== ================= ======================================
CH bond        310               k cal mol\ :sup:`-1` Å\ :sup:`-2`
CC bond        220               k cal mol\ :sup:`-1` Å\ :sup:`-2`
HCH bend       35                k cal mol\ :sup:`-1` rad\ :sup:`-2`
HCC bend       30                k cal mol\ :sup:`-1` rad\ :sup:`-2`
CCC bend       60                k cal mol\ :sup:`-1` rad\ :sup:`-2`
============== ================= ======================================

The program below uses an object-oriented approach to implement the force-field
model. Each energy term is conceived as an object of either the
``BondStretchTerm`` or the ``BendingAngleTerm`` class. They both derive from
the ``HarmonicEnergyTerm`` class where all the Hessian logic is implemented.
Each term object contains attributes for the force-field parameters and the atom
indexes that are involved in the internal coordinate. The term objects are kept
in a list in the ``ForceField`` class that has a method to compute the Hessian
for a given geometry.

File: ``examples/003_internal_coordinates/c_ff_hessian.py``

.. literalinclude:: ../../examples/003_internal_coordinates/c_ff_hessian.py


DFT Hessian in `internal coordinates`
-------------------------------------


Problems
~~~~~~~~


All bond CH-bond lengths in propane
-----------------------------------

Write a program that computes all CH-bond lengths in the propane molecule. Use
the same style as in the example that computes the bending angles in the
dopamine molecule.


Complete Force-field Hessian
----------------------------

Modify the program that computes the Hessian of propane in such a way that it
takes into account the contributions due to deviations from the rest value of
each force-field term.


Fitting force-constants
-----------------------
