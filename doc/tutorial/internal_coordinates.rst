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
Each internal coordinate function follows the same `API
<http://en.wikipedia.org/wiki/Application_programming_interface>`_ style, i.e.
::

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

File: ``molmod/examples/003_internal_coordinates/a_bond_length.py``

.. literalinclude:: ../../molmod/examples/003_internal_coordinates/a_bond_length.py


All bending angles in dopamine
------------------------------

We use the graph features to detect all bending angles in the dopamine molecule.
An overview is printed including the atomic indexes and elements involved in
each angle.

File: ``molmod/examples/003_internal_coordinates/b_bending_angles.py``

.. literalinclude:: ../../molmod/examples/003_internal_coordinates/b_bending_angles.py



Simple Force-field Hessian
--------------------------

In this example we do a simple computation of a force-field Hessian for the
propane molecule. The force-field in this example model has only bond-stretch
and bending angle terms. We further assume that all internal coordinates are at
the rest-value of the corresponding energy term in the force field. This
facilitates the program as we only have to take into account the force constants
of the energy terms.

The mathematical form of the FF model is:

.. math:: E_{\text{FF}} =
            \sum_{i=1}^{M_\text{bonds}} K_{i,\text{bond}} (b_i - b_{i,0})^2 +
            \sum_{i=1}^{M_\text{bends}} K_{i,\text{bend}} (\theta_i - \theta_{i,0})^2.

The Cartesian gradient of the force-field energy becomes:

.. math:: \frac{\partial E_{\text{FF}}}{\partial x_{j}} =
            \sum_{i=1}^{M_\text{bonds}} 2 K_{i,\text{bond}}
                (b_i - b_{i,0}) \frac{\partial b_i}{\partial x_{j}} +
            \sum_{i=1}^{M_\text{bends}} 2 K_{i,\text{bend}}
                (\theta_i - \theta_{i,0}) \frac{\partial \theta_i}{\partial x_{j}}

The Cartesian Hessian of the force-field becomes:

.. math:: \frac{\partial^2 E_{\text{FF}}}{\partial x_{j_1} \partial x_{j_2}} =
            \sum_{i=1}^{M_\text{bonds}} 2 K_{i,\text{bond}} \left[
                \frac{\partial b_i}{\partial x_{j_1}} \frac{\partial b_i}{\partial x_{j_2}} +
                (b_i - b_{i,0}) \frac{\partial^2 b}{\partial x_{j_1} \partial x_{j_2}}
            \right] +
            \sum_{i=1}^{M_\text{bends}} 2 K_{i,\text{bend}} \left[
                \frac{\partial \theta_i}{\partial x_{j_1}} \frac{\partial \theta_i}{\partial x_{j_2}} +
                (\theta_i - \theta_{i,0}) \frac{\partial^2 \theta_i}{\partial x_{j_1} \partial x_{j_2}} +
            \right]

For the sake of a simple example, we assume that all bond lengths and valence
angles are at their optimum such that the second derivatives of the internal
coordinates towards the Cartesian coordinates drop out of the expression for the
Hessian:

.. math:: \frac{\partial^2 E_{\text{FF}}}{\partial x_{j_1} \partial x_{j_2}} \approx
            \sum_{i=1}^{M_\text{bonds}} 2 K_{i,\text{bond}}
                \frac{\partial b_i}{\partial x_{j_1}} \frac{\partial b_i}{\partial x_{j_2}} +
            \sum_{i=1}^{M_\text{bends}} 2 K_{i,\text{bend}}
                \frac{\partial \theta_i}{\partial x_{j_1}} \frac{\partial \theta_i}{\partial x_{j_2}} +

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

File: ``molmod/examples/003_internal_coordinates/c_ff_hessian.py``

.. literalinclude:: ../../molmod/examples/003_internal_coordinates/c_ff_hessian.py


DFT Hessian in `internal coordinates` -- Stationary geometries
--------------------------------------------------------------

Given a DFT Hessian for a stationary molecular system, and a list of internal
coordinates, :math:`q_i`, one may construct a model for the Hessian based on a
complete second order expansion in terms of internal coordinates. This means
that one uses a force-field model with all possible cross terms as follows:

.. math:: E_{\text{FF}} = \sum_{i_1=1}^M
                          \sum_{i_2=1}^M
                          K_{i_1 i_2} (q_{i_1} - q_{i_1,0}) (q_{i_2} - q_{i_2,0})

One can define internal force constants, :math:`K_{ij}` such that the
force-field Hessian coincides with a given Hessian from a DFT computation. One
could call the matrix :math:`K` the Hessian in `internal coordinates`, but as
we will see below, this not a very strict definition. The purpose of such a
transformation is that the matrix :math:`K` is (or can be made) more diagonally
dominant than the original Hessian. It is a tool to get insight in the
contributions to the Hessian that one should include in a force-field model.
Similar models are used in redundant optimization methods.

Assuming that the geometry is at a stationary point of the potential energy
surface, one has to following expression for the force-field Hessian:

.. math:: \frac{\partial^2 E_{\text{FF}}}{\partial x_{j_1} \partial x_{j_2}}
            = \sum_{i_1=1}^M \sum_{i_2=1}^M
            K_{i_1 i_2} \left(
                \frac{\partial q_{i_1}}{\partial x_{j_1}}
                \frac{\partial q_{i_2}}{\partial x_{j_2}}
                +\frac{\partial q_{i_1}}{\partial x_{j_2}}
                \frac{\partial q_{i_2}}{\partial x_{j_1}}
            \right)
            = 2\sum_{i_1=1}^M \sum_{i_2=1}^M
            K_{i_1 i_2}
            \frac{\partial q_{i_1}}{\partial x_{j_1}}
            \frac{\partial q_{i_2}}{\partial x_{j_2}}

The last step assumes that :math:`K_{i_1 i_2} = K_{i_2 i_1}`. This can be
rewritten in matrix notation:

.. math:: H = J K J^T
    :label: matrix_hessian_int

with

.. math::
    :nowrap:

    \begin{align*}
        H_{j_1 j_2} & = \frac{\partial^2 E_{\text{FF}}}{\partial x_{j_1} \partial x_{j_2}} \\
        J_{j_1 i_i} & = \frac{\partial q_{i_1}}{\partial x_{j_1}}
    \end{align*}

In practice, the Jacobian matrix :math:`J` is rectangular, mainly because there
are often much more internal coordinates than Cartesian coordinates. Therefore
one can not simply invert :eq:`matrix_hessian_int` to obtain the Hessian in
internal coordinates. There are actually many matrices :math:`K` that solve
equation :eq:`matrix_hessian_int`, and one has to introduce some additional
criteria to fix :math:`K`.

For the sake of simplicity, we will use the `Moore-Penrose pseudoinverse
<http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse>`_ of the
Jacobian to invert :eq:`matrix_hessian_int`, but there may be better choices.
This may seem a unique choice, but in practice it is not. The problem is that
the columns of the Jacobian can have different units. Therefore the numerical
Jacobian and its pseudoinverse depend on the choice of the units. Again, we make
a simple choice here to solve this issue: all columns where the internal
coordinates are some sort of distance, are kept as is. All angles are converted
to a length unit with a fixed conversion factor: c = 1Å / 5°. Again, there may
be better choices, e.g. one may normalize the columns of the Jacobian.

The program below performs such an inversion on a DFT Hessian of the for the
dopamine molecule, using all bond lengths, bending angles and dihedral angles.
The script is written in an object-oriented style: each internal coordinate is
an object of the class ``BondLength``, ``BendingAngle`` or ``DihdralAngle``.
These three classes are derived from ``InternalCoordinates`` and share a common
`API <http://en.wikipedia.org/wiki/Application_programming_interface>`_ such
that the main program does not have to worry about the nature of the internal
coordinates.

The Hessian is computed with a Gaussian03 B3LYP/6-31G(d) frequency job. The
result is stored in the file ``dopamine.fchk``, which is a stripped version of the
file generated by the ``fromchk`` program that is part of Gaussian03.

File: ``molmod/examples/003_internal_coordinates/d_dft_hessian.py``

.. literalinclude:: ../../molmod/examples/003_internal_coordinates/d_dft_hessian.py

For the computation of the pseudo-inverse, there is one more gotcha's: some
internal coordinates are exactly redundant, i.e. in the case of dopamine, the
the dihedral angles in the aromatic ring are a bit problematic. It is comparable
to the situation of ethene, where one dihdral angle can always be written as a
simple linear function of the three other ones. The columns in the Jacobian
corresponding to these dihedrals are linearly dependent. One can in principle
leave out 6 columns in the case of dopamine. However, any selection of internal
coordinates to be removed would be a subjective choice. One can avoid such
subjective input by dropping the almost-zero singular values during he
computation of the generalized inverse. Therefore the second argument to
``pinv`` in the script is set to ``1e-5``, which means that all singular values
that are 100000 times smaller than the largest one, are treated as if they were
zeros.


Problems
~~~~~~~~


All bond CH-bond lengths in propane
-----------------------------------

Write a program that computes all CH-bond lengths in the propane molecule. Use
the same style as in the example that computes the bending angles in the
dopamine molecule.


More complete Force-field Hessian
---------------------------------

Modify the program that computes the Hessian of propane in such a way that it
takes into account the contributions due to deviations of internal coordinates
from the rest value of each force-field term.


DFT Hessian in `internal coordinates` -- Non-stationary geometries
------------------------------------------------------------------

Extend the program for the computation of the Hessian in internal coordinates
such that it also works for non-stationery points on the potential energy
surface. The starting point is a more general expression for the force-field
model

.. math:: E_{\text{FF}} =
            \sum_{i=1}^M A (q_{i} - q_{i,0}) +
            \sum_{i_1=1}^M \sum_{i_2=1}^M
            K_{i_1 i_2} (q_{i_1} - q_{i_1,0}) (q_{i_2} - q_{i_2,0})

where the rest values can be the current internal coordinates, or some reference
values. The vector :math:`A` and the matrix :math:`K` can be defined in such
a way that the force-field gradient and Hessian in Cartesian coordinates
coincide with their DFT counter-parts. First derive proper forms for :math:`A`
and :math:`K`.

Let the program still use the current internal coordinates as reference values
for the force field. Let it also print out the vector :math:`A`.


Fitting force-constants
-----------------------

Write a program that fits the force-constants of the propane molecule based on a
DFT Hessian computation on the ground-state geometry, for given values of the
rest lengths and angles. Test to what extent these force constants depend on the
choice of rest parameters.
