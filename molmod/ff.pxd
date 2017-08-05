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
# --


cdef extern from "ff.h":
    double ff_dm_quad(
      size_t natom, int periodic, double *cor, double *dm0, double *dmk,
      double amp, double *gradient, double *matrix, double *reciprocal
    )

    double ff_dm_reci(
      size_t natom, int periodic, double *cor, double *radii, long *dm0,
      double amp, double *gradient, double *matrix, double *reciprocal
    )

    double ff_bond_quad(
      size_t npair, int periodic, double *cor, long *pairs, double *lengths,
      double amp, double *gradient, double *matrix, double *reciprocal
    )

    double ff_bond_hyper(
      size_t npair, int periodic, double *cor, long *pairs, double *lengths,
      double scale, double amp, double *gradient, double *matrix, double *reciprocal
    )
