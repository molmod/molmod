// MolMod is a collection of molecular modelling tools for python.
// Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
// for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
// reserved unless otherwise stated.
//
// This file is part of MolMod.
//
// MolMod is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// MolMod is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
// --



#ifndef MOLMOD_COMMON_H
#define MOLMOD_COMMON_H


void dot_matrix_vector_ddd(double *matrix, double *in, double *out);
void dot_matrix_vector_did(double *matrix, long *in, double *out);
void dot_matrixT_vector_ddd(double *matrix, double *in, double *out);
void dot_matrixT_vector_did(double *matrix, long *in, double *out);
double distance(double *a, double *b);
double distance_periodic(double *a, double *b, double *matrix, double *reciprocal);
double distance_delta(double *a, double *b, double *delta);
double distance_delta_periodic(double *a, double *b, double *delta, double *matrix, double *reciprocal);
double norm(double *a);

#endif
