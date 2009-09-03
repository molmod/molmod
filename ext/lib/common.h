// MolModExt implements a few number crunching routines for the molmod package in C.
// Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
//
// This file is part of MolModExt.
//
// MolModExt is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// MolModExt is distributed in the hope that it will be useful,
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
void dot_matrix_vector_did(double *matrix, int *in, double *out);
void dot_matrixT_vector_ddd(double *matrix, double *in, double *out);
void dot_matrixT_vector_did(double *matrix, int *in, double *out);
double distance(double *a, double *b);
double norm(double *a);

#endif
