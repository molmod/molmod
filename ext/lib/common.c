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


#include <math.h>
#include "common.h"


void dot_matrix_vector_ddd(double *matrix, double *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[1]*in[1] + matrix[2]*in[2];
  out[1] = matrix[3]*in[0] + matrix[4]*in[1] + matrix[5]*in[2];
  out[2] = matrix[6]*in[0] + matrix[7]*in[1] + matrix[8]*in[2];
}

void dot_matrix_vector_did(double *matrix, int *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[1]*in[1] + matrix[2]*in[2];
  out[1] = matrix[3]*in[0] + matrix[4]*in[1] + matrix[5]*in[2];
  out[2] = matrix[6]*in[0] + matrix[7]*in[1] + matrix[8]*in[2];
}

void dot_matrixT_vector_ddd(double *matrix, double *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[3]*in[1] + matrix[6]*in[2];
  out[1] = matrix[1]*in[0] + matrix[4]*in[1] + matrix[7]*in[2];
  out[2] = matrix[2]*in[0] + matrix[5]*in[1] + matrix[8]*in[2];
}

void dot_matrixT_vector_did(double *matrix, int *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[3]*in[1] + matrix[6]*in[2];
  out[1] = matrix[1]*in[0] + matrix[4]*in[1] + matrix[7]*in[2];
  out[2] = matrix[2]*in[0] + matrix[5]*in[1] + matrix[8]*in[2];
}

double distance(double *a, double *b) {
  double tmp, dsq;
  tmp = a[0]-b[0];
  dsq = tmp*tmp;
  tmp = a[1]-b[1];
  dsq += tmp*tmp;
  tmp = a[2]-b[2];
  dsq += tmp*tmp;
  return sqrt(dsq);
}

double norm(double *a) {
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

