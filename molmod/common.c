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


#include <math.h>
#include "common.h"


void dot_matrix_vector_ddd(double *matrix, double *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[1]*in[1] + matrix[2]*in[2];
  out[1] = matrix[3]*in[0] + matrix[4]*in[1] + matrix[5]*in[2];
  out[2] = matrix[6]*in[0] + matrix[7]*in[1] + matrix[8]*in[2];
}

void dot_matrix_vector_did(double *matrix, long *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[1]*in[1] + matrix[2]*in[2];
  out[1] = matrix[3]*in[0] + matrix[4]*in[1] + matrix[5]*in[2];
  out[2] = matrix[6]*in[0] + matrix[7]*in[1] + matrix[8]*in[2];
}

void dot_matrixT_vector_ddd(double *matrix, double *in, double *out) {
  out[0] = matrix[0]*in[0] + matrix[3]*in[1] + matrix[6]*in[2];
  out[1] = matrix[1]*in[0] + matrix[4]*in[1] + matrix[7]*in[2];
  out[2] = matrix[2]*in[0] + matrix[5]*in[1] + matrix[8]*in[2];
}

void dot_matrixT_vector_did(double *matrix, long *in, double *out) {
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

double distance_periodic(double *a, double *b, double *matrix, double *reciprocal) {
  double delta[3];
  return distance_delta_periodic(a, b, delta, matrix, reciprocal);
}

double distance_delta(double *a, double *b, double *delta) {
  delta[0] = a[0] - b[0];
  delta[1] = a[1] - b[1];
  delta[2] = a[2] - b[2];
  return sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
}

double distance_delta_periodic(double *a, double *b, double *delta, double *matrix, double *reciprocal) {
  double fractional[3], delta_cell[3];
  delta[0] = a[0] - b[0];
  delta[1] = a[1] - b[1];
  delta[2] = a[2] - b[2];
  dot_matrixT_vector_ddd(reciprocal, delta, fractional);
  //printf("delta[0]=%f  delta[1]=%f  delta[2]=%f\n", delta[0], delta[1], delta[2]);
  //printf("fractional[0]=%f  fractional[1]=%f  fractional[2]=%f\n", fractional[0], fractional[1], fractional[2]);
  fractional[0] = floor(fractional[0] + 0.5);
  fractional[1] = floor(fractional[1] + 0.5);
  fractional[2] = floor(fractional[2] + 0.5);
  dot_matrix_vector_ddd(matrix, fractional, delta_cell);
  delta[0] -= delta_cell[0];
  delta[1] -= delta_cell[1];
  delta[2] -= delta_cell[2];
  //printf("delta[0]=%f  delta[1]=%f  delta[2]=%f\n", delta[0], delta[1], delta[2]);
  //printf("fractional[0]=%f  fractional[1]=%f  fractional[2]=%f\n\n", fractional[0], fractional[1], fractional[2]);
  return sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
}

double norm(double *a) {
  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
