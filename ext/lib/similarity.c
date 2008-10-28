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

void similarity_table_labels(int n, int *labels, int* labels_table) {
  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      (*labels_table) = labels[i];
      labels_table++;
      (*labels_table) = labels[j];
      labels_table++;
    }
  }
}

void similarity_table_distances(int n, double *distance_matrix, double *distances_table) {
  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      (*distances_table) = distance_matrix[i*n+j];
      distances_table++;
    }
  }
}

#define PAIR_LESS(l1,l2) (l1[0]==l2[0]?(l1[1]<l2[1]):(l1[0]<l2[0]))
#define PAIR_GREATER(l1,l2) (l1[0]==l2[0]?(l1[1]>l2[1]):(l1[0]>l2[0]))
#define PAIR_EQUAL(l1,l2) (l1[0]==l2[0]&&l1[1]==l2[1])

double similarity_measure(int n1, int *labels1, double *distances1, int n2, int *labels2, double *distances2, double margin, double cutoff) {
  double result, dav, delta;
  int i1, i2, c2;
  int *c_labels2;
  double *c_distances2;

  result = 0.0;

  i2 = 0;
  for (i1=0;i1<n1;i1++) {
    if (PAIR_LESS(labels1,labels2)) {
      goto next_iter; // goto hell if you don't like goto.
    }
    while (PAIR_GREATER(labels1,labels2) && i2 < n2) {
      labels2 += 2;
      distances2++;
      i2++;
    }
    if (i2 >= n2) {
      break; // end of the second table is reached.
    }
    if (PAIR_LESS(labels1,labels2)) {
      goto next_iter; // goto hell if you don't like goto.
    }
    //printf("Starting with (%i,%i) (%i,%i)\n", labels1[0], labels1[1], labels2[0], labels2[1]);

    c_labels2 = labels2;
    c_distances2 = distances2;
    c2 = i2;
    while (PAIR_EQUAL(labels1,c_labels2) && c2 < n2) {
      dav = 0.5*(*distances1 + *c_distances2);
      //printf("dav=%f    cutoff=%f\n", dav, cutoff);
      if (dav < cutoff) {
        delta = fabs(*distances1 - *c_distances2);
        //printf("delta=%f    rmargin=%f\n", delta, margin);
        if (delta < margin) {
          result += (1-dav/cutoff)*0.5*(cos(delta/margin/M_PI)+1);
          //printf("delta=%f    result=%f    dav=%f    margin=%f    cutoff=%f    M_PI=%f\n", delta, result, dav, margin, cutoff, M_PI);
        }
      }
      c_labels2 += 2;
      c_distances2++;
      c2++;
    }
    //printf("Did %i\n", c2-i2);

next_iter:
    // iterate to the next pair in table 1
    labels1 += 2;
    distances1++;
  }

  return result;
}
