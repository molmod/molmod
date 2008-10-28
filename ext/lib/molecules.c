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

void molecules_distance_matrix(int n, double *cor, double* dm) {
  int i,j;
  double d,tmp;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      tmp = cor[3*i  ] - cor[3*j  ];
      d = tmp*tmp;
      tmp = cor[3*i+1] - cor[3*j+1];
      d += tmp*tmp;
      tmp = cor[3*i+2] - cor[3*j+2];
      d += tmp*tmp;
      d = sqrt(d);
      dm[i*n+j] = d;
      dm[j*n+i] = d;
    }
  }
}

