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
#include <stdlib.h>

inline void add_grad(int i, int j, double s, double *cor, double *gradient) {
  gradient[i*3  ] += s*(cor[i*3  ]-cor[j*3  ]);
  gradient[j*3  ] -= s*(cor[i*3  ]-cor[j*3  ]);
  gradient[i*3+1] += s*(cor[i*3+1]-cor[j*3+1]);
  gradient[j*3+1] -= s*(cor[i*3+1]-cor[j*3+1]);
  gradient[i*3+2] += s*(cor[i*3+2]-cor[j*3+2]);
  gradient[j*3+2] -= s*(cor[i*3+2]-cor[j*3+2]);

}

double ff_dm_quad(int n, double *radii, double *cor, int *dm0, double amp, double *gradient) {
  int i,j,d0;
  double d, tmp, result;

  result = 0.0;
  //printf("n=%i\n", n);
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      d0 = dm0[i*n+j];
      //printf("i=%i  j=%i  d0=%i\n", i,j,d0);
      if (d0>0) {
        tmp = cor[i*3  ]-cor[j*3  ];
        d = tmp*tmp;
        tmp = cor[i*3+1]-cor[j*3+1];
        d += tmp*tmp;
        tmp = cor[i*3+2]-cor[j*3+2];
        d += tmp*tmp;
        d = sqrt(d);

        //printf("d=%f    d0=%i\n", d, d0);
        tmp = (d-d0*(radii[i]+radii[j]));
        result += amp*tmp*tmp*radii[i]*radii[j]/d0;
        if (gradient!=NULL) {
          tmp = 2*amp*tmp/d0*radii[i]*radii[j]/d;
          add_grad(i, j, tmp, cor, gradient);
        }
        //result += tmp*tmp;
      }
    }
  }
  //printf("result=%f\n", result);
  return result;
}


double ff_dm_reci(int n, double *radii, double *cor, int *dm0, double amp, double *gradient) {
  int i,j,d0,r0;
  double d, tmp, result;

  result = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      d0 = dm0[i*n+j];
      if (d0>1) {
        tmp = cor[i*3  ]-cor[j*3  ];
        d = tmp*tmp;
        tmp = cor[i*3+1]-cor[j*3+1];
        d += tmp*tmp;
        tmp = cor[i*3+2]-cor[j*3+2];
        d += tmp*tmp;
        d = sqrt(d);
        r0 = radii[i]+radii[j];
        if (d < r0) {
            d /= r0;
            result += amp*(d-1)*(d-1)/d/d0;
            if (gradient!=NULL) {
              tmp = amp*(1-1/d/d)/r0/d/d0;
              add_grad(i, j, tmp, cor, gradient);
            }
        }
      }
    }
  }
  return result;
}


double ff_bond_quad(int m, int n, double *cor, int *pairs, double *lengths, double amp, double *gradient) {
  int b, i, j;
  double result, d, tmp;

  result = 0.0;
  //printf("m=%i\n", m);
  for (b=0; b<m; b++) {
    i = pairs[2*b  ];
    j = pairs[2*b+1];
    tmp = cor[i*3  ]-cor[j*3  ];
    d = tmp*tmp;
    tmp = cor[i*3+1]-cor[j*3+1];
    d += tmp*tmp;
    tmp = cor[i*3+2]-cor[j*3+2];
    d += tmp*tmp;
    d = sqrt(d);

    tmp = d-lengths[b];
    result += amp*tmp*tmp;
    if (gradient!=NULL) {
      tmp = 2*amp*tmp/d;
      add_grad(i, j, tmp, cor, gradient);
    }
    //printf("result=%f\n", result);
  }
  return result;
}

double ff_bond_hyper(int m, int n, double *cor, int *pairs, double *lengths, double scale, double amp, double *gradient) {
  int b, i, j;
  double result, d, tmp;

  result = 0.0;
  for (b=0; b<m; b++) {
    i = pairs[2*b  ];
    j = pairs[2*b+1];
    tmp = cor[i*3  ]-cor[j*3  ];
    d = tmp*tmp;
    tmp = cor[i*3+1]-cor[j*3+1];
    d += tmp*tmp;
    tmp = cor[i*3+2]-cor[j*3+2];
    d += tmp*tmp;
    d = sqrt(d);

    tmp = d-lengths[b];
    result += amp*(cosh(scale*tmp)-1);
    if (gradient!=NULL) {
      tmp = amp*scale*sinh(scale*tmp)/d;
      add_grad(i, j, tmp, cor, gradient);
    }
  }
  return result;
}


