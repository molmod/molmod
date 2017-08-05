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


#include "ff.h"

#include <math.h>
#include <stdlib.h>
#include "common.h"

void add_grad(
  size_t i, size_t j, double s, double *cor, double *delta,
  double *gradient
) {
  gradient[i*3  ] += s*delta[0];
  gradient[j*3  ] -= s*delta[0];
  gradient[i*3+1] += s*delta[1];
  gradient[j*3+1] -= s*delta[1];
  gradient[i*3+2] += s*delta[2];
  gradient[j*3+2] -= s*delta[2];
}

double ff_dm_quad(
  size_t natom, int periodic, double *cor, double *dm0, double *dmk,
  double amp, double *gradient, double *matrix, double *reciprocal
) {
  size_t i, j;
  double delta[3], d, d0, k, tmp, result;

  result = 0.0;
  //printf("natom=%i\n", natom);
  for (i=0; i<natom; i++) {
    for (j=0; j<i; j++) {
      d0 = dm0[i*natom+j];
      k = dmk[i*natom+j];
      //printf("i=%i  j=%i  d0=%i\n", i,j,d0);
      if (d0>0) {
        if (periodic) {
          d = distance_delta_periodic(cor + 3*i, cor + 3*j, delta, matrix, reciprocal);
        } else {
          d = distance_delta(cor + 3*i, cor + 3*j, delta);
        }
        tmp = (d-d0);
        result += amp*k*tmp*tmp;
        if (gradient!=NULL) {
          //tmp = 2*amp*tmp/d0*radii[i]*radii[j]/d;
          tmp = 2*amp*k*tmp/d;
          add_grad(i, j, tmp, cor, delta, gradient);
        }
        //result += tmp*tmp;
      }
    }
  }
  //printf("result=%f\n", result);
  return result;
}


double ff_dm_reci(
  size_t natom, int periodic, double *cor, double *radii, long *dm0,
  double amp, double *gradient, double *matrix, double *reciprocal
) {
  size_t i, j;
  double delta[3], d, r0, tmp, result;

  result = 0.0;
  for (i=0; i<natom; i++) {
    for (j=0; j<i; j++) {
      if (dm0[i*natom+j]>1) {
        if (periodic) {
          d = distance_delta_periodic(cor + 3*i, cor + 3*j, delta, matrix, reciprocal);
        } else {
          d = distance_delta(cor + 3*i, cor + 3*j, delta);
        }
        r0 = radii[i]+radii[j];
        if (d < r0) {
            d /= r0;
            result += amp*(d-1)*(d-1)/d;
            if (gradient!=NULL) {
              tmp = amp*(1-1/d/d)/r0/d/r0;
              add_grad(i, j, tmp, cor, delta, gradient);
            }
        }
      }
    }
  }
  return result;
}


double ff_bond_quad(
  size_t npair, int periodic, double *cor, long *pairs, double *lengths,
  double amp, double *gradient, double *matrix, double *reciprocal
) {
  size_t b, i, j;
  double delta[3], result, d, tmp;

  result = 0.0;
  for (b=0; b<npair; b++) {
    i = pairs[2*b  ];
    j = pairs[2*b+1];
    if (periodic) {
      d = distance_delta_periodic(cor + 3*i, cor + 3*j, delta, matrix, reciprocal);
    } else {
      d = distance_delta(cor + 3*i, cor + 3*j, delta);
    }

    tmp = d-lengths[b];
    result += amp*tmp*tmp;
    if (gradient!=NULL) {
      tmp = 2*amp*tmp/d;
      add_grad(i, j, tmp, cor, delta, gradient);
    }
    //printf("result=%f\n", result);
  }
  return result;
}

double ff_bond_hyper(
  size_t npair, int periodic, double *cor, long *pairs, double *lengths,
  double scale, double amp, double *gradient, double *matrix, double *reciprocal
) {
  size_t b, i, j;
  double delta[3], result, d, tmp;

  result = 0.0;
  for (b=0; b<npair; b++) {
    i = pairs[2*b  ];
    j = pairs[2*b+1];
    if (periodic) {
      d = distance_delta_periodic(cor + 3*i, cor + 3*j, delta, matrix, reciprocal);
    } else {
      d = distance_delta(cor + 3*i, cor + 3*j, delta);
    }

    tmp = d-lengths[b];
    result += amp*(cosh(scale*tmp)-1);
    if (gradient!=NULL) {
      tmp = amp*scale*sinh(scale*tmp)/d;
      add_grad(i, j, tmp, cor, delta, gradient);
    }
  }
  return result;
}
