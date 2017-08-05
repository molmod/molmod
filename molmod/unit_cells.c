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


#include "unit_cells.h"

#include <math.h>
#include <stdlib.h>

#include "common.h"


size_t unit_cell_get_radius_indexes(double *matrix, double *reciprocal, double radius,
                                    long *max_ranges, long *indexes) {
  long lowlim[3], uplim[3], index[3], side[3];
  long i0, i1, i2, j0, j1, j2;
  size_t counter, k, sum;
  double center[3], corner[3], pos[3], frac[3];
  double scale, d;

  for (k = 0; k < 3; k++) {
    /* See UnitCell.get_radius_ranges(self, radius) in the python code */
    uplim[k] = ceil(radius*sqrt(
        reciprocal[k  ]*reciprocal[k  ] +
        reciprocal[k+3]*reciprocal[k+3] +
        reciprocal[k+6]*reciprocal[k+6]
    ));
    lowlim[k] = -uplim[k];
    /* reduce the ranges a priori */
    if ((max_ranges[k] > 0) && (uplim[k]*2+1 > max_ranges[k])) {
      if (max_ranges[k]%2 == 0) {
        uplim[k] = max_ranges[k]/2-1;
        lowlim[k] = -max_ranges[k]/2;
      } else {
        uplim[k] = max_ranges[k]/2;
        lowlim[k] = -max_ranges[k]/2;
      }
      /*printf(" k=%i   max_ranges[k]=%i   lowlim[k]=%i   uplim[k]=%i\n", k, max_ranges[k], lowlim[k], uplim[k]);*/
    }
  }
  counter = 0;
  /* iterator over all neigboring cells in the range: i0, i1, i2 */
  for (i0 = lowlim[0]; i0 <= uplim[0]; i0++) {
    index[0] = i0;
    for (i1 = lowlim[1]; i1 <= uplim[1]; i1++) {
      index[1] = i1;
      for (i2 = lowlim[2]; i2 <= 0; i2++) {
        index[2] = i2;
        /* The center of a neighboring cell */
        dot_matrix_vector_did(matrix, index, center);

        /* Find the shortest distance between a point in the central cell
           and the point in a neighboring cell.

           The problem can be simplified to finding the shortest distance
           from the origin to a point in a double unit cell box centered
           around the relative vector between the two unit cell boxes under
           scrunity. Therefore the relative vector is called 'center' in
           this routine.

           This is a constrained optimization problem that would in principle
           be solved with the active set algorithm. We prefer to keep the
           implementation simple and will try all reasonable combinations of
           constraints.

           There are 6 constraints, i.e. the six faces of the box, which
           naively leads to 2**6 possible combinations of turning constraints
           of and on. Constrainging the solution to two oposite faces at the
           same time does not make sense. Consequently, there are only 27
           reasonable combinations of constraints.

           For each combination, only those that have a solution in, or on
           the edge of the box are considered. From the latter, we test if
           it is within the radius. */

        /* Iterator over all possible combinations of constraints: j0, j1, j2

           ji = -1: constraint to face i on negative side of the box
           ji =  0: no constraint to face i
           ji = +1: constraint to face i on negative side of the box
        */
        for (j0 = -1; j0 <= 1; j0++) {
          side[0] = j0;
          for (j1 = -1; j1 <= 1; j1++) {
            side[1] = j1;
            for (j2 = -1; j2 <= 1; j2++) {
              side[2] = j2;
              /* the point on the boundary of the double unit cell which is
                 common to all active constraints. */
              dot_matrix_vector_did(matrix, side, corner);
              corner[0] += center[0];
              corner[1] += center[1];
              corner[2] += center[2];
              /* the number of active constraints */
              sum = abs(j0) + abs(j1) + abs(j2);
              if (sum == 3) {
                /* The constraints fully determine the solution. */
                pos[0] = corner[0];
                pos[1] = corner[1];
                pos[2] = corner[2];
              } else if (sum == 2) {
                /* Two constraints, line */
                k = 0; /* column k of matrix the unconstrained direction) */
                if (j1 == 0) k = 1;
                if (j2 == 0) k = 2;
                /* the closest point to the origin on this line */
                scale = (
                    (corner[0]*matrix[k] + corner[1]*matrix[k+3] + corner[2]*matrix[k+6])/
                    (matrix[k]*matrix[k] + matrix[k+3]*matrix[k+3] + matrix[k+6]*matrix[k+6])
                );
                pos[0] = corner[0] - matrix[k]*scale;
                pos[1] = corner[1] - matrix[k+3]*scale;
                pos[2] = corner[2] - matrix[k+6]*scale;
              } else if (sum == 1) {
                /* One constraint, plane */
                k = 0; /* column k with the unconstrained direction */
                if (j1 != 0) k = 1;
                if (j2 != 0) k = 2;
                /* the closest point to the origin on this plane */
                scale = (
                    (corner[0]*reciprocal[k] + corner[1]*reciprocal[k+3] + corner[2]*reciprocal[k+6])/
                    (reciprocal[k]*reciprocal[k] + reciprocal[k+3]*reciprocal[k+3] + reciprocal[k+6]*reciprocal[k+6])
                );
                pos[0] = reciprocal[k]*scale;
                pos[1] = reciprocal[k+3]*scale;
                pos[2] = reciprocal[k+6]*scale;
              } else {
                /* No constraints */
                pos[0] = 0.0;
                pos[1] = 0.0;
                pos[2] = 0.0;
              }
              /* The solution in fractional coordinates, to check if the
                 solution is insde or on the edge of the box */
              dot_matrixT_vector_ddd(reciprocal, pos, frac);
              k = 0; /* k = 1 means within radius, recycle variable */
              if ((abs(frac[0] - i0) < 1.000001) && (abs(frac[1] - i1) < 1.000001) && (abs(frac[2] - i2) < 1.000001)) {
                d = norm(pos);
                k = d < radius;
              }
              if (k) break;
            } /* j2 */
            if (k) break;
          } /* j1 */
          if (k) break;
        } /* j0 */
        if (k) {
          indexes[counter] = i0;
          indexes[counter+1] = i1;
          indexes[counter+2] = i2;
          counter += 3;
          if ((i2 != 0) && (-i2 <= uplim[2])) {
            /* since the resulting array contents is point symmetric and we
               have only half a loop over i2 */
            indexes[counter] = -i0;
            indexes[counter+1] = -i1;
            indexes[counter+2] = -i2;
            counter += 3;
          }
        }
      } /* i2 */
    } /* i1 */
  } /* i0 */
  return counter/3;
}
