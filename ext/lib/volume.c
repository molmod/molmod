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


static inline double distance_sq(double *p1, double *p2) {
  double tmp, d2;
  tmp = p1[0]-p2[0];
  d2 = tmp*tmp;
  tmp = p1[1]-p2[1];
  d2 += tmp*tmp;
  tmp = p1[2]-p2[2];
  d2 += tmp*tmp;
  return d2;
}


static inline double distance(double *p1, double *p2) {
  return sqrt(distance_sq(p1, p2));
}


int in_spheres(int n, double* probe, double probe_radius, double *spheres, double *sphere_radii) {
  int i;
  double tmp;
  for (i=n-1; i>=0; i--) {
    tmp = (*sphere_radii + probe_radius);
    tmp = distance_sq(spheres, probe) - tmp*tmp;
    if (tmp < -1e-10) {
      return -1;
    }
    spheres += 3;
    sphere_radii++;
  }
  return 0;
}


int in_spheres_all(int n, double* probe, double probe_radius, double *spheres, double *sphere_radii, int *hits) {
  int i, counter;
  double tmp;

  /* collect all overlapping atoms */
  counter = 0;
  for (i=n-1; i>=0; i--) {
    tmp = (*sphere_radii + probe_radius);
    tmp = distance_sq(spheres, probe) - tmp*tmp;
    if (tmp < -1e-10) {
      counter += 1;
      *hits = n-i-1;
      hits++;
    }
    spheres += 3;
    sphere_radii++;
  }

  return counter;
}


void center_ses1(double *probe, double *close1, double close1_radius, double *center) {
  /* compute the center of the sphere that obeys the following conditions:
       - the probe sphere with the new center is tangent to the sphere close1
       - the new center is colinear with the original probe and the center of
         the sphere close1.
     The implementation below is based on the solution of an equivalent but
     simpler problem: find the point on the sphere with a radius equal to the
     sum of probe and close1 radius that also obeys the colinearity condition.
  */

  double ratio;
  ratio = close1_radius/distance(probe, close1);
  center[0] = close1[0] + (probe[0] - close1[0])*ratio;
  center[1] = close1[1] + (probe[1] - close1[1])*ratio;
  center[2] = close1[2] + (probe[2] - close1[2])*ratio;
}


int center_ses2(double *probe, double *close1, double close1_radius, double *close2, double close2_radius, double *center) {
  /* compute the center of the sphere that obeys the following conditions:
       - the probe sphere with the new center is tangent to the sphere close1
         and tangent to the sphere close2
       - the new center is coplanar with the original probe and the centers of
         close1 and close2.
     The implementation below is based on the solution of an equivalent but
     simpler problem, in analogy with center_ses1.
  */
  double r, r1, dot, para[3], ortho[3];

  r = distance(close1, close2);
  if (r > close1_radius + close2_radius) goto fail;
  if (r < fabs(close1_radius - close2_radius)) goto fail;

  /* the intersection between the line through centers of close1 and close2 and
     the plane that contains the intersection of the spheres close1 and close2. */
  r1 = (r*r + close1_radius*close1_radius - close2_radius*close2_radius)/(2*r);
  dot = r1/r;
  center[0] = close1[0] + dot*(close2[0] - close1[0]);
  center[1] = close1[1] + dot*(close2[1] - close1[1]);
  center[2] = close1[2] + dot*(close2[2] - close1[2]);

  /* a unit vector parallel to the unit plane and parallel to the plane formed
     by the two centers close1 and close 2 and the probe. */
  para[0] = probe[0] - close1[0];
  para[1] = probe[1] - close1[1];
  para[2] = probe[2] - close1[2];
  ortho[0] = (close2[0] - close1[0])/r;
  ortho[1] = (close2[1] - close1[1])/r;
  ortho[2] = (close2[2] - close1[2])/r;
  dot = para[0]*ortho[0] + para[1]*ortho[1] + para[2]*ortho[2];
  para[0] -= dot*ortho[0];
  para[1] -= dot*ortho[1];
  para[2] -= dot*ortho[2];
  dot = sqrt(para[0]*para[0] + para[1]*para[1] + para[2]*para[2]);
  para[0] /= dot;
  para[1] /= dot;
  para[2] /= dot;

  r1 = sqrt(close1_radius*close1_radius - r1*r1);
  center[0] += r1*para[0];
  center[1] += r1*para[1];
  center[2] += r1*para[2];
  return 0;
fail:
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;
  return -1;
}


static int get_plane(double *p1, double r1, double *p2, double r2, double *normal, double *r, double *c) {
  normal[0] = p2[0] - p1[0];
  normal[1] = p2[1] - p1[1];
  normal[2] = p2[2] - p1[2];
  *r = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  if (*r <= fabs(r1 - r2)) return -1;
  if (*r > r1+r2) return -1;
  normal[0] /= *r;
  normal[1] /= *r;
  normal[2] /= *r;
  *c = normal[0]*p1[0] +
       normal[1]*p1[1] +
       normal[2]*p1[2] +
       (*r*(*r) + r1*r1 - r2*r2)/(2*(*r));
  /*printf("normal = [%f, %f, %f], norm = %f, coeff = %f\n", normal[0], normal[1], normal[2], *r, *c);*/
  return 0;
}


int center_ses3(double *probe, double *close1, double close1_radius, double *close2, double close2_radius, double *close3, double close3_radius, double *center) {
  /* compute the center of the sphere that obeys the following conditions:
       - the probe sphere with the new center is tangent to the sphere close1,
         tangent to the sphere close2 and tangent to sphere close3.
     The implementation below is based on the solution of an equivalent but
     simpler problem, in analogy with center_ses1.
  */

  double normal12[3], normal23[3], normal[3];
  double r12, r23, r, c12, c23, c, det;

  /* the two intersection planes: close1-close2 and close2-close3 */
  if (get_plane(close1, close1_radius, close2, close2_radius, normal12, &r12, &c12) == -1) goto fail;
  if (get_plane(close2, close2_radius, close3, close3_radius, normal23, &r23, &c23) == -1) goto fail;
  /* the plane through the three sphere centers */
  normal[0] = normal12[2]*normal23[1] - normal12[1]*normal23[2];
  normal[1] = normal12[0]*normal23[2] - normal12[2]*normal23[0];
  normal[2] = normal12[1]*normal23[0] - normal12[0]*normal23[1];
  c = normal[0]*close2[0] + normal[1]*close2[1] + normal[2]*close2[2];

  /* the intersection point of the three planes */
  det = normal12[0]*(normal23[1]*normal[2] - normal23[2]*normal[1]) +
        normal12[1]*(normal23[2]*normal[0] - normal23[0]*normal[2]) +
        normal12[2]*(normal23[0]*normal[1] - normal23[1]*normal[0]);
  if (det == 0) goto fail;
  /*printf("det = %f\n", det);*/
  center[0] = ((normal23[1]*normal[2] - normal23[2]*normal[1])*c12 + (normal12[2]*normal[1] - normal12[1]*normal[2])*c23 + (normal12[1]*normal23[2] - normal12[2]*normal23[1])*c)/det;
  center[1] = ((normal23[2]*normal[0] - normal23[0]*normal[2])*c12 + (normal12[0]*normal[2] - normal12[2]*normal[0])*c23 + (normal12[2]*normal23[0] - normal12[0]*normal23[2])*c)/det;
  center[2] = ((normal23[0]*normal[1] - normal23[1]*normal[0])*c12 + (normal12[1]*normal[0] - normal12[0]*normal[1])*c23 + (normal12[0]*normal23[1] - normal12[1]*normal23[0])*c)/det;
  /*printf("center = [%f, %f, %f]\n", center[0], center[1], center[2]);*/

  /* the remaining distance */
  r = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  if (normal[0]*(probe[0]-center[0]) + normal[1]*(probe[1]-center[1]) + normal[2]*(probe[2]-center[2]) < 0) {
    r *= -1; /* make sure it points towards the probe */
  }
  normal[0] /= r;
  normal[1] /= r;
  normal[2] /= r;
  /*printf("normal = [%f, %f, %f]\n", normal[0], normal[1], normal[2]);*/

  /* displace the center to the final position */
  r = close2_radius*close2_radius - distance_sq(close2, center);
  /*printf("difference = %f\n", r);*/
  if (r<0) {
    /*printf("radius1-distance=%f\n", close1_radius - distance(close1, center));
    printf("radius2-distance=%f\n", close2_radius - distance(close2, center));
    printf("radius3-distance=%f\n", close3_radius - distance(close3, center));
    printf("radius1+radius2-distance=%f\n", close1_radius + close2_radius - distance(close1, close2));
    printf("radius2+radius3-distance=%f\n", close2_radius + close3_radius - distance(close2, close3));
    printf("radius3+radius1-distance=%f\n", close3_radius + close1_radius - distance(close3, close1));

    printf("distance1c=%f\n", distance(close1, center));
    printf("distance2c=%f\n", distance(close2, center));
    printf("distance3c=%f\n", distance(close3, center));
    printf("distance12=%f\n", distance(close1, close2));
    printf("distance23=%f\n", distance(close2, close3));
    printf("distance31=%f\n", distance(close3, close1));
    printf("error of last kind.\n");*/
    goto fail;
  }
  r = sqrt(r);
  center[0] += normal[0]*r;
  center[1] += normal[1]*r;
  center[2] += normal[2]*r;

  return 0;
fail:
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;
  return -1;
}


double monte_carlo_volumes(int n, double probe_radius, double *spheres, double *sphere_radii, int num_iter, int bigbox, int ses_iter, long int *counts) {
  /* use monte carlo sampling the integrate the van der waals volume, the
     volume enclosed by the accessible surface and the solvent excluded surface.
     */

  double low[3], size[3], probe[3], center[3], tmp;
  double close1_radius, close2_radius, close3_radius;
  double *sphere, *close1, *close2, *close3, *radius;
  int i, j1, j2, j3, num_hits;
  int *hits;

  /*** 0. allocate work memory ***/
  hits = malloc(n*sizeof(int));
  if (hits==NULL) return -1.0;

  /*** A. construct a bounding box ***/
  low[0] = low[1] = low[2] = size[0] = size[1] = size[2] = 0;
  sphere = spheres;
  radius = sphere_radii;
  for (i=n-1; i>=0; i--) {
    tmp = sphere[0] - *radius; if (low[0] > tmp) low[0] = tmp;
    tmp = sphere[1] - *radius; if (low[1] > tmp) low[1] = tmp;
    tmp = sphere[2] - *radius; if (low[2] > tmp) low[2] = tmp;
    tmp = sphere[0] + *radius; if (size[0] < tmp) size[0] = tmp;
    tmp = sphere[1] + *radius; if (size[1] < tmp) size[1] = tmp;
    tmp = sphere[2] + *radius; if (size[2] < tmp) size[2] = tmp;
    radius += 1;
    sphere += 3;
  }
  if (bigbox!=0) bigbox = 1;
  low[0] -= bigbox*probe_radius;
  low[1] -= bigbox*probe_radius;
  low[2] -= bigbox*probe_radius;
  size[0] += bigbox*probe_radius - low[0];
  size[1] += bigbox*probe_radius - low[1];
  size[2] += bigbox*probe_radius - low[2];

  /*** B. Monte carlo ***/
  for (i=num_iter-1; i>=0; i--) {
    probe[0] = (double)rand() / (double)(RAND_MAX) * size[0] + low[0];
    probe[1] = (double)rand() / (double)(RAND_MAX) * size[1] + low[1];
    probe[2] = (double)rand() / (double)(RAND_MAX) * size[2] + low[2];
    if (in_spheres(n, probe, 0.0, spheres, sphere_radii)) {
      /* in the vdw volume */
      counts[0]++;
      counts[1]++;
      counts[2]++;
    } else {
      num_hits = in_spheres_all(n, probe, probe_radius, spheres, sphere_radii, hits);
      if (num_hits > 0) {
        for (j1=0; j1 < num_hits; j1++) {
          close1 = spheres + 3*hits[j1];
          close1_radius = sphere_radii[hits[j1]] + probe_radius;
          center_ses1(probe, close1, close1_radius, center);
          if (in_spheres(n, center, probe_radius, spheres, sphere_radii)==0) goto next; /* not in the ses volume */
        }
        for (j1=0; j1<num_hits; j1++) {
          close1 = spheres + 3*hits[j1];
          close1_radius = sphere_radii[hits[j1]] + probe_radius;
          for (j2=0; j2<j1; j2++) {
            close2 = spheres + 3*hits[j2];
            close2_radius = sphere_radii[hits[j2]] + probe_radius;
            if (center_ses2(probe, close1, close1_radius, close2, close2_radius, center)==-1) {
              printf("This should never happen 1\n");
              return -1.0;
            }
            if (distance(probe, center) - probe_radius < -1e-10) {
              if (in_spheres(n, center, probe_radius, spheres, sphere_radii)==0) goto next; /* not in the ses volume */
            }
          }
        }
        for (j1=0; j1<num_hits; j1++) {
          close1 = spheres + 3*hits[j1];
          close1_radius = sphere_radii[hits[j1]] + probe_radius;
          for (j2=0; j2<j1; j2++) {
            close2 = spheres + 3*hits[j2];
            close2_radius = sphere_radii[hits[j2]] + probe_radius;
            for (j3=0; j3<j2; j3++) {
              close3 = spheres + 3*hits[j3];
              close3_radius = sphere_radii[hits[j3]] + probe_radius;
              if (center_ses3(probe, close1, close1_radius, close2, close2_radius, close3, close3_radius, center)==-1) {
                /* sometimes three overlapping spheres have no common
                   intersection point, but an entire common volume. In such case
                   there is no useful solutions. Just continue */
                continue;
              }
              if (distance(probe, center) - probe_radius < -1e-10) {
                if (in_spheres(n, center, probe_radius, spheres, sphere_radii)==0) goto next; /* not in the ses volume */
              }
            }
          }
        }
        /* if we get here, the point lies in the ses volume. */
        counts[2]++;
        /*printf("radius=%8.3f   vdw=%8.3f    over=%8.3f   z=%8.3f\n",
          sqrt(probe[0]*probe[0]+probe[1]*probe[1]),
          sphere_radii[0],
          distance(probe, center) - probe_radius,
          probe[2]
        );*/
next:
        /* in the sas volume */
        counts[1]++;
      }
    }
  }
  counts[4] += num_iter;

  /*** -1. free the memory ***/
  free(hits);

  return size[0]*size[1]*size[2];
}
