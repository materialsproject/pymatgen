/* Copyright (C) 2015 Atsushi Togo */
/* All rights reserved. */

/* This file is part of niggli. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the phonopy project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "niggli.h"

static double A, B, C, eta, xi, zeta, eps;
static int l, m, n;
static double *tmat = NULL;
static double *lattice = NULL;

static int initialize(const double *lattice_, const double eps_);
static void finalize(double *lattice_);
static void reset(void);
static void step0(void);
static int step1(void);
static int step2(void);
static int step3(void);
static int step4(void);
static int step5(void);
static int step6(void);
static int step7(void);
static int step8(void);
static void set_parameters(void);
static void set_angle_types(void);
static double * get_transpose(const double *M);
static double * get_metric(const double *M);
static double * multiply_matrices(const double *A, const double *B);

#ifdef NIGGLI_DEBUG
#define debug_print(...) printf(__VA_ARGS__)
static void debug_show(void);
static void debug_show(void)
{
  int i;
  printf("%f %f %f %f %f %f\n", A, B, C, xi, eta, zeta);
  printf("%d %d %d\n", l, m, n);
  
  for (i = 0; i < 3; i++) {
    printf("%f %f %f\n", lattice[i * 3], lattice[i * 3 + 1], lattice[i * 3 + 2]);
  }
}
#else
#define debug_print(...)
#define debug_show(...)
#endif

#ifdef NIGGLI_WARNING
#include <stdio.h>
#define warning_print(...) fprintf(stderr,__VA_ARGS__)
#else
#define warning_print(...)
#endif

/* return 0 if failed */
int niggli_reduce(double *lattice_, const double eps_)
{
  int i;

  if (! initialize(lattice_, eps_)) {
    return 0;
  }

  step0();
  
  for (i = 0; i < 10; i++) {
    if (step1()) {
      debug_print("step1\n");
      debug_show();
      debug_print("\n");
    }

    if (step2()) {
      debug_print("step2\n");
      debug_show();
      debug_print("\n");
      continue;
    }

    if (step3()) {
      debug_print("step3\n");
      debug_show();
      debug_print("\n");
    }

    if (step4()) {
      debug_print("step4\n");
      debug_show();
      debug_print("\n");
    }

    if (step5()) {
      debug_print("step5\n");
      debug_show();
      debug_print("\n");
      continue;
    }

    if (step6()) {
      debug_print("step6\n");
      debug_show();
      debug_print("\n");
      continue;
    }

    if (step7()) {
      debug_print("step7\n");
      debug_show();
      debug_print("\n");
      continue;
    }

    if (step8()) {
      debug_print("step7\n");
      debug_show();
      debug_print("\n");
      continue;
    }

    break;
  }

  finalize(lattice_);
  return 1;
}

static int initialize(const double *lattice_, const double eps_)
{
  if ((tmat = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    return 0;
  }

  eps = eps_;
  if ((lattice = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    free(tmat);
    tmat = NULL;
    return 0;
  }
  memcpy(lattice, lattice_, sizeof(double) * 9);

  return 1;
}

static void finalize(double *lattice_)
{
  free(tmat);
  tmat = NULL;
  memcpy(lattice_, lattice, sizeof(double) * 9);
  free(lattice);
  lattice = NULL;
}

static void reset(void)
{
  double *lat_tmp;

  lat_tmp = NULL;
  
  lat_tmp = multiply_matrices(lattice, tmat);

  memcpy(lattice, lat_tmp, sizeof(double) * 9);
  step0();
  free(lat_tmp);
  lat_tmp = NULL;
}

static void step0(void)
{
  set_parameters();
  set_angle_types();
}

static int step1(void)
{
  if (A > B + eps ||
      (! (fabs(A -B) > eps) && fabs(xi) > fabs(eta) + eps)) {
    tmat[0] = 0,  tmat[1] = -1, tmat[2] = 0;
    tmat[3] = -1, tmat[4] = 0,  tmat[5] = 0;
    tmat[6] = 0,  tmat[7] = 0,  tmat[8] = -1;
    reset();
    return 1;
  }
  else {return 0;}
}

static int step2(void)
{
  if (B > C + eps ||
      (! (fabs(B - C) > eps) && fabs(eta) > fabs(zeta) + eps)) {
    tmat[0] = -1, tmat[1] = 0,  tmat[2] = 0;
    tmat[3] = 0,  tmat[4] = 0,  tmat[5] = -1;
    tmat[6] = 0,  tmat[7] = -1, tmat[8] = 0;
    reset();
    return 1;
  }
  else {return 0;}
}

static int step3(void)
{
  int i, j, k;
  if (l * m * n == 1) {
    if (l == -1) {i = -1;} else {i = 1;}
    if (m == -1) {j = -1;} else {j = 1;}
    if (n == -1) {k = -1;} else {k = 1;}
    tmat[0] = i, tmat[1] = 0, tmat[2] = 0;
    tmat[3] = 0, tmat[4] = j, tmat[5] = 0;
    tmat[6] = 0, tmat[7] = 0, tmat[8] = k;
    reset();
    return 1;
  }
  else {return 0;}
}

static int step4(void)
{
  int i, j, k;
  if (l * m * n == 0 || l * m * n == -1) {
    if (l == -1) {i = -1;} else {i = 1;}
    if (m == -1) {j = -1;} else {j = 1;}
    if (n == -1) {k = -1;} else {k = 1;}

    if (i * j * k == -1) {
      if (l == 0) {i = -1;}
      if (m == 0) {j = -1;}
      if (n == 0) {k = -1;}
    }
    
    tmat[0] = i, tmat[1] = 0, tmat[2] = 0;
    tmat[3] = 0, tmat[4] = j, tmat[5] = 0;
    tmat[6] = 0, tmat[7] = 0, tmat[8] = k;
    reset();
    return 1;
  }
  else {return 0;}
}

static int step5(void)
{
  if (fabs(xi) > B + eps ||
      (! (fabs(B - xi) > eps) && 2 * eta < zeta - eps) ||
      (! (fabs(B + xi) > eps) && zeta < -eps)) {
    tmat[0] = 1, tmat[1] = 0, tmat[2] = 0;
    tmat[3] = 0, tmat[4] = 1, tmat[5] = 0;
    tmat[6] = 0, tmat[7] = 0, tmat[8] = 1;
    if (xi > 0) {tmat[5] = -1;}
    if (xi < 0) {tmat[5] = 1;}
    reset();
    return 1;
  }
  else {return 0;}
}

static int step6(void)
{
  if (fabs(eta) > A + eps ||
      (! (fabs(A - eta) > eps) && 2 * xi < zeta - eps) ||
      (! (fabs(A + eta) > eps) && zeta < -eps)) {
    tmat[0] = 1, tmat[1] = 0, tmat[2] = 0;
    tmat[3] = 0, tmat[4] = 1, tmat[5] = 0;
    tmat[6] = 0, tmat[7] = 0, tmat[8] = 1;
    if (eta > 0) {tmat[2] = -1;}
    if (eta < 0) {tmat[2] = 1;}
    reset();
    return 1;
  }
  else {return 0;}
}

static int step7(void)
{
  if (fabs(zeta) > A + eps ||
      (! (fabs(A - zeta) > eps) && 2 * xi < eta - eps) ||
      (! (fabs(A + zeta) > eps) && eta < -eps)) {
    tmat[0] = 1, tmat[1] = 0, tmat[2] = 0;
    tmat[3] = 0, tmat[4] = 1, tmat[5] = 0;
    tmat[6] = 0, tmat[7] = 0, tmat[8] = 1;
    if (zeta > 0) {tmat[1] = -1;}
    if (zeta < 0) {tmat[1] = 1;}
    reset();
    return 1;
  }
  else {return 0;}
}

static int step8(void)
{
  if (xi + eta + zeta + A + B < -eps ||
      (! (fabs(xi + eta + zeta + A + B) > eps) && 2 * (A + eta) + zeta > eps)) {
    tmat[0] = 1, tmat[1] = 0, tmat[2] = 1;
    tmat[3] = 0, tmat[4] = 1, tmat[5] = 1;
    tmat[6] = 0, tmat[7] = 0, tmat[8] = 1;
    reset();
    return 1;
  }
  else {return 0;}
}

static void set_angle_types(void)
{
  l = 0, m = 0, n = 0;
  if (xi < -eps) {l = -1;}
  if (xi > eps) {l = 1;}
  if (eta < -eps) {m = -1;}
  if (eta > eps) {m = 1;}
  if (zeta < -eps) {n = -1;}
  if (zeta > eps) {n = 1;}
}

static void set_parameters(void)
{
  double *G;

  G = NULL;

  G = get_metric(lattice);

  A = G[0];
  B = G[4];
  C = G[8];
  xi = G[5] * 2;
  eta = G[2] * 2;
  zeta = G[1] * 2;

  free(G);
  G = NULL;
}

static double * get_transpose(const double *M)
{
  int i, j;
  double *M_T;

  M_T = NULL;

  if ((M_T = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    return NULL;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M_T[i * 3 + j] = M[j * 3 + i];
    }
  }
  return M_T;
}

static double * get_metric(const double *M)
{
  double *G, *M_T;

  G = NULL;
  M_T = NULL;

  M_T = get_transpose(M);

  G = multiply_matrices(M_T, M);
  free(M_T);
  M_T = NULL;

  return G;
}

static double * multiply_matrices(const double *L, const double *R)
{
  int i, j, k;
  double *M;

  M = NULL;

  if ((M = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    return NULL;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M[i * 3 + j] = 0;
      for (k = 0; k < 3; k++) {
	M[i * 3 + j] += L[i * 3 + k] * R[k * 3 + j];
      }
    }
  }
  return M;
}
