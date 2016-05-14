/* Copyright (C) 2015 Atsushi Togo */
/* All rights reserved. */

/* This file was originally part of spglib and is part of kspclib. */

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

#ifndef __kgrid_H__
#define __kgrid_H__

/* #define GRID_ORDER_XYZ */
/* This changes behaviour of index order of address. */
/* Without GRID_ORDER_XYZ, left most element of address runs first. */
/* grid_address (e.g. 4x4x4 mesh, unless GRID_ORDER_XYZ is defined) */
/*    [[ 0  0  0]                                                   */
/*     [ 1  0  0]                                                   */
/*     [ 2  0  0]                                                   */
/*     [-1  0  0]                                                   */
/*     [ 0  1  0]                                                   */
/*     [ 1  1  0]                                                   */
/*     [ 2  1  0]                                                   */
/*     [-1  1  0]                                                   */
/*     ....      ]                                                  */
/*                                                                  */
/* With GRID_ORDER_XYZ, right most element of address runs first.   */
/* grid_address (e.g. 4x4x4 mesh, if GRID_ORDER_XYZ is defined)     */
/*    [[ 0  0  0]                                                   */
/*     [ 0  0  1]                                                   */
/*     [ 0  0  2]                                                   */
/*     [ 0  0 -1]                                                   */
/*     [ 0  1  0]                                                   */
/*     [ 0  1  1]                                                   */
/*     [ 0  1  2]                                                   */
/*     [ 0  1 -1]                                                   */
/*     ....      ]                                                  */

/* #define GRID_BOUNDARY_AS_NEGATIVE */
/* This changes the behaviour of address elements on the surface of  */
/* parallelepiped. */
/* For odd mesh number, this affects nothing, e.g., [-2, -1, 0, 1, 2]. */
/* regardless of with and without GRID_BOUNDARY_AS_NEGATIVE. */
/* For even mesh number, this affects as follows: */
/* without GRID_BOUNDARY_AS_NEGATIVE, e.g., [-2, -1, 0, 1, 2, 3]. */
/* with GRID_BOUNDARY_AS_NEGATIVE, e.g., [-3, -2, -1, 0, 1, 2]. */

void kgd_get_all_grid_addresses(int grid_address[][3], const int mesh[3]);
int kgd_get_grid_point_double_mesh(const int address_double[3],
				   const int mesh[3]);
void kgd_get_grid_address_double_mesh(int address_double[3],
				      const int address[3],
				      const int mesh[3],
				      const int is_shift[3]);

#endif
