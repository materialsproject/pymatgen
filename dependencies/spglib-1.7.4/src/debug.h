/* debug.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __debug_H__
#define __debug_H__

#ifdef SPGDEBUG
#define debug_print(...) printf(__VA_ARGS__)
#define debug_print_matrix_d3( a ) dbg_print_matrix_d3( a )
#define debug_print_matrix_i3( a ) dbg_print_matrix_i3( a )
#define debug_print_vectors_d3(...) dbg_print_vectors_d3(__VA_ARGS__)
#define debug_print_vectors_with_label(...) dbg_print_vectors_with_label(__VA_ARGS__)

void dbg_print_matrix_d3(double a[3][3]);
void dbg_print_matrix_i3(int a[3][3]);
void dbg_print_vectors_d3(double a[][3], int size);
void dbg_print_vectors_with_label(double a[][3], int b[], int size);

#else
#define debug_print(...)
#define debug_print_matrix_d3( a )
#define debug_print_matrix_i3( a )
#define debug_print_vectors_d3(...)
#define debug_print_vectors_with_label(...)
#endif

#ifdef SPGWARNING
#include <stdio.h>
#define warning_print(...) fprintf(stderr,__VA_ARGS__)
#else
#define warning_print(...)
#endif

#endif
