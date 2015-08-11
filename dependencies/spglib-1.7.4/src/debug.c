/* debug.c */
/* Copyright (C) 2008 Atsushi Togo */

#ifdef SPGDEBUG
#include <stdio.h>
#include "debug.h"

void dbg_print_matrix_d3(double a[3][3])
{
    int i;
    for (i = 0; i < 3; i++) {
        printf("%f %f %f\n", a[i][0], a[i][1], a[i][2]);
    }
}

void dbg_print_matrix_i3(int a[3][3])
{
    int i;
    for (i = 0; i < 3; i++) {
        printf("%d %d %d\n", a[i][0], a[i][1], a[i][2]);
    }
}

void dbg_print_vectors_d3(double a[][3], int size)
{
    int i;
    for (i = 0; i < size; i++) {
        printf("%d: %f %f %f\n", i + 1, a[i][0], a[i][1], a[i][2]);
    }
}

void dbg_print_vectors_with_label(double a[][3], int b[], int size)
{
    int i;
    for (i = 0; i < size; i++) {
        printf("%d: %f %f %f\n", b[i], a[i][0], a[i][1], a[i][2]);
    }
}

#endif
