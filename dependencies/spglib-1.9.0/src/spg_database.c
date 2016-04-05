/* Copyright (C) 2010 Atsushi Togo */
/* All rights reserved. */

/* This file is part of spglib. */

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

#include <stdlib.h>
#include "spg_database.h"

/* In Hall symbols (3rd column), '=' is used instead of '"'. */
static const SpacegroupType spacegroup_types[] = {
  {  0, "      ", "                ", "                               ", "                   ", "          ", "     ", CENTERING_ERROR,  0 }, /*   0 */
  {  1, "C1^1  ", "P 1             ", "P 1                            ", "P 1                ", "P1        ", "     ",       PRIMITIVE,  1 }, /*   1 */
  {  2, "Ci^1  ", "-P 1            ", "P -1                           ", "P -1               ", "P-1       ", "     ",       PRIMITIVE,  2 }, /*   2 */
  {  3, "C2^1  ", "P 2y            ", "P 2 = P 1 2 1                  ", "P 1 2 1            ", "P2        ", "b    ",       PRIMITIVE,  3 }, /*   3 */
  {  3, "C2^1  ", "P 2             ", "P 2 = P 1 1 2                  ", "P 1 1 2            ", "P2        ", "c    ",       PRIMITIVE,  3 }, /*   4 */
  {  3, "C2^1  ", "P 2x            ", "P 2 = P 2 1 1                  ", "P 2 1 1            ", "P2        ", "a    ",       PRIMITIVE,  3 }, /*   5 */
  {  4, "C2^2  ", "P 2yb           ", "P 2_1 = P 1 2_1 1              ", "P 1 2_1 1          ", "P2_1      ", "b    ",       PRIMITIVE,  3 }, /*   6 */
  {  4, "C2^2  ", "P 2c            ", "P 2_1 = P 1 1 2_1              ", "P 1 1 2_1          ", "P2_1      ", "c    ",       PRIMITIVE,  3 }, /*   7 */
  {  4, "C2^2  ", "P 2xa           ", "P 2_1 = P 2_1 1 1              ", "P 2_1 1 1          ", "P2_1      ", "a    ",       PRIMITIVE,  3 }, /*   8 */
  {  5, "C2^3  ", "C 2y            ", "C 2 = C 1 2 1                  ", "C 1 2 1            ", "C2        ", "b1   ",          C_FACE,  3 }, /*   9 */
  {  5, "C2^3  ", "A 2y            ", "C 2 = A 1 2 1                  ", "A 1 2 1            ", "C2        ", "b2   ",          A_FACE,  3 }, /*  10 */
  {  5, "C2^3  ", "I 2y            ", "C 2 = I 1 2 1                  ", "I 1 2 1            ", "C2        ", "b3   ",            BODY,  3 }, /*  11 */
  {  5, "C2^3  ", "A 2             ", "C 2 = A 1 1 2                  ", "A 1 1 2            ", "C2        ", "c1   ",          A_FACE,  3 }, /*  12 */
  {  5, "C2^3  ", "B 2             ", "C 2 = B 1 1 2 = B 2            ", "B 1 1 2            ", "C2        ", "c2   ",       PRIMITIVE,  3 }, /*  13 */
  {  5, "C2^3  ", "I 2             ", "C 2 = I 1 1 2                  ", "I 1 1 2            ", "C2        ", "c3   ",            BODY,  3 }, /*  14 */
  {  5, "C2^3  ", "B 2x            ", "C 2 = B 2 1 1                  ", "B 2 1 1            ", "C2        ", "a1   ",       PRIMITIVE,  3 }, /*  15 */
  {  5, "C2^3  ", "C 2x            ", "C 2 = C 2 1 1                  ", "C 2 1 1            ", "C2        ", "a2   ",          C_FACE,  3 }, /*  16 */
  {  5, "C2^3  ", "I 2x            ", "C 2 = I 2 1 1                  ", "I 2 1 1            ", "C2        ", "a3   ",            BODY,  3 }, /*  17 */
  {  6, "Cs^1  ", "P -2y           ", "P m = P 1 m 1                  ", "P 1 m 1            ", "Pm        ", "b    ",       PRIMITIVE,  4 }, /*  18 */
  {  6, "Cs^1  ", "P -2            ", "P m = P 1 1 m                  ", "P 1 1 m            ", "Pm        ", "c    ",       PRIMITIVE,  4 }, /*  19 */
  {  6, "Cs^1  ", "P -2x           ", "P m = P m 1 1                  ", "P m 1 1            ", "Pm        ", "a    ",       PRIMITIVE,  4 }, /*  20 */
  {  7, "Cs^2  ", "P -2yc          ", "P c = P 1 c 1                  ", "P 1 c 1            ", "Pc        ", "b1   ",       PRIMITIVE,  4 }, /*  21 */
  {  7, "Cs^2  ", "P -2yac         ", "P c = P 1 n 1                  ", "P 1 n 1            ", "Pc        ", "b2   ",       PRIMITIVE,  4 }, /*  22 */
  {  7, "Cs^2  ", "P -2ya          ", "P c = P 1 a 1                  ", "P 1 a 1            ", "Pc        ", "b3   ",       PRIMITIVE,  4 }, /*  23 */
  {  7, "Cs^2  ", "P -2a           ", "P c = P 1 1 a                  ", "P 1 1 a            ", "Pc        ", "c1   ",       PRIMITIVE,  4 }, /*  24 */
  {  7, "Cs^2  ", "P -2ab          ", "P c = P 1 1 n                  ", "P 1 1 n            ", "Pc        ", "c2   ",       PRIMITIVE,  4 }, /*  25 */
  {  7, "Cs^2  ", "P -2b           ", "P c = P 1 1 b = P b            ", "P 1 1 b            ", "Pc        ", "c3   ",       PRIMITIVE,  4 }, /*  26 */
  {  7, "Cs^2  ", "P -2xb          ", "P c = P b 1 1                  ", "P b 1 1            ", "Pc        ", "a1   ",       PRIMITIVE,  4 }, /*  27 */
  {  7, "Cs^2  ", "P -2xbc         ", "P c = P n 1 1                  ", "P n 1 1            ", "Pc        ", "a2   ",       PRIMITIVE,  4 }, /*  28 */
  {  7, "Cs^2  ", "P -2xc          ", "P c = P c 1 1                  ", "P c 1 1            ", "Pc        ", "a3   ",       PRIMITIVE,  4 }, /*  29 */
  {  8, "Cs^3  ", "C -2y           ", "C m = C 1 m 1                  ", "C 1 m 1            ", "Cm        ", "b1   ",          C_FACE,  4 }, /*  30 */
  {  8, "Cs^3  ", "A -2y           ", "C m = A 1 m 1                  ", "A 1 m 1            ", "Cm        ", "b2   ",          A_FACE,  4 }, /*  31 */
  {  8, "Cs^3  ", "I -2y           ", "C m = I 1 m 1                  ", "I 1 m 1            ", "Cm        ", "b3   ",            BODY,  4 }, /*  32 */
  {  8, "Cs^3  ", "A -2            ", "C m = A 1 1 m                  ", "A 1 1 m            ", "Cm        ", "c1   ",          A_FACE,  4 }, /*  33 */
  {  8, "Cs^3  ", "B -2            ", "C m = B 1 1 m = B m            ", "B 1 1 m            ", "Cm        ", "c2   ",       PRIMITIVE,  4 }, /*  34 */
  {  8, "Cs^3  ", "I -2            ", "C m = I 1 1 m                  ", "I 1 1 m            ", "Cm        ", "c3   ",            BODY,  4 }, /*  35 */
  {  8, "Cs^3  ", "B -2x           ", "C m = B m 1 1                  ", "B m 1 1            ", "Cm        ", "a1   ",       PRIMITIVE,  4 }, /*  36 */
  {  8, "Cs^3  ", "C -2x           ", "C m = C m 1 1                  ", "C m 1 1            ", "Cm        ", "a2   ",          C_FACE,  4 }, /*  37 */
  {  8, "Cs^3  ", "I -2x           ", "C m = I m 1 1                  ", "I m 1 1            ", "Cm        ", "a3   ",            BODY,  4 }, /*  38 */
  {  9, "Cs^4  ", "C -2yc          ", "C c = C 1 c 1                  ", "C 1 c 1            ", "Cc        ", "b1   ",          C_FACE,  4 }, /*  39 */
  {  9, "Cs^4  ", "A -2yac         ", "C c = A 1 n 1                  ", "A 1 n 1            ", "Cc        ", "b2   ",          A_FACE,  4 }, /*  40 */
  {  9, "Cs^4  ", "I -2ya          ", "C c = I 1 a 1                  ", "I 1 a 1            ", "Cc        ", "b3   ",            BODY,  4 }, /*  41 */
  {  9, "Cs^4  ", "A -2ya          ", "C c = A 1 a 1                  ", "A 1 a 1            ", "Cc        ", "-b1  ",          A_FACE,  4 }, /*  42 */
  {  9, "Cs^4  ", "C -2ybc         ", "C c = C 1 n 1                  ", "C 1 n 1            ", "Cc        ", "-b2  ",          C_FACE,  4 }, /*  43 */
  {  9, "Cs^4  ", "I -2yc          ", "C c = I 1 c 1                  ", "I 1 c 1            ", "Cc        ", "-b3  ",            BODY,  4 }, /*  44 */
  {  9, "Cs^4  ", "A -2a           ", "C c = A 1 1 a                  ", "A 1 1 a            ", "Cc        ", "c1   ",          A_FACE,  4 }, /*  45 */
  {  9, "Cs^4  ", "B -2bc          ", "C c = B 1 1 n                  ", "B 1 1 n            ", "Cc        ", "c2   ",       PRIMITIVE,  4 }, /*  46 */
  {  9, "Cs^4  ", "I -2b           ", "C c = I 1 1 b                  ", "I 1 1 b            ", "Cc        ", "c3   ",            BODY,  4 }, /*  47 */
  {  9, "Cs^4  ", "B -2b           ", "C c = B 1 1 b = B b            ", "B 1 1 b            ", "Cc        ", "-c1  ",       PRIMITIVE,  4 }, /*  48 */
  {  9, "Cs^4  ", "A -2ac          ", "C c = A 1 1 n                  ", "A 1 1 n            ", "Cc        ", "-c2  ",          A_FACE,  4 }, /*  49 */
  {  9, "Cs^4  ", "I -2a           ", "C c = I 1 1 a                  ", "I 1 1 a            ", "Cc        ", "-c3  ",            BODY,  4 }, /*  50 */
  {  9, "Cs^4  ", "B -2xb          ", "C c = B b 1 1                  ", "B b 1 1            ", "Cc        ", "a1   ",       PRIMITIVE,  4 }, /*  51 */
  {  9, "Cs^4  ", "C -2xbc         ", "C c = C n 1 1                  ", "C n 1 1            ", "Cc        ", "a2   ",          C_FACE,  4 }, /*  52 */
  {  9, "Cs^4  ", "I -2xc          ", "C c = I c 1 1                  ", "I c 1 1            ", "Cc        ", "a3   ",            BODY,  4 }, /*  53 */
  {  9, "Cs^4  ", "C -2xc          ", "C c = C c 1 1                  ", "C c 1 1            ", "Cc        ", "-a1  ",          C_FACE,  4 }, /*  54 */
  {  9, "Cs^4  ", "B -2xbc         ", "C c = B n 1 1                  ", "B n 1 1            ", "Cc        ", "-a2  ",       PRIMITIVE,  4 }, /*  55 */
  {  9, "Cs^4  ", "I -2xb          ", "C c = I b 1 1                  ", "I b 1 1            ", "Cc        ", "-a3  ",            BODY,  4 }, /*  56 */
  { 10, "C2h^1 ", "-P 2y           ", "P 2/m = P 1 2/m 1              ", "P 1 2/m 1          ", "P2/m      ", "b    ",       PRIMITIVE,  5 }, /*  57 */
  { 10, "C2h^1 ", "-P 2            ", "P 2/m = P 1 1 2/m              ", "P 1 1 2/m          ", "P2/m      ", "c    ",       PRIMITIVE,  5 }, /*  58 */
  { 10, "C2h^1 ", "-P 2x           ", "P 2/m = P 2/m 1 1              ", "P 2/m 1 1          ", "P2/m      ", "a    ",       PRIMITIVE,  5 }, /*  59 */
  { 11, "C2h^2 ", "-P 2yb          ", "P 2_1/m = P 1 2_1/m 1          ", "P 1 2_1/m 1        ", "P2_1/m    ", "b    ",       PRIMITIVE,  5 }, /*  60 */
  { 11, "C2h^2 ", "-P 2c           ", "P 2_1/m = P 1 1 2_1/m          ", "P 1 1 2_1/m        ", "P2_1/m    ", "c    ",       PRIMITIVE,  5 }, /*  61 */
  { 11, "C2h^2 ", "-P 2xa          ", "P 2_1/m = P 2_1/m 1 1          ", "P 2_1/m 1 1        ", "P2_1/m    ", "a    ",       PRIMITIVE,  5 }, /*  62 */
  { 12, "C2h^3 ", "-C 2y           ", "C 2/m = C 1 2/m 1              ", "C 1 2/m 1          ", "C2/m      ", "b1   ",          C_FACE,  5 }, /*  63 */
  { 12, "C2h^3 ", "-A 2y           ", "C 2/m = A 1 2/m 1              ", "A 1 2/m 1          ", "C2/m      ", "b2   ",          A_FACE,  5 }, /*  64 */
  { 12, "C2h^3 ", "-I 2y           ", "C 2/m = I 1 2/m 1              ", "I 1 2/m 1          ", "C2/m      ", "b3   ",            BODY,  5 }, /*  65 */
  { 12, "C2h^3 ", "-A 2            ", "C 2/m = A 1 1 2/m              ", "A 1 1 2/m          ", "C2/m      ", "c1   ",          A_FACE,  5 }, /*  66 */
  { 12, "C2h^3 ", "-B 2            ", "C 2/m = B 1 1 2/m = B 2/m      ", "B 1 1 2/m          ", "C2/m      ", "c2   ",       PRIMITIVE,  5 }, /*  67 */
  { 12, "C2h^3 ", "-I 2            ", "C 2/m = I 1 1 2/m              ", "I 1 1 2/m          ", "C2/m      ", "c3   ",            BODY,  5 }, /*  68 */
  { 12, "C2h^3 ", "-B 2x           ", "C 2/m = B 2/m 1 1              ", "B 2/m 1 1          ", "C2/m      ", "a1   ",       PRIMITIVE,  5 }, /*  69 */
  { 12, "C2h^3 ", "-C 2x           ", "C 2/m = C 2/m 1 1              ", "C 2/m 1 1          ", "C2/m      ", "a2   ",          C_FACE,  5 }, /*  70 */
  { 12, "C2h^3 ", "-I 2x           ", "C 2/m = I 2/m 1 1              ", "I 2/m 1 1          ", "C2/m      ", "a3   ",            BODY,  5 }, /*  71 */
  { 13, "C2h^4 ", "-P 2yc          ", "P 2/c = P 1 2/c 1              ", "P 1 2/c 1          ", "P2/c      ", "b1   ",       PRIMITIVE,  5 }, /*  72 */
  { 13, "C2h^4 ", "-P 2yac         ", "P 2/c = P 1 2/n 1              ", "P 1 2/n 1          ", "P2/c      ", "b2   ",       PRIMITIVE,  5 }, /*  73 */
  { 13, "C2h^4 ", "-P 2ya          ", "P 2/c = P 1 2/a 1              ", "P 1 2/a 1          ", "P2/c      ", "b3   ",       PRIMITIVE,  5 }, /*  74 */
  { 13, "C2h^4 ", "-P 2a           ", "P 2/c = P 1 1 2/a              ", "P 1 1 2/a          ", "P2/c      ", "c1   ",       PRIMITIVE,  5 }, /*  75 */
  { 13, "C2h^4 ", "-P 2ab          ", "P 2/c = P 1 1 2/n              ", "P 1 1 2/n          ", "P2/c      ", "c2   ",       PRIMITIVE,  5 }, /*  76 */
  { 13, "C2h^4 ", "-P 2b           ", "P 2/c = P 1 1 2/b = P 2/b      ", "P 1 1 2/b          ", "P2/c      ", "c3   ",       PRIMITIVE,  5 }, /*  77 */
  { 13, "C2h^4 ", "-P 2xb          ", "P 2/c = P 2/b 1 1              ", "P 2/b 1 1          ", "P2/c      ", "a1   ",       PRIMITIVE,  5 }, /*  78 */
  { 13, "C2h^4 ", "-P 2xbc         ", "P 2/c = P 2/n 1 1              ", "P 2/n 1 1          ", "P2/c      ", "a2   ",       PRIMITIVE,  5 }, /*  79 */
  { 13, "C2h^4 ", "-P 2xc          ", "P 2/c = P 2/c 1 1              ", "P 2/c 1 1          ", "P2/c      ", "a3   ",       PRIMITIVE,  5 }, /*  80 */
  { 14, "C2h^5 ", "-P 2ybc         ", "P 2_1/c = P 1 2_1/c 1          ", "P 1 2_1/c 1        ", "P2_1/c    ", "b1   ",       PRIMITIVE,  5 }, /*  81 */
  { 14, "C2h^5 ", "-P 2yn          ", "P 2_1/c = P 1 2_1/n 1          ", "P 1 2_1/n 1        ", "P2_1/c    ", "b2   ",       PRIMITIVE,  5 }, /*  82 */
  { 14, "C2h^5 ", "-P 2yab         ", "P 2_1/c = P 1 2_1/a 1          ", "P 1 2_1/a 1        ", "P2_1/c    ", "b3   ",       PRIMITIVE,  5 }, /*  83 */
  { 14, "C2h^5 ", "-P 2ac          ", "P 2_1/c = P 1 1 2_1/a          ", "P 1 1 2_1/a        ", "P2_1/c    ", "c1   ",       PRIMITIVE,  5 }, /*  84 */
  { 14, "C2h^5 ", "-P 2n           ", "P 2_1/c = P 1 1 2_1/n          ", "P 1 1 2_1/n        ", "P2_1/c    ", "c2   ",       PRIMITIVE,  5 }, /*  85 */
  { 14, "C2h^5 ", "-P 2bc          ", "P 2_1/c = P 1 1 2_1/b = P 2_1/b", "P 1 1 2_1/b        ", "P2_1/c    ", "c3   ",       PRIMITIVE,  5 }, /*  86 */
  { 14, "C2h^5 ", "-P 2xab         ", "P 2_1/c = P 2_1/b 1 1          ", "P 2_1/b 1 1        ", "P2_1/c    ", "a1   ",       PRIMITIVE,  5 }, /*  87 */
  { 14, "C2h^5 ", "-P 2xn          ", "P 2_1/c = P 2_1/n 1 1          ", "P 2_1/n 1 1        ", "P2_1/c    ", "a2   ",       PRIMITIVE,  5 }, /*  88 */
  { 14, "C2h^5 ", "-P 2xac         ", "P 2_1/c = P 2_1/c 1 1          ", "P 2_1/c 1 1        ", "P2_1/c    ", "a3   ",       PRIMITIVE,  5 }, /*  89 */
  { 15, "C2h^6 ", "-C 2yc          ", "C 2/c = C 1 2/c 1              ", "C 1 2/c 1          ", "C2/c      ", "b1   ",          C_FACE,  5 }, /*  90 */
  { 15, "C2h^6 ", "-A 2yac         ", "C 2/c = A 1 2/n 1              ", "A 1 2/n 1          ", "C2/c      ", "b2   ",          A_FACE,  5 }, /*  91 */
  { 15, "C2h^6 ", "-I 2ya          ", "C 2/c = I 1 2/a 1              ", "I 1 2/a 1          ", "C2/c      ", "b3   ",            BODY,  5 }, /*  92 */
  { 15, "C2h^6 ", "-A 2ya          ", "C 2/c = A 1 2/a 1              ", "A 1 2/a 1          ", "C2/c      ", "-b1  ",          A_FACE,  5 }, /*  93 */
  { 15, "C2h^6 ", "-C 2ybc         ", "C 2/c = C 1 2/n 1              ", "C 1 2/n 1          ", "C2/c      ", "-b2  ",          C_FACE,  5 }, /*  94 */
  { 15, "C2h^6 ", "-I 2yc          ", "C 2/c = I 1 2/c 1              ", "I 1 2/c 1          ", "C2/c      ", "-b3  ",            BODY,  5 }, /*  95 */
  { 15, "C2h^6 ", "-A 2a           ", "C 2/c = A 1 1 2/a              ", "A 1 1 2/a          ", "C2/c      ", "c1   ",          A_FACE,  5 }, /*  96 */
  { 15, "C2h^6 ", "-B 2bc          ", "C 2/c = B 1 1 2/n              ", "B 1 1 2/n          ", "C2/c      ", "c2   ",       PRIMITIVE,  5 }, /*  97 */
  { 15, "C2h^6 ", "-I 2b           ", "C 2/c = I 1 1 2/b              ", "I 1 1 2/b          ", "C2/c      ", "c3   ",            BODY,  5 }, /*  98 */
  { 15, "C2h^6 ", "-B 2b           ", "C 2/c = B 1 1 2/b = B 2/b      ", "B 1 1 2/b          ", "C2/c      ", "-c1  ",       PRIMITIVE,  5 }, /*  99 */
  { 15, "C2h^6 ", "-A 2ac          ", "C 2/c = A 1 1 2/n              ", "A 1 1 2/n          ", "C2/c      ", "-c2  ",          A_FACE,  5 }, /* 100 */
  { 15, "C2h^6 ", "-I 2a           ", "C 2/c = I 1 1 2/a              ", "I 1 1 2/a          ", "C2/c      ", "-c3  ",            BODY,  5 }, /* 101 */
  { 15, "C2h^6 ", "-B 2xb          ", "C 2/c = B 2/b 1 1              ", "B 2/b 1 1          ", "C2/c      ", "a1   ",       PRIMITIVE,  5 }, /* 102 */
  { 15, "C2h^6 ", "-C 2xbc         ", "C 2/c = C 2/n 1 1              ", "C 2/n 1 1          ", "C2/c      ", "a2   ",          C_FACE,  5 }, /* 103 */
  { 15, "C2h^6 ", "-I 2xc          ", "C 2/c = I 2/c 1 1              ", "I 2/c 1 1          ", "C2/c      ", "a3   ",            BODY,  5 }, /* 104 */
  { 15, "C2h^6 ", "-C 2xc          ", "C 2/c = C 2/c 1 1              ", "C 2/c 1 1          ", "C2/c      ", "-a1  ",          C_FACE,  5 }, /* 105 */
  { 15, "C2h^6 ", "-B 2xbc         ", "C 2/c = B 2/n 1 1              ", "B 2/n 1 1          ", "C2/c      ", "-a2  ",       PRIMITIVE,  5 }, /* 106 */
  { 15, "C2h^6 ", "-I 2xb          ", "C 2/c = I 2/b 1 1              ", "I 2/b 1 1          ", "C2/c      ", "-a3  ",            BODY,  5 }, /* 107 */
  { 16, "D2^1  ", "P 2 2           ", "P 2 2 2                        ", "P 2 2 2            ", "P222      ", "     ",       PRIMITIVE,  6 }, /* 108 */
  { 17, "D2^2  ", "P 2c 2          ", "P 2 2 2_1                      ", "P 2 2 2_1          ", "P222_1    ", "     ",       PRIMITIVE,  6 }, /* 109 */
  { 17, "D2^2  ", "P 2a 2a         ", "P 2_1 2 2                      ", "P 2_1 2 2          ", "P2_122    ", "cab  ",       PRIMITIVE,  6 }, /* 110 */
  { 17, "D2^2  ", "P 2 2b          ", "P 2 2_1 2                      ", "P 2 2_1 2          ", "P22_12    ", "bca  ",       PRIMITIVE,  6 }, /* 111 */
  { 18, "D2^3  ", "P 2 2ab         ", "P 2_1 2_1 2                    ", "P 2_1 2_1 2        ", "P2_12_12  ", "     ",       PRIMITIVE,  6 }, /* 112 */
  { 18, "D2^3  ", "P 2bc 2         ", "P 2 2_1 2_1                    ", "P 2 2_1 2_1        ", "P22_12_1  ", "cab  ",       PRIMITIVE,  6 }, /* 113 */
  { 18, "D2^3  ", "P 2ac 2ac       ", "P 2_1 2 2_1                    ", "P 2_1 2 2_1        ", "P2_122_1  ", "bca  ",       PRIMITIVE,  6 }, /* 114 */
  { 19, "D2^4  ", "P 2ac 2ab       ", "P 2_1 2_1 2_1                  ", "P 2_1 2_1 2_1      ", "P2_12_12_1", "     ",       PRIMITIVE,  6 }, /* 115 */
  { 20, "D2^5  ", "C 2c 2          ", "C 2 2 2_1                      ", "C 2 2 2_1          ", "C222_1    ", "     ",          C_FACE,  6 }, /* 116 */
  { 20, "D2^5  ", "A 2a 2a         ", "A 2_1 2 2                      ", "A 2_1 2 2          ", "A2_122    ", "cab  ",          A_FACE,  6 }, /* 117 */
  { 20, "D2^5  ", "B 2 2b          ", "B 2 2_1 2                      ", "B 2 2_1 2          ", "B22_12    ", "bca  ",       PRIMITIVE,  6 }, /* 118 */
  { 21, "D2^6  ", "C 2 2           ", "C 2 2 2                        ", "C 2 2 2            ", "C222      ", "     ",          C_FACE,  6 }, /* 119 */
  { 21, "D2^6  ", "A 2 2           ", "A 2 2 2                        ", "A 2 2 2            ", "A222      ", "cab  ",          A_FACE,  6 }, /* 120 */
  { 21, "D2^6  ", "B 2 2           ", "B 2 2 2                        ", "B 2 2 2            ", "B222      ", "bca  ",       PRIMITIVE,  6 }, /* 121 */
  { 22, "D2^7  ", "F 2 2           ", "F 2 2 2                        ", "F 2 2 2            ", "F222      ", "     ",            FACE,  6 }, /* 122 */
  { 23, "D2^8  ", "I 2 2           ", "I 2 2 2                        ", "I 2 2 2            ", "I222      ", "     ",            BODY,  6 }, /* 123 */
  { 24, "D2^9  ", "I 2b 2c         ", "I 2_1 2_1 2_1                  ", "I 2_1 2_1 2_1      ", "I2_12_12_1", "     ",            BODY,  6 }, /* 124 */
  { 25, "C2v^1 ", "P 2 -2          ", "P m m 2                        ", "P m m 2            ", "Pmm2      ", "     ",       PRIMITIVE,  7 }, /* 125 */
  { 25, "C2v^1 ", "P -2 2          ", "P 2 m m                        ", "P 2 m m            ", "P2mm      ", "cab  ",       PRIMITIVE,  7 }, /* 126 */
  { 25, "C2v^1 ", "P -2 -2         ", "P m 2 m                        ", "P m 2 m            ", "Pm2m      ", "bca  ",       PRIMITIVE,  7 }, /* 127 */
  { 26, "C2v^2 ", "P 2c -2         ", "P m c 2_1                      ", "P m c 2_1          ", "Pmc2_1    ", "     ",       PRIMITIVE,  7 }, /* 128 */
  { 26, "C2v^2 ", "P 2c -2c        ", "P c m 2_1                      ", "P c m 2_1          ", "Pcm2_1    ", "ba-c ",       PRIMITIVE,  7 }, /* 129 */
  { 26, "C2v^2 ", "P -2a 2a        ", "P 2_1 m a                      ", "P 2_1 m a          ", "P2_1ma    ", "cab  ",       PRIMITIVE,  7 }, /* 130 */
  { 26, "C2v^2 ", "P -2 2a         ", "P 2_1 a m                      ", "P 2_1 a m          ", "P2_1am    ", "-cba ",       PRIMITIVE,  7 }, /* 131 */
  { 26, "C2v^2 ", "P -2 -2b        ", "P b 2_1 m                      ", "P b 2_1 m          ", "Pb2_1m    ", "bca  ",       PRIMITIVE,  7 }, /* 132 */
  { 26, "C2v^2 ", "P -2b -2        ", "P m 2_1 b                      ", "P m 2_1 b          ", "Pm2_1b    ", "a-cb ",       PRIMITIVE,  7 }, /* 133 */
  { 27, "C2v^3 ", "P 2 -2c         ", "P c c 2                        ", "P c c 2            ", "Pcc2      ", "     ",       PRIMITIVE,  7 }, /* 134 */
  { 27, "C2v^3 ", "P -2a 2         ", "P 2 a a                        ", "P 2 a a            ", "P2aa      ", "cab  ",       PRIMITIVE,  7 }, /* 135 */
  { 27, "C2v^3 ", "P -2b -2b       ", "P b 2 b                        ", "P b 2 b            ", "Pb2b      ", "bca  ",       PRIMITIVE,  7 }, /* 136 */
  { 28, "C2v^4 ", "P 2 -2a         ", "P m a 2                        ", "P m a 2            ", "Pma2      ", "     ",       PRIMITIVE,  7 }, /* 137 */
  { 28, "C2v^4 ", "P 2 -2b         ", "P b m 2                        ", "P b m 2            ", "Pbm2      ", "ba-c ",       PRIMITIVE,  7 }, /* 138 */
  { 28, "C2v^4 ", "P -2b 2         ", "P 2 m b                        ", "P 2 m b            ", "P2mb      ", "cab  ",       PRIMITIVE,  7 }, /* 139 */
  { 28, "C2v^4 ", "P -2c 2         ", "P 2 c m                        ", "P 2 c m            ", "P2cm      ", "-cba ",       PRIMITIVE,  7 }, /* 140 */
  { 28, "C2v^4 ", "P -2c -2c       ", "P c 2 m                        ", "P c 2 m            ", "Pc2m      ", "bca  ",       PRIMITIVE,  7 }, /* 141 */
  { 28, "C2v^4 ", "P -2a -2a       ", "P m 2 a                        ", "P m 2 a            ", "Pm2a      ", "a-cb ",       PRIMITIVE,  7 }, /* 142 */
  { 29, "C2v^5 ", "P 2c -2ac       ", "P c a 2_1                      ", "P c a 2_1          ", "Pca2_1    ", "     ",       PRIMITIVE,  7 }, /* 143 */
  { 29, "C2v^5 ", "P 2c -2b        ", "P b c 2_1                      ", "P b c 2_1          ", "Pbc2_1    ", "ba-c ",       PRIMITIVE,  7 }, /* 144 */
  { 29, "C2v^5 ", "P -2b 2a        ", "P 2_1 a b                      ", "P 2_1 a b          ", "P2_1ab    ", "cab  ",       PRIMITIVE,  7 }, /* 145 */
  { 29, "C2v^5 ", "P -2ac 2a       ", "P 2_1 c a                      ", "P 2_1 c a          ", "P2_1ca    ", "-cba ",       PRIMITIVE,  7 }, /* 146 */
  { 29, "C2v^5 ", "P -2bc -2c      ", "P c 2_1 b                      ", "P c 2_1 b          ", "Pc2_1b    ", "bca  ",       PRIMITIVE,  7 }, /* 147 */
  { 29, "C2v^5 ", "P -2a -2ab      ", "P b 2_1 a                      ", "P b 2_1 a          ", "Pb2_1a    ", "a-cb ",       PRIMITIVE,  7 }, /* 148 */
  { 30, "C2v^6 ", "P 2 -2bc        ", "P n c 2                        ", "P n c 2            ", "Pnc2      ", "     ",       PRIMITIVE,  7 }, /* 149 */
  { 30, "C2v^6 ", "P 2 -2ac        ", "P c n 2                        ", "P c n 2            ", "Pcn2      ", "ba-c ",       PRIMITIVE,  7 }, /* 150 */
  { 30, "C2v^6 ", "P -2ac 2        ", "P 2 n a                        ", "P 2 n a            ", "P2na      ", "cab  ",       PRIMITIVE,  7 }, /* 151 */
  { 30, "C2v^6 ", "P -2ab 2        ", "P 2 a n                        ", "P 2 a n            ", "P2an      ", "-cba ",       PRIMITIVE,  7 }, /* 152 */
  { 30, "C2v^6 ", "P -2ab -2ab     ", "P b 2 n                        ", "P b 2 n            ", "Pb2n      ", "bca  ",       PRIMITIVE,  7 }, /* 153 */
  { 30, "C2v^6 ", "P -2bc -2bc     ", "P n 2 b                        ", "P n 2 b            ", "Pn2b      ", "a-cb ",       PRIMITIVE,  7 }, /* 154 */
  { 31, "C2v^7 ", "P 2ac -2        ", "P m n 2_1                      ", "P m n 2_1          ", "Pmn2_1    ", "     ",       PRIMITIVE,  7 }, /* 155 */
  { 31, "C2v^7 ", "P 2bc -2bc      ", "P n m 2_1                      ", "P n m 2_1          ", "Pnm2_1    ", "ba-c ",       PRIMITIVE,  7 }, /* 156 */
  { 31, "C2v^7 ", "P -2ab 2ab      ", "P 2_1 m n                      ", "P 2_1 m n          ", "P2_1mn    ", "cab  ",       PRIMITIVE,  7 }, /* 157 */
  { 31, "C2v^7 ", "P -2 2ac        ", "P 2_1 n m                      ", "P 2_1 n m          ", "P2_1nm    ", "-cba ",       PRIMITIVE,  7 }, /* 158 */
  { 31, "C2v^7 ", "P -2 -2bc       ", "P n 2_1 m                      ", "P n 2_1 m          ", "Pn2_1m    ", "bca  ",       PRIMITIVE,  7 }, /* 159 */
  { 31, "C2v^7 ", "P -2ab -2       ", "P m 2_1 n                      ", "P m 2_1 n          ", "Pm2_1n    ", "a-cb ",       PRIMITIVE,  7 }, /* 160 */
  { 32, "C2v^8 ", "P 2 -2ab        ", "P b a 2                        ", "P b a 2            ", "Pba2      ", "     ",       PRIMITIVE,  7 }, /* 161 */
  { 32, "C2v^8 ", "P -2bc 2        ", "P 2 c b                        ", "P 2 c b            ", "P2cb      ", "cab  ",       PRIMITIVE,  7 }, /* 162 */
  { 32, "C2v^8 ", "P -2ac -2ac     ", "P c 2 a                        ", "P c 2 a            ", "Pc2a      ", "bca  ",       PRIMITIVE,  7 }, /* 163 */
  { 33, "C2v^9 ", "P 2c -2n        ", "P n a 2_1                      ", "P n a 2_1          ", "Pna2_1    ", "     ",       PRIMITIVE,  7 }, /* 164 */
  { 33, "C2v^9 ", "P 2c -2ab       ", "P b n 2_1                      ", "P b n 2_1          ", "Pbn2_1    ", "ba-c ",       PRIMITIVE,  7 }, /* 165 */
  { 33, "C2v^9 ", "P -2bc 2a       ", "P 2_1 n b                      ", "P 2_1 n b          ", "P2_1nb    ", "cab  ",       PRIMITIVE,  7 }, /* 166 */
  { 33, "C2v^9 ", "P -2n 2a        ", "P 2_1 c n                      ", "P 2_1 c n          ", "P2_1cn    ", "-cba ",       PRIMITIVE,  7 }, /* 167 */
  { 33, "C2v^9 ", "P -2n -2ac      ", "P c 2_1 n                      ", "P c 2_1 n          ", "Pc2_1n    ", "bca  ",       PRIMITIVE,  7 }, /* 168 */
  { 33, "C2v^9 ", "P -2ac -2n      ", "P n 2_1 a                      ", "P n 2_1 a          ", "Pn2_1a    ", "a-cb ",       PRIMITIVE,  7 }, /* 169 */
  { 34, "C2v^10", "P 2 -2n         ", "P n n 2                        ", "P n n 2            ", "Pnn2      ", "     ",       PRIMITIVE,  7 }, /* 170 */
  { 34, "C2v^10", "P -2n 2         ", "P 2 n n                        ", "P 2 n n            ", "P2nn      ", "cab  ",       PRIMITIVE,  7 }, /* 171 */
  { 34, "C2v^10", "P -2n -2n       ", "P n 2 n                        ", "P n 2 n            ", "Pn2n      ", "bca  ",       PRIMITIVE,  7 }, /* 172 */
  { 35, "C2v^11", "C 2 -2          ", "C m m 2                        ", "C m m 2            ", "Cmm2      ", "     ",          C_FACE,  7 }, /* 173 */
  { 35, "C2v^11", "A -2 2          ", "A 2 m m                        ", "A 2 m m            ", "A2mm      ", "cab  ",          A_FACE,  7 }, /* 174 */
  { 35, "C2v^11", "B -2 -2         ", "B m 2 m                        ", "B m 2 m            ", "Bm2m      ", "bca  ",       PRIMITIVE,  7 }, /* 175 */
  { 36, "C2v^12", "C 2c -2         ", "C m c 2_1                      ", "C m c 2_1          ", "Cmc2_1    ", "     ",          C_FACE,  7 }, /* 176 */
  { 36, "C2v^12", "C 2c -2c        ", "C c m 2_1                      ", "C c m 2_1          ", "Ccm2_1    ", "ba-c ",          C_FACE,  7 }, /* 177 */
  { 36, "C2v^12", "A -2a 2a        ", "A 2_1 m a                      ", "A 2_1 m a          ", "A2_1ma    ", "cab  ",          A_FACE,  7 }, /* 178 */
  { 36, "C2v^12", "A -2 2a         ", "A 2_1 a m                      ", "A 2_1 a m          ", "A2_1am    ", "-cba ",          A_FACE,  7 }, /* 179 */
  { 36, "C2v^12", "B -2 -2b        ", "B b 2_1 m                      ", "B b 2_1 m          ", "Bb2_1m    ", "bca  ",       PRIMITIVE,  7 }, /* 180 */
  { 36, "C2v^12", "B -2b -2        ", "B m 2_1 b                      ", "B m 2_1 b          ", "Bm2_1b    ", "a-cb ",       PRIMITIVE,  7 }, /* 181 */
  { 37, "C2v^13", "C 2 -2c         ", "C c c 2                        ", "C c c 2            ", "Ccc2      ", "     ",          C_FACE,  7 }, /* 182 */
  { 37, "C2v^13", "A -2a 2         ", "A 2 a a                        ", "A 2 a a            ", "A2aa      ", "cab  ",          A_FACE,  7 }, /* 183 */
  { 37, "C2v^13", "B -2b -2b       ", "B b 2 b                        ", "B b 2 b            ", "Bb2b      ", "bca  ",       PRIMITIVE,  7 }, /* 184 */
  { 38, "C2v^14", "A 2 -2          ", "A m m 2                        ", "A m m 2            ", "Amm2      ", "     ",          A_FACE,  7 }, /* 185 */
  { 38, "C2v^14", "B 2 -2          ", "B m m 2                        ", "B m m 2            ", "Bmm2      ", "ba-c ",       PRIMITIVE,  7 }, /* 186 */
  { 38, "C2v^14", "B -2 2          ", "B 2 m m                        ", "B 2 m m            ", "B2mm      ", "cab  ",       PRIMITIVE,  7 }, /* 187 */
  { 38, "C2v^14", "C -2 2          ", "C 2 m m                        ", "C 2 m m            ", "C2mm      ", "-cba ",          C_FACE,  7 }, /* 188 */
  { 38, "C2v^14", "C -2 -2         ", "C m 2 m                        ", "C m 2 m            ", "Cm2m      ", "bca  ",          C_FACE,  7 }, /* 189 */
  { 38, "C2v^14", "A -2 -2         ", "A m 2 m                        ", "A m 2 m            ", "Am2m      ", "a-cb ",          A_FACE,  7 }, /* 190 */
  { 39, "C2v^15", "A 2 -2c         ", "A e m 2                        ", "A e m 2            ", "Aem2      ", "     ",          A_FACE,  7 }, /* 191 */
  { 39, "C2v^15", "B 2 -2c         ", "B m e 2                        ", "B m e 2            ", "Bme2      ", "ba-c ",       PRIMITIVE,  7 }, /* 192 */
  { 39, "C2v^15", "B -2c 2         ", "B 2 e m                        ", "B 2 e m            ", "B2em      ", "cab  ",       PRIMITIVE,  7 }, /* 193 */
  { 39, "C2v^15", "C -2b 2         ", "C 2 m e                        ", "C 2 m e            ", "C2me      ", "-cba ",          C_FACE,  7 }, /* 194 */
  { 39, "C2v^15", "C -2b -2b       ", "C m 2 e                        ", "C m 2 e            ", "Cm2e      ", "bca  ",          C_FACE,  7 }, /* 195 */
  { 39, "C2v^15", "A -2c -2c       ", "A e 2 m                        ", "A e 2 m            ", "Ae2m      ", "a-cb ",          A_FACE,  7 }, /* 196 */
  { 40, "C2v^16", "A 2 -2a         ", "A m a 2                        ", "A m a 2            ", "Ama2      ", "     ",          A_FACE,  7 }, /* 197 */
  { 40, "C2v^16", "B 2 -2b         ", "B b m 2                        ", "B b m 2            ", "Bbm2      ", "ba-c ",       PRIMITIVE,  7 }, /* 198 */
  { 40, "C2v^16", "B -2b 2         ", "B 2 m b                        ", "B 2 m b            ", "B2mb      ", "cab  ",       PRIMITIVE,  7 }, /* 199 */
  { 40, "C2v^16", "C -2c 2         ", "C 2 c m                        ", "C 2 c m            ", "C2cm      ", "-cba ",          C_FACE,  7 }, /* 200 */
  { 40, "C2v^16", "C -2c -2c       ", "C c 2 m                        ", "C c 2 m            ", "Cc2m      ", "bca  ",          C_FACE,  7 }, /* 201 */
  { 40, "C2v^16", "A -2a -2a       ", "A m 2 a                        ", "A m 2 a            ", "Am2a      ", "a-cb ",          A_FACE,  7 }, /* 202 */
  { 41, "C2v^17", "A 2 -2ac        ", "A e a 2                        ", "A e a 2            ", "Aea2      ", "     ",          A_FACE,  7 }, /* 203 */
  { 41, "C2v^17", "B 2 -2bc        ", "B b e 2                        ", "B b e 2            ", "Bbe2      ", "ba-c ",       PRIMITIVE,  7 }, /* 204 */
  { 41, "C2v^17", "B -2bc 2        ", "B 2 e b                        ", "B 2 e b            ", "B2eb      ", "cab  ",       PRIMITIVE,  7 }, /* 205 */
  { 41, "C2v^17", "C -2bc 2        ", "C 2 c e                        ", "C 2 c e            ", "C2ce      ", "-cba ",          C_FACE,  7 }, /* 206 */
  { 41, "C2v^17", "C -2bc -2bc     ", "C c 2 e                        ", "C c 2 e            ", "Cc2e      ", "bca  ",          C_FACE,  7 }, /* 207 */
  { 41, "C2v^17", "A -2ac -2ac     ", "A e 2 a                        ", "A e 2 a            ", "Ae2a      ", "a-cb ",          A_FACE,  7 }, /* 208 */
  { 42, "C2v^18", "F 2 -2          ", "F m m 2                        ", "F m m 2            ", "Fmm2      ", "     ",            FACE,  7 }, /* 209 */
  { 42, "C2v^18", "F -2 2          ", "F 2 m m                        ", "F 2 m m            ", "F2mm      ", "cab  ",            FACE,  7 }, /* 210 */
  { 42, "C2v^18", "F -2 -2         ", "F m 2 m                        ", "F m 2 m            ", "Fm2m      ", "bca  ",            FACE,  7 }, /* 211 */
  { 43, "C2v^19", "F 2 -2d         ", "F d d 2                        ", "F d d 2            ", "Fdd2      ", "     ",            FACE,  7 }, /* 212 */
  { 43, "C2v^19", "F -2d 2         ", "F 2 d d                        ", "F 2 d d            ", "F2dd      ", "cab  ",            FACE,  7 }, /* 213 */
  { 43, "C2v^19", "F -2d -2d       ", "F d 2 d                        ", "F d 2 d            ", "Fd2d      ", "bca  ",            FACE,  7 }, /* 214 */
  { 44, "C2v^20", "I 2 -2          ", "I m m 2                        ", "I m m 2            ", "Imm2      ", "     ",            BODY,  7 }, /* 215 */
  { 44, "C2v^20", "I -2 2          ", "I 2 m m                        ", "I 2 m m            ", "I2mm      ", "cab  ",            BODY,  7 }, /* 216 */
  { 44, "C2v^20", "I -2 -2         ", "I m 2 m                        ", "I m 2 m            ", "Im2m      ", "bca  ",            BODY,  7 }, /* 217 */
  { 45, "C2v^21", "I 2 -2c         ", "I b a 2                        ", "I b a 2            ", "Iba2      ", "     ",            BODY,  7 }, /* 218 */
  { 45, "C2v^21", "I -2a 2         ", "I 2 c b                        ", "I 2 c b            ", "I2cb      ", "cab  ",            BODY,  7 }, /* 219 */
  { 45, "C2v^21", "I -2b -2b       ", "I c 2 a                        ", "I c 2 a            ", "Ic2a      ", "bca  ",            BODY,  7 }, /* 220 */
  { 46, "C2v^22", "I 2 -2a         ", "I m a 2                        ", "I m a 2            ", "Ima2      ", "     ",            BODY,  7 }, /* 221 */
  { 46, "C2v^22", "I 2 -2b         ", "I b m 2                        ", "I b m 2            ", "Ibm2      ", "ba-c ",            BODY,  7 }, /* 222 */
  { 46, "C2v^22", "I -2b 2         ", "I 2 m b                        ", "I 2 m b            ", "I2mb      ", "cab  ",            BODY,  7 }, /* 223 */
  { 46, "C2v^22", "I -2c 2         ", "I 2 c m                        ", "I 2 c m            ", "I2cm      ", "-cba ",            BODY,  7 }, /* 224 */
  { 46, "C2v^22", "I -2c -2c       ", "I c 2 m                        ", "I c 2 m            ", "Ic2m      ", "bca  ",            BODY,  7 }, /* 225 */
  { 46, "C2v^22", "I -2a -2a       ", "I m 2 a                        ", "I m 2 a            ", "Im2a      ", "a-cb ",            BODY,  7 }, /* 226 */
  { 47, "D2h^1 ", "-P 2 2          ", "P m m m                        ", "P 2/m 2/m 2/m      ", "Pmmm      ", "     ",       PRIMITIVE,  8 }, /* 227 */
  { 48, "D2h^2 ", "P 2 2 -1n       ", "P n n n                        ", "P 2/n 2/n 2/n      ", "Pnnn      ", "1    ",       PRIMITIVE,  8 }, /* 228 */
  { 48, "D2h^2 ", "-P 2ab 2bc      ", "P n n n                        ", "P 2/n 2/n 2/n      ", "Pnnn      ", "2    ",       PRIMITIVE,  8 }, /* 229 */
  { 49, "D2h^3 ", "-P 2 2c         ", "P c c m                        ", "P 2/c 2/c 2/m      ", "Pccm      ", "     ",       PRIMITIVE,  8 }, /* 230 */
  { 49, "D2h^3 ", "-P 2a 2         ", "P m a a                        ", "P 2/m 2/a 2/a      ", "Pmaa      ", "cab  ",       PRIMITIVE,  8 }, /* 231 */
  { 49, "D2h^3 ", "-P 2b 2b        ", "P b m b                        ", "P 2/b 2/m 2/b      ", "Pbmb      ", "bca  ",       PRIMITIVE,  8 }, /* 232 */
  { 50, "D2h^4 ", "P 2 2 -1ab      ", "P b a n                        ", "P 2/b 2/a 2/n      ", "Pban      ", "1    ",       PRIMITIVE,  8 }, /* 233 */
  { 50, "D2h^4 ", "-P 2ab 2b       ", "P b a n                        ", "P 2/b 2/a 2/n      ", "Pban      ", "2    ",       PRIMITIVE,  8 }, /* 234 */
  { 50, "D2h^4 ", "P 2 2 -1bc      ", "P n c b                        ", "P 2/n 2/c 2/b      ", "Pncb      ", "1cab ",       PRIMITIVE,  8 }, /* 235 */
  { 50, "D2h^4 ", "-P 2b 2bc       ", "P n c b                        ", "P 2/n 2/c 2/b      ", "Pncb      ", "2cab ",       PRIMITIVE,  8 }, /* 236 */
  { 50, "D2h^4 ", "P 2 2 -1ac      ", "P c n a                        ", "P 2/c 2/n 2/a      ", "Pcna      ", "1bca ",       PRIMITIVE,  8 }, /* 237 */
  { 50, "D2h^4 ", "-P 2a 2c        ", "P c n a                        ", "P 2/c 2/n 2/a      ", "Pcna      ", "2bca ",       PRIMITIVE,  8 }, /* 238 */
  { 51, "D2h^5 ", "-P 2a 2a        ", "P m m a                        ", "P 2_1/m 2/m 2/a    ", "Pmma      ", "     ",       PRIMITIVE,  8 }, /* 239 */
  { 51, "D2h^5 ", "-P 2b 2         ", "P m m b                        ", "P 2/m 2_1/m 2/b    ", "Pmmb      ", "ba-c ",       PRIMITIVE,  8 }, /* 240 */
  { 51, "D2h^5 ", "-P 2 2b         ", "P b m m                        ", "P 2/b 2_1/m 2/m    ", "Pbmm      ", "cab  ",       PRIMITIVE,  8 }, /* 241 */
  { 51, "D2h^5 ", "-P 2c 2c        ", "P c m m                        ", "P 2/c 2/m 2_1/m    ", "Pcmm      ", "-cba ",       PRIMITIVE,  8 }, /* 242 */
  { 51, "D2h^5 ", "-P 2c 2         ", "P m c m                        ", "P 2/m 2/c 2_1/m    ", "Pmcm      ", "bca  ",       PRIMITIVE,  8 }, /* 243 */
  { 51, "D2h^5 ", "-P 2 2a         ", "P m a m                        ", "P 2_1/m 2/a 2/m    ", "Pmam      ", "a-cb ",       PRIMITIVE,  8 }, /* 244 */
  { 52, "D2h^6 ", "-P 2a 2bc       ", "P n n a                        ", "P 2/n 2_1/n 2/a    ", "Pnna      ", "     ",       PRIMITIVE,  8 }, /* 245 */
  { 52, "D2h^6 ", "-P 2b 2n        ", "P n n b                        ", "P 2_1/n 2/n 2/b    ", "Pnnb      ", "ba-c ",       PRIMITIVE,  8 }, /* 246 */
  { 52, "D2h^6 ", "-P 2n 2b        ", "P b n n                        ", "P 2/b 2/n 2_1/n    ", "Pbnn      ", "cab  ",       PRIMITIVE,  8 }, /* 247 */
  { 52, "D2h^6 ", "-P 2ab 2c       ", "P c n n                        ", "P 2/c 2_1/n 2/n    ", "Pcnn      ", "-cba ",       PRIMITIVE,  8 }, /* 248 */
  { 52, "D2h^6 ", "-P 2ab 2n       ", "P n c n                        ", "P 2_1/n 2/c 2/n    ", "Pncn      ", "bca  ",       PRIMITIVE,  8 }, /* 249 */
  { 52, "D2h^6 ", "-P 2n 2bc       ", "P n a n                        ", "P 2/n 2/a 2_1/n    ", "Pnan      ", "a-cb ",       PRIMITIVE,  8 }, /* 250 */
  { 53, "D2h^7 ", "-P 2ac 2        ", "P m n a                        ", "P 2/m 2/n 2_1/a    ", "Pmna      ", "     ",       PRIMITIVE,  8 }, /* 251 */
  { 53, "D2h^7 ", "-P 2bc 2bc      ", "P n m b                        ", "P 2/n 2/m 2_1/b    ", "Pnmb      ", "ba-c ",       PRIMITIVE,  8 }, /* 252 */
  { 53, "D2h^7 ", "-P 2ab 2ab      ", "P b m n                        ", "P 2_1/b 2/m 2/n    ", "Pbmn      ", "cab  ",       PRIMITIVE,  8 }, /* 253 */
  { 53, "D2h^7 ", "-P 2 2ac        ", "P c n m                        ", "P 2_1/c 2/n 2/m    ", "Pcnm      ", "-cba ",       PRIMITIVE,  8 }, /* 254 */
  { 53, "D2h^7 ", "-P 2 2bc        ", "P n c m                        ", "P 2/n 2_1/c 2/m    ", "Pncm      ", "bca  ",       PRIMITIVE,  8 }, /* 255 */
  { 53, "D2h^7 ", "-P 2ab 2        ", "P m a n                        ", "P 2/m 2_1/a 2/n    ", "Pman      ", "a-cb ",       PRIMITIVE,  8 }, /* 256 */
  { 54, "D2h^8 ", "-P 2a 2ac       ", "P c c a                        ", "P 2_1/c 2/c 2/a    ", "Pcca      ", "     ",       PRIMITIVE,  8 }, /* 257 */
  { 54, "D2h^8 ", "-P 2b 2c        ", "P c c b                        ", "P 2/c 2_1/c 2/b    ", "Pccb      ", "ba-c ",       PRIMITIVE,  8 }, /* 258 */
  { 54, "D2h^8 ", "-P 2a 2b        ", "P b a a                        ", "P 2/b 2_1/a 2/a    ", "Pbaa      ", "cab  ",       PRIMITIVE,  8 }, /* 259 */
  { 54, "D2h^8 ", "-P 2ac 2c       ", "P c a a                        ", "P 2/c 2/a 2_1/a    ", "Pcaa      ", "-cba ",       PRIMITIVE,  8 }, /* 260 */
  { 54, "D2h^8 ", "-P 2bc 2b       ", "P b c b                        ", "P 2/b 2/c 2_1/b    ", "Pbcb      ", "bca  ",       PRIMITIVE,  8 }, /* 261 */
  { 54, "D2h^8 ", "-P 2b 2ab       ", "P b a b                        ", "P 2_1/b 2/a 2/b    ", "Pbab      ", "a-cb ",       PRIMITIVE,  8 }, /* 262 */
  { 55, "D2h^9 ", "-P 2 2ab        ", "P b a m                        ", "P 2_1/b 2_1/a 2/m  ", "Pbam      ", "     ",       PRIMITIVE,  8 }, /* 263 */
  { 55, "D2h^9 ", "-P 2bc 2        ", "P m c b                        ", "P 2/m 2_1/c 2_1/b  ", "Pmcb      ", "cab  ",       PRIMITIVE,  8 }, /* 264 */
  { 55, "D2h^9 ", "-P 2ac 2ac      ", "P c m a                        ", "P 2_1/c 2/m 2_1/a  ", "Pcma      ", "bca  ",       PRIMITIVE,  8 }, /* 265 */
  { 56, "D2h^10", "-P 2ab 2ac      ", "P c c n                        ", "P 2_1/c 2_1/c 2/n  ", "Pccn      ", "     ",       PRIMITIVE,  8 }, /* 266 */
  { 56, "D2h^10", "-P 2ac 2bc      ", "P n a a                        ", "P 2/n 2_1/a 2_1/a  ", "Pnaa      ", "cab  ",       PRIMITIVE,  8 }, /* 267 */
  { 56, "D2h^10", "-P 2bc 2ab      ", "P b n b                        ", "P 2_1/b 2/n 2_1/b  ", "Pbnb      ", "bca  ",       PRIMITIVE,  8 }, /* 268 */
  { 57, "D2h^11", "-P 2c 2b        ", "P b c m                        ", "P 2/b 2_1/c 2_1/m  ", "Pbcm      ", "     ",       PRIMITIVE,  8 }, /* 269 */
  { 57, "D2h^11", "-P 2c 2ac       ", "P c a m                        ", "P 2_1/c 2/a 2_1/m  ", "Pcam      ", "ba-c ",       PRIMITIVE,  8 }, /* 270 */
  { 57, "D2h^11", "-P 2ac 2a       ", "P m c a                        ", "P 2_1/m 2/c 2_1/a  ", "Pmca      ", "cab  ",       PRIMITIVE,  8 }, /* 271 */
  { 57, "D2h^11", "-P 2b 2a        ", "P m a b                        ", "P 2_1/m 2_1/a 2/b  ", "Pmab      ", "-cba ",       PRIMITIVE,  8 }, /* 272 */
  { 57, "D2h^11", "-P 2a 2ab       ", "P b m a                        ", "P 2_1/b 2_1/m 2/a  ", "Pbma      ", "bca  ",       PRIMITIVE,  8 }, /* 273 */
  { 57, "D2h^11", "-P 2bc 2c       ", "P c m b                        ", "P 2/c 2_1/m 2_1/b  ", "Pcmb      ", "a-cb ",       PRIMITIVE,  8 }, /* 274 */
  { 58, "D2h^12", "-P 2 2n         ", "P n n m                        ", "P 2_1/n 2_1/n 2/m  ", "Pnnm      ", "     ",       PRIMITIVE,  8 }, /* 275 */
  { 58, "D2h^12", "-P 2n 2         ", "P m n n                        ", "P 2/m 2_1/n 2_1/n  ", "Pmnn      ", "cab  ",       PRIMITIVE,  8 }, /* 276 */
  { 58, "D2h^12", "-P 2n 2n        ", "P n m n                        ", "P 2_1/n 2/m 2_1/n  ", "Pnmn      ", "bca  ",       PRIMITIVE,  8 }, /* 277 */
  { 59, "D2h^13", "P 2 2ab -1ab    ", "P m m n                        ", "P 2_1/m 2_1/m 2/n  ", "Pmmn      ", "1    ",       PRIMITIVE,  8 }, /* 278 */
  { 59, "D2h^13", "-P 2ab 2a       ", "P m m n                        ", "P 2_1/m 2_1/m 2/n  ", "Pmmn      ", "2    ",       PRIMITIVE,  8 }, /* 279 */
  { 59, "D2h^13", "P 2bc 2 -1bc    ", "P n m m                        ", "P 2/n 2_1/m 2_1/m  ", "Pnmm      ", "1cab ",       PRIMITIVE,  8 }, /* 280 */
  { 59, "D2h^13", "-P 2c 2bc       ", "P n m m                        ", "P 2/n 2_1/m 2_1/m  ", "Pnmm      ", "2cab ",       PRIMITIVE,  8 }, /* 281 */
  { 59, "D2h^13", "P 2ac 2ac -1ac  ", "P m n m                        ", "P 2_1/m 2/n 2_1/m  ", "Pmnm      ", "1bca ",       PRIMITIVE,  8 }, /* 282 */
  { 59, "D2h^13", "-P 2c 2a        ", "P m n m                        ", "P 2_1/m 2/n 2_1/m  ", "Pmnm      ", "2bca ",       PRIMITIVE,  8 }, /* 283 */
  { 60, "D2h^14", "-P 2n 2ab       ", "P b c n                        ", "P 2_1/b 2/c 2_1/n  ", "Pbcn      ", "     ",       PRIMITIVE,  8 }, /* 284 */
  { 60, "D2h^14", "-P 2n 2c        ", "P c a n                        ", "P 2/c 2_1/a 2_1/n  ", "Pcan      ", "ba-c ",       PRIMITIVE,  8 }, /* 285 */
  { 60, "D2h^14", "-P 2a 2n        ", "P n c a                        ", "P 2_1/n 2_1/c 2/a  ", "Pnca      ", "cab  ",       PRIMITIVE,  8 }, /* 286 */
  { 60, "D2h^14", "-P 2bc 2n       ", "P n a b                        ", "P 2_1/n 2/a 2_1/b  ", "Pnab      ", "-cba ",       PRIMITIVE,  8 }, /* 287 */
  { 60, "D2h^14", "-P 2ac 2b       ", "P b n a                        ", "P 2/b 2_1/n 2_1/a  ", "Pbna      ", "bca  ",       PRIMITIVE,  8 }, /* 288 */
  { 60, "D2h^14", "-P 2b 2ac       ", "P c n b                        ", "P 2_1/c 2_1/n 2/b  ", "Pcnb      ", "a-cb ",       PRIMITIVE,  8 }, /* 289 */
  { 61, "D2h^15", "-P 2ac 2ab      ", "P b c a                        ", "P 2_1/b 2_1/c 2_1/a", "Pbca      ", "     ",       PRIMITIVE,  8 }, /* 290 */
  { 61, "D2h^15", "-P 2bc 2ac      ", "P c a b                        ", "P 2_1/c 2_1/a 2_1/b", "Pcab      ", "ba-c ",       PRIMITIVE,  8 }, /* 291 */
  { 62, "D2h^16", "-P 2ac 2n       ", "P n m a                        ", "P 2_1/n 2_1/m 2_1/a", "Pnma      ", "     ",       PRIMITIVE,  8 }, /* 292 */
  { 62, "D2h^16", "-P 2bc 2a       ", "P m n b                        ", "P 2_1/m 2_1/n 2_1/b", "Pmnb      ", "ba-c ",       PRIMITIVE,  8 }, /* 293 */
  { 62, "D2h^16", "-P 2c 2ab       ", "P b n m                        ", "P 2_1/b 2_1/n 2_1/m", "Pbnm      ", "cab  ",       PRIMITIVE,  8 }, /* 294 */
  { 62, "D2h^16", "-P 2n 2ac       ", "P c m n                        ", "P 2_1/c 2_1/m 2_1/n", "Pcmn      ", "-cba ",       PRIMITIVE,  8 }, /* 295 */
  { 62, "D2h^16", "-P 2n 2a        ", "P m c n                        ", "P 2_1/m 2_1/c 2_1/n", "Pmcn      ", "bca  ",       PRIMITIVE,  8 }, /* 296 */
  { 62, "D2h^16", "-P 2c 2n        ", "P n a m                        ", "P 2_1/n 2_1/a 2_1/m", "Pnam      ", "a-cb ",       PRIMITIVE,  8 }, /* 297 */
  { 63, "D2h^17", "-C 2c 2         ", "C m c m                        ", "C 2/m 2/c 2_1/m    ", "Cmcm      ", "     ",          C_FACE,  8 }, /* 298 */
  { 63, "D2h^17", "-C 2c 2c        ", "C c m m                        ", "C 2/c 2/m 2_1/m    ", "Ccmm      ", "ba-c ",          C_FACE,  8 }, /* 299 */
  { 63, "D2h^17", "-A 2a 2a        ", "A m m a                        ", "A 2_1/m 2/m 2/a    ", "Amma      ", "cab  ",          A_FACE,  8 }, /* 300 */
  { 63, "D2h^17", "-A 2 2a         ", "A m a m                        ", "A 2_1/m 2/a 2/m    ", "Amam      ", "-cba ",          A_FACE,  8 }, /* 301 */
  { 63, "D2h^17", "-B 2 2b         ", "B b m m                        ", "B 2/b 2_1/m 2/m    ", "Bbmm      ", "bca  ",       PRIMITIVE,  8 }, /* 302 */
  { 63, "D2h^17", "-B 2b 2         ", "B m m b                        ", "B 2/m 2_1/m 2/b    ", "Bmmb      ", "a-cb ",       PRIMITIVE,  8 }, /* 303 */
  { 64, "D2h^18", "-C 2bc 2        ", "C m c e                        ", "C 2/m 2/c 2_1/e    ", "Cmce      ", "     ",          C_FACE,  8 }, /* 304 */
  { 64, "D2h^18", "-C 2bc 2bc      ", "C c m e                        ", "C 2/c 2/m 2_1/e    ", "Ccme      ", "ba-c ",          C_FACE,  8 }, /* 305 */
  { 64, "D2h^18", "-A 2ac 2ac      ", "A e m a                        ", "A 2_1/e 2/m 2/a    ", "Aema      ", "cab  ",          A_FACE,  8 }, /* 306 */
  { 64, "D2h^18", "-A 2 2ac        ", "A e a m                        ", "A 2_1/e 2/a 2/m    ", "Aeam      ", "-cba ",          A_FACE,  8 }, /* 307 */
  { 64, "D2h^18", "-B 2 2bc        ", "B b e m                        ", "B 2/b 2_1/e 2/m    ", "Bbem      ", "bca  ",       PRIMITIVE,  8 }, /* 308 */
  { 64, "D2h^18", "-B 2bc 2        ", "B m e b                        ", "B 2/m 2_1/e 2/b    ", "Bmeb      ", "a-cb ",       PRIMITIVE,  8 }, /* 309 */
  { 65, "D2h^19", "-C 2 2          ", "C m m m                        ", "C 2/m 2/m 2/m      ", "Cmmm      ", "     ",          C_FACE,  8 }, /* 310 */
  { 65, "D2h^19", "-A 2 2          ", "A m m m                        ", "A 2/m 2/m 2/m      ", "Ammm      ", "cab  ",          A_FACE,  8 }, /* 311 */
  { 65, "D2h^19", "-B 2 2          ", "B m m m                        ", "B 2/m 2/m 2/m      ", "Bmmm      ", "bca  ",       PRIMITIVE,  8 }, /* 312 */
  { 66, "D2h^20", "-C 2 2c         ", "C c c m                        ", "C 2/c 2/c 2/m      ", "Cccm      ", "     ",          C_FACE,  8 }, /* 313 */
  { 66, "D2h^20", "-A 2a 2         ", "A m a a                        ", "A 2/m 2/a 2/a      ", "Amaa      ", "cab  ",          A_FACE,  8 }, /* 314 */
  { 66, "D2h^20", "-B 2b 2b        ", "B b m b                        ", "B 2/b 2/m 2/b      ", "Bbmb      ", "bca  ",       PRIMITIVE,  8 }, /* 315 */
  { 67, "D2h^21", "-C 2b 2         ", "C m m e                        ", "C 2/m 2/m 2/e      ", "Cmme      ", "     ",          C_FACE,  8 }, /* 316 */
  { 67, "D2h^21", "-C 2b 2b        ", "C m m e                        ", "C 2/m 2/m 2/e      ", "Cmme      ", "ba-c ",          C_FACE,  8 }, /* 317 */
  { 67, "D2h^21", "-A 2c 2c        ", "A e m m                        ", "A 2/e 2/m 2/m      ", "Aemm      ", "cab  ",          A_FACE,  8 }, /* 318 */
  { 67, "D2h^21", "-A 2 2c         ", "A e m m                        ", "A 2/e 2/m 2/m      ", "Aemm      ", "-cba ",          A_FACE,  8 }, /* 319 */
  { 67, "D2h^21", "-B 2 2c         ", "B m e m                        ", "B 2/m 2/e 2/m      ", "Bmem      ", "bca  ",       PRIMITIVE,  8 }, /* 320 */
  { 67, "D2h^21", "-B 2c 2         ", "B m e m                        ", "B 2/m 2/e 2/m      ", "Bmem      ", "a-cb ",       PRIMITIVE,  8 }, /* 321 */
  { 68, "D2h^22", "C 2 2 -1bc      ", "C c c e                        ", "C 2/c 2/c 2/e      ", "Ccce      ", "1    ",          C_FACE,  8 }, /* 322 */
  { 68, "D2h^22", "-C 2b 2bc       ", "C c c e                        ", "C 2/c 2/c 2/e      ", "Ccce      ", "2    ",          C_FACE,  8 }, /* 323 */
  { 68, "D2h^22", "C 2 2 -1bc      ", "C c c e                        ", "C 2/c 2/c 2/e      ", "Ccce      ", "1ba-c",          C_FACE,  8 }, /* 324 */
  { 68, "D2h^22", "-C 2b 2c        ", "C c c e                        ", "C 2/c 2/c 2/e      ", "Ccce      ", "2ba-c",          C_FACE,  8 }, /* 325 */
  { 68, "D2h^22", "A 2 2 -1ac      ", "A e a a                        ", "A 2/e 2/a 2/a      ", "Aeaa      ", "1cab ",          A_FACE,  8 }, /* 326 */
  { 68, "D2h^22", "-A 2a 2c        ", "A e a a                        ", "A 2/e 2/a 2/a      ", "Aeaa      ", "2cab ",          A_FACE,  8 }, /* 327 */
  { 68, "D2h^22", "A 2 2 -1ac      ", "A e a a                        ", "A 2/e 2/a 2/a      ", "Aeaa      ", "1-cba",          A_FACE,  8 }, /* 328 */
  { 68, "D2h^22", "-A 2ac 2c       ", "A e a a                        ", "A 2/e 2/a 2/a      ", "Aeaa      ", "2-cba",          A_FACE,  8 }, /* 329 */
  { 68, "D2h^22", "B 2 2 -1bc      ", "B b e b                        ", "B 2/b 2/e 2/b      ", "Bbeb      ", "1bca ",       PRIMITIVE,  8 }, /* 330 */
  { 68, "D2h^22", "-B 2bc 2b       ", "B b c b                        ", "B 2/b 2/e 2/b      ", "Bbcb      ", "2bca ",       PRIMITIVE,  8 }, /* 331 */
  { 68, "D2h^22", "B 2 2 -1bc      ", "B b e b                        ", "B 2/b 2/e 2/b      ", "Bbeb      ", "1a-cb",       PRIMITIVE,  8 }, /* 332 */
  { 68, "D2h^22", "-B 2b 2bc       ", "B b e b                        ", "B 2/b 2/e 2/b      ", "Bbeb      ", "2a-cb",       PRIMITIVE,  8 }, /* 333 */
  { 69, "D2h^23", "-F 2 2          ", "F m m m                        ", "F 2/m 2/m 2/m      ", "Fmmm      ", "     ",            FACE,  8 }, /* 334 */
  { 70, "D2h^24", "F 2 2 -1d       ", "F d d d                        ", "F 2/d 2/d 2/d      ", "Fddd      ", "1    ",            FACE,  8 }, /* 335 */
  { 70, "D2h^24", "-F 2uv 2vw      ", "F d d d                        ", "F 2/d 2/d 2/d      ", "Fddd      ", "2    ",            FACE,  8 }, /* 336 */
  { 71, "D2h^25", "-I 2 2          ", "I m m m                        ", "I 2/m 2/m 2/m      ", "Immm      ", "     ",            BODY,  8 }, /* 337 */
  { 72, "D2h^26", "-I 2 2c         ", "I b a m                        ", "I 2/b 2/a 2/m      ", "Ibam      ", "     ",            BODY,  8 }, /* 338 */
  { 72, "D2h^26", "-I 2a 2         ", "I m c b                        ", "I 2/m 2/c 2/b      ", "Imcb      ", "cab  ",            BODY,  8 }, /* 339 */
  { 72, "D2h^26", "-I 2b 2b        ", "I c m a                        ", "I 2/c 2/m 2/a      ", "Icma      ", "bca  ",            BODY,  8 }, /* 340 */
  { 73, "D2h^27", "-I 2b 2c        ", "I b c a                        ", "I 2/b 2/c 2/a      ", "Ibca      ", "     ",            BODY,  8 }, /* 341 */
  { 73, "D2h^27", "-I 2a 2b        ", "I c a b                        ", "I 2/c 2/a 2/b      ", "Icab      ", "ba-c ",            BODY,  8 }, /* 342 */
  { 74, "D2h^28", "-I 2b 2         ", "I m m a                        ", "I 2/m 2/m 2/a      ", "Imma      ", "     ",            BODY,  8 }, /* 343 */
  { 74, "D2h^28", "-I 2a 2a        ", "I m m b                        ", "I 2/m 2/m 2/b      ", "Immb      ", "ba-c ",            BODY,  8 }, /* 344 */
  { 74, "D2h^28", "-I 2c 2c        ", "I b m m                        ", "I 2/b 2/m 2/m      ", "Ibmm      ", "cab  ",            BODY,  8 }, /* 345 */
  { 74, "D2h^28", "-I 2 2b         ", "I c m m                        ", "I 2/c 2/m 2/m      ", "Icmm      ", "-cba ",            BODY,  8 }, /* 346 */
  { 74, "D2h^28", "-I 2 2a         ", "I m c m                        ", "I 2/m 2/c 2/m      ", "Imcm      ", "bca  ",            BODY,  8 }, /* 347 */
  { 74, "D2h^28", "-I 2c 2         ", "I m a m                        ", "I 2/m 2/a 2/m      ", "Imam      ", "a-cb ",            BODY,  8 }, /* 348 */
  { 75, "C4^1  ", "P 4             ", "P 4                            ", "P 4                ", "P4        ", "     ",       PRIMITIVE,  9 }, /* 349 */
  { 76, "C4^2  ", "P 4w            ", "P 4_1                          ", "P 4_1              ", "P4_1      ", "     ",       PRIMITIVE,  9 }, /* 350 */
  { 77, "C4^3  ", "P 4c            ", "P 4_2                          ", "P 4_2              ", "P4_2      ", "     ",       PRIMITIVE,  9 }, /* 351 */
  { 78, "C4^4  ", "P 4cw           ", "P 4_3                          ", "P 4_3              ", "P4_3      ", "     ",       PRIMITIVE,  9 }, /* 352 */
  { 79, "C4^5  ", "I 4             ", "I 4                            ", "I 4                ", "I4        ", "     ",            BODY,  9 }, /* 353 */
  { 80, "C4^6  ", "I 4bw           ", "I 4_1                          ", "I 4_1              ", "I4_1      ", "     ",            BODY,  9 }, /* 354 */
  { 81, "S4^1  ", "P -4            ", "P -4                           ", "P -4               ", "P-4       ", "     ",       PRIMITIVE, 10 }, /* 355 */
  { 82, "S4^2  ", "I -4            ", "I -4                           ", "I -4               ", "I-4       ", "     ",            BODY, 10 }, /* 356 */
  { 83, "C4h^1 ", "-P 4            ", "P 4/m                          ", "P 4/m              ", "P4/m      ", "     ",       PRIMITIVE, 11 }, /* 357 */
  { 84, "C4h^2 ", "-P 4c           ", "P 4_2/m                        ", "P 4_2/m            ", "P4_2/m    ", "     ",       PRIMITIVE, 11 }, /* 358 */
  { 85, "C4h^3 ", "P 4ab -1ab      ", "P 4/n                          ", "P 4/n              ", "P4/n      ", "1    ",       PRIMITIVE, 11 }, /* 359 */
  { 85, "C4h^3 ", "-P 4a           ", "P 4/n                          ", "P 4/n              ", "P4/n      ", "2    ",       PRIMITIVE, 11 }, /* 360 */
  { 86, "C4h^4 ", "P 4n -1n        ", "P 4_2/n                        ", "P 4_2/n            ", "P4_2/n    ", "1    ",       PRIMITIVE, 11 }, /* 361 */
  { 86, "C4h^4 ", "-P 4bc          ", "P 4_2/n                        ", "P 4_2/n            ", "P4_2/n    ", "2    ",       PRIMITIVE, 11 }, /* 362 */
  { 87, "C4h^5 ", "-I 4            ", "I 4/m                          ", "I 4/m              ", "I4/m      ", "     ",            BODY, 11 }, /* 363 */
  { 88, "C4h^6 ", "I 4bw -1bw      ", "I 4_1/a                        ", "I 4_1/a            ", "I4_1/a    ", "1    ",            BODY, 11 }, /* 364 */
  { 88, "C4h^6 ", "-I 4ad          ", "I 4_1/a                        ", "I 4_1/a            ", "I4_1/a    ", "2    ",            BODY, 11 }, /* 365 */
  { 89, "D4^1  ", "P 4 2           ", "P 4 2 2                        ", "P 4 2 2            ", "P422      ", "     ",       PRIMITIVE, 12 }, /* 366 */
  { 90, "D4^2  ", "P 4ab 2ab       ", "P 4 2_1 2                      ", "P 4 2_1 2          ", "P42_12    ", "     ",       PRIMITIVE, 12 }, /* 367 */
  { 91, "D4^3  ", "P 4w 2c         ", "P 4_1 2 2                      ", "P 4_1 2 2          ", "P4_122    ", "     ",       PRIMITIVE, 12 }, /* 368 */
  { 92, "D4^4  ", "P 4abw 2nw      ", "P 4_1 2_1 2                    ", "P 4_1 2_1 2        ", "P4_12_12  ", "     ",       PRIMITIVE, 12 }, /* 369 */
  { 93, "D4^5  ", "P 4c 2          ", "P 4_2 2 2                      ", "P 4_2 2 2          ", "P4_222    ", "     ",       PRIMITIVE, 12 }, /* 370 */
  { 94, "D4^6  ", "P 4n 2n         ", "P 4_2 2_1 2                    ", "P 4_2 2_1 2        ", "P4_22_12  ", "     ",       PRIMITIVE, 12 }, /* 371 */
  { 95, "D4^7  ", "P 4cw 2c        ", "P 4_3 2 2                      ", "P 4_3 2 2          ", "P4_322    ", "     ",       PRIMITIVE, 12 }, /* 372 */
  { 96, "D4^8  ", "P 4nw 2abw      ", "P 4_3 2_1 2                    ", "P 4_3 2_1 2        ", "P4_32_12  ", "     ",       PRIMITIVE, 12 }, /* 373 */
  { 97, "D4^9  ", "I 4 2           ", "I 4 2 2                        ", "I 4 2 2            ", "I422      ", "     ",            BODY, 12 }, /* 374 */
  { 98, "D4^10 ", "I 4bw 2bw       ", "I 4_1 2 2                      ", "I 4_1 2 2          ", "I4_122    ", "     ",            BODY, 12 }, /* 375 */
  { 99, "C4v^1 ", "P 4 -2          ", "P 4 m m                        ", "P 4 m m            ", "P4mm      ", "     ",       PRIMITIVE, 13 }, /* 376 */
  {100, "C4v^2 ", "P 4 -2ab        ", "P 4 b m                        ", "P 4 b m            ", "P4bm      ", "     ",       PRIMITIVE, 13 }, /* 377 */
  {101, "C4v^3 ", "P 4c -2c        ", "P 4_2 c m                      ", "P 4_2 c m          ", "P4_2cm    ", "     ",       PRIMITIVE, 13 }, /* 378 */
  {102, "C4v^4 ", "P 4n -2n        ", "P 4_2 n m                      ", "P 4_2 n m          ", "P4_2nm    ", "     ",       PRIMITIVE, 13 }, /* 379 */
  {103, "C4v^5 ", "P 4 -2c         ", "P 4 c c                        ", "P 4 c c            ", "P4cc      ", "     ",       PRIMITIVE, 13 }, /* 380 */
  {104, "C4v^6 ", "P 4 -2n         ", "P 4 n c                        ", "P 4 n c            ", "P4nc      ", "     ",       PRIMITIVE, 13 }, /* 381 */
  {105, "C4v^7 ", "P 4c -2         ", "P 4_2 m c                      ", "P 4_2 m c          ", "P4_2mc    ", "     ",       PRIMITIVE, 13 }, /* 382 */
  {106, "C4v^8 ", "P 4c -2ab       ", "P 4_2 b c                      ", "P 4_2 b c          ", "P4_2bc    ", "     ",       PRIMITIVE, 13 }, /* 383 */
  {107, "C4v^9 ", "I 4 -2          ", "I 4 m m                        ", "I 4 m m            ", "I4mm      ", "     ",            BODY, 13 }, /* 384 */
  {108, "C4v^10", "I 4 -2c         ", "I 4 c m                        ", "I 4 c m            ", "I4cm      ", "     ",            BODY, 13 }, /* 385 */
  {109, "C4v^11", "I 4bw -2        ", "I 4_1 m d                      ", "I 4_1 m d          ", "I4_1md    ", "     ",            BODY, 13 }, /* 386 */
  {110, "C4v^12", "I 4bw -2c       ", "I 4_1 c d                      ", "I 4_1 c d          ", "I4_1cd    ", "     ",            BODY, 13 }, /* 387 */
  {111, "D2d^1 ", "P -4 2          ", "P -4 2 m                       ", "P -4 2 m           ", "P-42m     ", "     ",       PRIMITIVE, 14 }, /* 388 */
  {112, "D2d^2 ", "P -4 2c         ", "P -4 2 c                       ", "P -4 2 c           ", "P-42c     ", "     ",       PRIMITIVE, 14 }, /* 389 */
  {113, "D2d^3 ", "P -4 2ab        ", "P -4 2_1 m                     ", "P -4 2_1 m         ", "P-42_1m   ", "     ",       PRIMITIVE, 14 }, /* 390 */
  {114, "D2d^4 ", "P -4 2n         ", "P -4 2_1 c                     ", "P -4 2_1 c         ", "P-42_1c   ", "     ",       PRIMITIVE, 14 }, /* 391 */
  {115, "D2d^5 ", "P -4 -2         ", "P -4 m 2                       ", "P -4 m 2           ", "P-4m2     ", "     ",       PRIMITIVE, 14 }, /* 392 */
  {116, "D2d^6 ", "P -4 -2c        ", "P -4 c 2                       ", "P -4 c 2           ", "P-4c2     ", "     ",       PRIMITIVE, 14 }, /* 393 */
  {117, "D2d^7 ", "P -4 -2ab       ", "P -4 b 2                       ", "P -4 b 2           ", "P-4b2     ", "     ",       PRIMITIVE, 14 }, /* 394 */
  {118, "D2d^8 ", "P -4 -2n        ", "P -4 n 2                       ", "P -4 n 2           ", "P-4n2     ", "     ",       PRIMITIVE, 14 }, /* 395 */
  {119, "D2d^9 ", "I -4 -2         ", "I -4 m 2                       ", "I -4 m 2           ", "I-4m2     ", "     ",            BODY, 14 }, /* 396 */
  {120, "D2d^10", "I -4 -2c        ", "I -4 c 2                       ", "I -4 c 2           ", "I-4c2     ", "     ",            BODY, 14 }, /* 397 */
  {121, "D2d^11", "I -4 2          ", "I -4 2 m                       ", "I -4 2 m           ", "I-42m     ", "     ",            BODY, 14 }, /* 398 */
  {122, "D2d^12", "I -4 2bw        ", "I -4 2 d                       ", "I -4 2 d           ", "I-42d     ", "     ",            BODY, 14 }, /* 399 */
  {123, "D4h^1 ", "-P 4 2          ", "P 4/m m m                      ", "P 4/m 2/m 2/m      ", "P4/mmm    ", "     ",       PRIMITIVE, 15 }, /* 400 */
  {124, "D4h^2 ", "-P 4 2c         ", "P 4/m c c                      ", "P 4/m 2/c 2/c      ", "P4/mcc    ", "     ",       PRIMITIVE, 15 }, /* 401 */
  {125, "D4h^3 ", "P 4 2 -1ab      ", "P 4/n b m                      ", "P 4/n 2/b 2/m      ", "P4/nbm    ", "1    ",       PRIMITIVE, 15 }, /* 402 */
  {125, "D4h^3 ", "-P 4a 2b        ", "P 4/n b m                      ", "P 4/n 2/b 2/m      ", "P4/nbm    ", "2    ",       PRIMITIVE, 15 }, /* 403 */
  {126, "D4h^4 ", "P 4 2 -1n       ", "P 4/n n c                      ", "P 4/n 2/n 2/c      ", "P4/nnc    ", "1    ",       PRIMITIVE, 15 }, /* 404 */
  {126, "D4h^4 ", "-P 4a 2bc       ", "P 4/n n c                      ", "P 4/n 2/n 2/c      ", "P4/nnc    ", "2    ",       PRIMITIVE, 15 }, /* 405 */
  {127, "D4h^5 ", "-P 4 2ab        ", "P 4/m b m                      ", "P 4/m 2_1/b m      ", "P4/mbm    ", "     ",       PRIMITIVE, 15 }, /* 406 */
  {128, "D4h^6 ", "-P 4 2n         ", "P 4/m n c                      ", "P 4/m 2_1/n c      ", "P4/mnc    ", "     ",       PRIMITIVE, 15 }, /* 407 */
  {129, "D4h^7 ", "P 4ab 2ab -1ab  ", "P 4/n m m                      ", "P 4/n 2_1/m m      ", "P4/nmm    ", "1    ",       PRIMITIVE, 15 }, /* 408 */
  {129, "D4h^7 ", "-P 4a 2a        ", "P 4/n m m                      ", "P 4/n 2_1/m m      ", "P4/nmm    ", "2    ",       PRIMITIVE, 15 }, /* 409 */
  {130, "D4h^8 ", "P 4ab 2n -1ab   ", "P 4/n c c                      ", "P 4/n 2_1/c c      ", "P4/ncc    ", "1    ",       PRIMITIVE, 15 }, /* 410 */
  {130, "D4h^8 ", "-P 4a 2ac       ", "P 4/n c c                      ", "P 4/n 2_1/c c      ", "P4/ncc    ", "2    ",       PRIMITIVE, 15 }, /* 411 */
  {131, "D4h^9 ", "-P 4c 2         ", "P 4_2/m m c                    ", "P 4_2/m 2/m 2/c    ", "P4_2/mmc  ", "     ",       PRIMITIVE, 15 }, /* 412 */
  {132, "D4h^10", "-P 4c 2c        ", "P 4_2/m c m                    ", "P 4_2/m 2/c 2/m    ", "P4_2/mcm  ", "     ",       PRIMITIVE, 15 }, /* 413 */
  {133, "D4h^11", "P 4n 2c -1n     ", "P 4_2/n b c                    ", "P 4_2/n 2/b 2/c    ", "P4_2/nbc  ", "1    ",       PRIMITIVE, 15 }, /* 414 */
  {133, "D4h^11", "-P 4ac 2b       ", "P 4_2/n b c                    ", "P 4_2/n 2/b 2/c    ", "P4_2/nbc  ", "2    ",       PRIMITIVE, 15 }, /* 415 */
  {134, "D4h^12", "P 4n 2 -1n      ", "P 4_2/n n m                    ", "P 4_2/n 2/n 2/m    ", "P4_2/nnm  ", "1    ",       PRIMITIVE, 15 }, /* 416 */
  {134, "D4h^12", "-P 4ac 2bc      ", "P 4_2/n n m                    ", "P 4_2/n 2/n 2/m    ", "P4_2/nnm  ", "2    ",       PRIMITIVE, 15 }, /* 417 */
  {135, "D4h^13", "-P 4c 2ab       ", "P 4_2/m b c                    ", "P 4_2/m 2_1/b 2/c  ", "P4_2/mbc  ", "     ",       PRIMITIVE, 15 }, /* 418 */
  {136, "D4h^14", "-P 4n 2n        ", "P 4_2/m n m                    ", "P 4_2/m 2_1/n 2/m  ", "P4_2/mnm  ", "     ",       PRIMITIVE, 15 }, /* 419 */
  {137, "D4h^15", "P 4n 2n -1n     ", "P 4_2/n m c                    ", "P 4_2/n 2_1/m 2/c  ", "P4_2/nmc  ", "1    ",       PRIMITIVE, 15 }, /* 420 */
  {137, "D4h^15", "-P 4ac 2a       ", "P 4_2/n m c                    ", "P 4_2/n 2_1/m 2/c  ", "P4_2/nmc  ", "2    ",       PRIMITIVE, 15 }, /* 421 */
  {138, "D4h^16", "P 4n 2ab -1n    ", "P 4_2/n c m                    ", "P 4_2/n 2_1/c 2/m  ", "P4_2/ncm  ", "1    ",       PRIMITIVE, 15 }, /* 422 */
  {138, "D4h^16", "-P 4ac 2ac      ", "P 4_2/n c m                    ", "P 4_2/n 2_1/c 2/m  ", "P4_2/ncm  ", "2    ",       PRIMITIVE, 15 }, /* 423 */
  {139, "D4h^17", "-I 4 2          ", "I 4/m m m                      ", "I 4/m 2/m 2/m      ", "I4/mmm    ", "     ",            BODY, 15 }, /* 424 */
  {140, "D4h^18", "-I 4 2c         ", "I 4/m c m                      ", "I 4/m 2/c 2/m      ", "I4/mcm    ", "     ",            BODY, 15 }, /* 425 */
  {141, "D4h^19", "I 4bw 2bw -1bw  ", "I 4_1/a m d                    ", "I 4_1/a 2/m 2/d    ", "I4_1/amd  ", "1    ",            BODY, 15 }, /* 426 */
  {141, "D4h^19", "-I 4bd 2        ", "I 4_1/a m d                    ", "I 4_1/a 2/m 2/d    ", "I4_1/amd  ", "2    ",            BODY, 15 }, /* 427 */
  {142, "D4h^20", "I 4bw 2aw -1bw  ", "I 4_1/a c d                    ", "I 4_1/a 2/c 2/d    ", "I4_1/acd  ", "1    ",            BODY, 15 }, /* 428 */
  {142, "D4h^20", "-I 4bd 2c       ", "I 4_1/a c d                    ", "I 4_1/a 2/c 2/d    ", "I4_1/acd  ", "2    ",            BODY, 15 }, /* 429 */
  {143, "C3^1  ", "P 3             ", "P 3                            ", "P 3                ", "P3        ", "     ",       PRIMITIVE, 16 }, /* 430 */
  {144, "C3^2  ", "P 31            ", "P 3_1                          ", "P 3_1              ", "P3_1      ", "     ",       PRIMITIVE, 16 }, /* 431 */
  {145, "C3^3  ", "P 32            ", "P 3_2                          ", "P 3_2              ", "P3_2      ", "     ",       PRIMITIVE, 16 }, /* 432 */
  {146, "C3^4  ", "R 3             ", "R 3                            ", "R 3                ", "R3        ", "H    ",       PRIMITIVE, 16 }, /* 433 */
  {146, "C3^4  ", "P 3*            ", "R 3                            ", "R 3                ", "R3        ", "R    ",        R_CENTER, 16 }, /* 434 */
  {147, "C3i^1 ", "-P 3            ", "P -3                           ", "P -3               ", "P-3       ", "     ",       PRIMITIVE, 17 }, /* 435 */
  {148, "C3i^2 ", "-R 3            ", "R -3                           ", "R -3               ", "R-3       ", "H    ",       PRIMITIVE, 17 }, /* 436 */
  {148, "C3i^2 ", "-P 3*           ", "R -3                           ", "R -3               ", "R-3       ", "R    ",        R_CENTER, 17 }, /* 437 */
  {149, "D3^1  ", "P 3 2           ", "P 3 1 2                        ", "P 3 1 2            ", "P312      ", "     ",       PRIMITIVE, 18 }, /* 438 */
  {150, "D3^2  ", "P 3 2=          ", "P 3 2 1                        ", "P 3 2 1            ", "P321      ", "     ",       PRIMITIVE, 18 }, /* 439 */
  {151, "D3^3  ", "P 31 2c (0 0 1) ", "P 3_1 1 2                      ", "P 3_1 1 2          ", "P3_112    ", "     ",       PRIMITIVE, 18 }, /* 440 */
  {152, "D3^4  ", "P 31 2=         ", "P 3_1 2 1                      ", "P 3_1 2 1          ", "P3_121    ", "     ",       PRIMITIVE, 18 }, /* 441 */
  {153, "D3^5  ", "P 32 2c (0 0 -1)", "P 3_2 1 2                      ", "P 3_2 1 2          ", "P3_212    ", "     ",       PRIMITIVE, 18 }, /* 442 */
  {154, "D3^6  ", "P 32 2=         ", "P 3_2 2 1                      ", "P 3_2 2 1          ", "P3_221    ", "     ",       PRIMITIVE, 18 }, /* 443 */
  {155, "D3^7  ", "R 3 2=          ", "R 3 2                          ", "R 3 2              ", "R32       ", "H    ",       PRIMITIVE, 18 }, /* 444 */
  {155, "D3^7  ", "P 3* 2          ", "R 3 2                          ", "R 3 2              ", "R32       ", "R    ",        R_CENTER, 18 }, /* 445 */
  {156, "C3v^1 ", "P 3 -2=         ", "P 3 m 1                        ", "P 3 m 1            ", "P3m1      ", "     ",       PRIMITIVE, 19 }, /* 446 */
  {157, "C3v^2 ", "P 3 -2          ", "P 3 1 m                        ", "P 3 1 m            ", "P31m      ", "     ",       PRIMITIVE, 19 }, /* 447 */
  {158, "C3v^3 ", "P 3 -2=c        ", "P 3 c 1                        ", "P 3 c 1            ", "P3c1      ", "     ",       PRIMITIVE, 19 }, /* 448 */
  {159, "C3v^4 ", "P 3 -2c         ", "P 3 1 c                        ", "P 3 1 c            ", "P31c      ", "     ",       PRIMITIVE, 19 }, /* 449 */
  {160, "C3v^5 ", "R 3 -2=         ", "R 3 m                          ", "R 3 m              ", "R3m       ", "H    ",       PRIMITIVE, 19 }, /* 450 */
  {160, "C3v^5 ", "P 3* -2         ", "R 3 m                          ", "R 3 m              ", "R3m       ", "R    ",        R_CENTER, 19 }, /* 451 */
  {161, "C3v^6 ", "R 3 -2=c        ", "R 3 c                          ", "R 3 c              ", "R3c       ", "H    ",       PRIMITIVE, 19 }, /* 452 */
  {161, "C3v^6 ", "P 3* -2n        ", "R 3 c                          ", "R 3 c              ", "R3c       ", "R    ",        R_CENTER, 19 }, /* 453 */
  {162, "D3d^1 ", "-P 3 2          ", "P -3 1 m                       ", "P -3 1 2/m         ", "P-31m     ", "     ",       PRIMITIVE, 20 }, /* 454 */
  {163, "D3d^2 ", "-P 3 2c         ", "P -3 1 c                       ", "P -3 1 2/c         ", "P-31c     ", "     ",       PRIMITIVE, 20 }, /* 455 */
  {164, "D3d^3 ", "-P 3 2=         ", "P -3 m 1                       ", "P -3 2/m 1         ", "P-3m1     ", "     ",       PRIMITIVE, 20 }, /* 456 */
  {165, "D3d^4 ", "-P 3 2=c        ", "P -3 c 1                       ", "P -3 2/c 1         ", "P-3c1     ", "     ",       PRIMITIVE, 20 }, /* 457 */
  {166, "D3d^5 ", "-R 3 2=         ", "R -3 m                         ", "R -3 2/m           ", "R-3m      ", "H    ",       PRIMITIVE, 20 }, /* 458 */
  {166, "D3d^5 ", "-P 3* 2         ", "R -3 m                         ", "R -3 2/m           ", "R-3m      ", "R    ",        R_CENTER, 20 }, /* 459 */
  {167, "D3d^6 ", "-R 3 2=c        ", "R -3 c                         ", "R -3 2/c           ", "R-3c      ", "H    ",       PRIMITIVE, 20 }, /* 460 */
  {167, "D3d^6 ", "-P 3* 2n        ", "R -3 c                         ", "R -3 2/c           ", "R-3c      ", "R    ",        R_CENTER, 20 }, /* 461 */
  {168, "C6^1  ", "P 6             ", "P 6                            ", "P 6                ", "P6        ", "     ",       PRIMITIVE, 21 }, /* 462 */
  {169, "C6^2  ", "P 61            ", "P 6_1                          ", "P 6_1              ", "P6_1      ", "     ",       PRIMITIVE, 21 }, /* 463 */
  {170, "C6^3  ", "P 65            ", "P 6_5                          ", "P 6_5              ", "P6_5      ", "     ",       PRIMITIVE, 21 }, /* 464 */
  {171, "C6^4  ", "P 62            ", "P 6_2                          ", "P 6_2              ", "P6_2      ", "     ",       PRIMITIVE, 21 }, /* 465 */
  {172, "C6^5  ", "P 64            ", "P 6_4                          ", "P 6_4              ", "P6_4      ", "     ",       PRIMITIVE, 21 }, /* 466 */
  {173, "C6^6  ", "P 6c            ", "P 6_3                          ", "P 6_3              ", "P6_3      ", "     ",       PRIMITIVE, 21 }, /* 467 */
  {174, "C3h^1 ", "P -6            ", "P -6                           ", "P -6               ", "P-6       ", "     ",       PRIMITIVE, 22 }, /* 468 */
  {175, "C6h^1 ", "-P 6            ", "P 6/m                          ", "P 6/m              ", "P6/m      ", "     ",       PRIMITIVE, 23 }, /* 469 */
  {176, "C6h^2 ", "-P 6c           ", "P 6_3/m                        ", "P 6_3/m            ", "P6_3/m    ", "     ",       PRIMITIVE, 23 }, /* 470 */
  {177, "D6^1  ", "P 6 2           ", "P 6 2 2                        ", "P 6 2 2            ", "P622      ", "     ",       PRIMITIVE, 24 }, /* 471 */
  {178, "D6^2  ", "P 61 2 (0 0 -1) ", "P 6_1 2 2                      ", "P 6_1 2 2          ", "P6_122    ", "     ",       PRIMITIVE, 24 }, /* 472 */
  {179, "D6^3  ", "P 65 2 (0 0 1)  ", "P 6_5 2 2                      ", "P 6_5 2 2          ", "P6_522    ", "     ",       PRIMITIVE, 24 }, /* 473 */
  {180, "D6^4  ", "P 62 2c (0 0 1) ", "P 6_2 2 2                      ", "P 6_2 2 2          ", "P6_222    ", "     ",       PRIMITIVE, 24 }, /* 474 */
  {181, "D6^5  ", "P 64 2c (0 0 -1)", "P 6_4 2 2                      ", "P 6_4 2 2          ", "P6_422    ", "     ",       PRIMITIVE, 24 }, /* 475 */
  {182, "D6^6  ", "P 6c 2c         ", "P 6_3 2 2                      ", "P 6_3 2 2          ", "P6_322    ", "     ",       PRIMITIVE, 24 }, /* 476 */
  {183, "C6v^1 ", "P 6 -2          ", "P 6 m m                        ", "P 6 m m            ", "P6mm      ", "     ",       PRIMITIVE, 25 }, /* 477 */
  {184, "C6v^2 ", "P 6 -2c         ", "P 6 c c                        ", "P 6 c c            ", "P6cc      ", "     ",       PRIMITIVE, 25 }, /* 478 */
  {185, "C6v^3 ", "P 6c -2         ", "P 6_3 c m                      ", "P 6_3 c m          ", "P6_3cm    ", "     ",       PRIMITIVE, 25 }, /* 479 */
  {186, "C6v^4 ", "P 6c -2c        ", "P 6_3 m c                      ", "P 6_3 m c          ", "P6_3mc    ", "     ",       PRIMITIVE, 25 }, /* 480 */
  {187, "D3h^1 ", "P -6 2          ", "P -6 m 2                       ", "P -6 m 2           ", "P-6m2     ", "     ",       PRIMITIVE, 26 }, /* 481 */
  {188, "D3h^2 ", "P -6c 2         ", "P -6 c 2                       ", "P -6 c 2           ", "P-6c2     ", "     ",       PRIMITIVE, 26 }, /* 482 */
  {189, "D3h^3 ", "P -6 -2         ", "P -6 2 m                       ", "P -6 2 m           ", "P-62m     ", "     ",       PRIMITIVE, 26 }, /* 483 */
  {190, "D3h^4 ", "P -6c -2c       ", "P -6 2 c                       ", "P -6 2 c           ", "P-62c     ", "     ",       PRIMITIVE, 26 }, /* 484 */
  {191, "D6h^1 ", "-P 6 2          ", "P 6/m m m                      ", "P 6/m 2/m 2/m      ", "P6/mmm    ", "     ",       PRIMITIVE, 27 }, /* 485 */
  {192, "D6h^2 ", "-P 6 2c         ", "P 6/m c c                      ", "P 6/m 2/c 2/c      ", "P6/mcc    ", "     ",       PRIMITIVE, 27 }, /* 486 */
  {193, "D6h^3 ", "-P 6c 2         ", "P 6_3/m c m                    ", "P 6_3/m 2/c 2/m    ", "P6_3/mcm  ", "     ",       PRIMITIVE, 27 }, /* 487 */
  {194, "D6h^4 ", "-P 6c 2c        ", "P 6_3/m m c                    ", "P 6_3/m 2/m 2/c    ", "P6_3/mmc  ", "     ",       PRIMITIVE, 27 }, /* 488 */
  {195, "T^1   ", "P 2 2 3         ", "P 2 3                          ", "P 2 3              ", "P23       ", "     ",       PRIMITIVE, 28 }, /* 489 */
  {196, "T^2   ", "F 2 2 3         ", "F 2 3                          ", "F 2 3              ", "F23       ", "     ",            FACE, 28 }, /* 490 */
  {197, "T^3   ", "I 2 2 3         ", "I 2 3                          ", "I 2 3              ", "I23       ", "     ",            BODY, 28 }, /* 491 */
  {198, "T^4   ", "P 2ac 2ab 3     ", "P 2_1 3                        ", "P 2_1 3            ", "P2_13     ", "     ",       PRIMITIVE, 28 }, /* 492 */
  {199, "T^5   ", "I 2b 2c 3       ", "I 2_1 3                        ", "I 2_1 3            ", "I2_13     ", "     ",            BODY, 28 }, /* 493 */
  {200, "Th^1  ", "-P 2 2 3        ", "P m 3                          ", "P 2/m -3           ", "Pm3       ", "     ",       PRIMITIVE, 29 }, /* 494 */
  {201, "Th^2  ", "P 2 2 3 -1n     ", "P n 3                          ", "P 2/n -3           ", "Pn3       ", "1    ",       PRIMITIVE, 29 }, /* 495 */
  {201, "Th^2  ", "-P 2ab 2bc 3    ", "P n 3                          ", "P 2/n -3           ", "Pn3       ", "2    ",       PRIMITIVE, 29 }, /* 496 */
  {202, "Th^3  ", "-F 2 2 3        ", "F m 3                          ", "F 2/m -3           ", "Fm3       ", "     ",            FACE, 29 }, /* 497 */
  {203, "Th^4  ", "F 2 2 3 -1d     ", "F d 3                          ", "F 2/d -3           ", "Fd3       ", "1    ",            FACE, 29 }, /* 498 */
  {203, "Th^4  ", "-F 2uv 2vw 3    ", "F d 3                          ", "F 2/d -3           ", "Fd3       ", "2    ",            FACE, 29 }, /* 499 */
  {204, "Th^5  ", "-I 2 2 3        ", "I m 3                          ", "I 2/m -3           ", "Im3       ", "     ",            BODY, 29 }, /* 500 */
  {205, "Th^6  ", "-P 2ac 2ab 3    ", "P a 3                          ", "P 2_1/a -3         ", "Pa3       ", "     ",       PRIMITIVE, 29 }, /* 501 */
  {206, "Th^7  ", "-I 2b 2c 3      ", "I a 3                          ", "I 2_1/a -3         ", "Ia3       ", "     ",            BODY, 29 }, /* 502 */
  {207, "O^1   ", "P 4 2 3         ", "P 4 3 2                        ", "P 4 3 2            ", "P432      ", "     ",       PRIMITIVE, 30 }, /* 503 */
  {208, "O^2   ", "P 4n 2 3        ", "P 4_2 3 2                      ", "P 4_2 3 2          ", "P4_232    ", "     ",       PRIMITIVE, 30 }, /* 504 */
  {209, "O^3   ", "F 4 2 3         ", "F 4 3 2                        ", "F 4 3 2            ", "F432      ", "     ",            FACE, 30 }, /* 505 */
  {210, "O^4   ", "F 4d 2 3        ", "F 4_1 3 2                      ", "F 4_1 3 2          ", "F4_132    ", "     ",            FACE, 30 }, /* 506 */
  {211, "O^5   ", "I 4 2 3         ", "I 4 3 2                        ", "I 4 3 2            ", "I432      ", "     ",            BODY, 30 }, /* 507 */
  {212, "O^6   ", "P 4acd 2ab 3    ", "P 4_3 3 2                      ", "P 4_3 3 2          ", "P4_332    ", "     ",       PRIMITIVE, 30 }, /* 508 */
  {213, "O^7   ", "P 4bd 2ab 3     ", "P 4_1 3 2                      ", "P 4_1 3 2          ", "P4_132    ", "     ",       PRIMITIVE, 30 }, /* 509 */
  {214, "O^8   ", "I 4bd 2c 3      ", "I 4_1 3 2                      ", "I 4_1 3 2          ", "I4_132    ", "     ",            BODY, 30 }, /* 510 */
  {215, "Td^1  ", "P -4 2 3        ", "P -4 3 m                       ", "P -4 3 m           ", "P-43m     ", "     ",       PRIMITIVE, 31 }, /* 511 */
  {216, "Td^2  ", "F -4 2 3        ", "F -4 3 m                       ", "F -4 3 m           ", "F-43m     ", "     ",            FACE, 31 }, /* 512 */
  {217, "Td^3  ", "I -4 2 3        ", "I -4 3 m                       ", "I -4 3 m           ", "I-43m     ", "     ",            BODY, 31 }, /* 513 */
  {218, "Td^4  ", "P -4n 2 3       ", "P -4 3 n                       ", "P -4 3 n           ", "P-43n     ", "     ",       PRIMITIVE, 31 }, /* 514 */
  {219, "Td^5  ", "F -4c 2 3       ", "F -4 3 c                       ", "F -4 3 c           ", "F-43c     ", "     ",            FACE, 31 }, /* 515 */
  {220, "Td^6  ", "I -4bd 2c 3     ", "I -4 3 d                       ", "I -4 3 d           ", "I-43d     ", "     ",            BODY, 31 }, /* 516 */
  {221, "Oh^1  ", "-P 4 2 3        ", "P m -3 m                       ", "P 4/m -3 2/m       ", "Pm-3m     ", "     ",       PRIMITIVE, 32 }, /* 517 */
  {222, "Oh^2  ", "P 4 2 3 -1n     ", "P n -3 n                       ", "P 4/n -3 2/n       ", "Pn-3n     ", "1    ",       PRIMITIVE, 32 }, /* 518 */
  {222, "Oh^2  ", "-P 4a 2bc 3     ", "P n -3 n                       ", "P 4/n -3 2/n       ", "Pn-3n     ", "2    ",       PRIMITIVE, 32 }, /* 519 */
  {223, "Oh^3  ", "-P 4n 2 3       ", "P m -3 n                       ", "P 4_2/m -3 2/n     ", "Pm-3n     ", "     ",       PRIMITIVE, 32 }, /* 520 */
  {224, "Oh^4  ", "P 4n 2 3 -1n    ", "P n -3 m                       ", "P 4_2/n -3 2/m     ", "Pn-3m     ", "1    ",       PRIMITIVE, 32 }, /* 521 */
  {224, "Oh^4  ", "-P 4bc 2bc 3    ", "P n -3 m                       ", "P 4_2/n -3 2/m     ", "Pn-3m     ", "2    ",       PRIMITIVE, 32 }, /* 522 */
  {225, "Oh^5  ", "-F 4 2 3        ", "F m -3 m                       ", "F 4/m -3 2/m       ", "Fm-3m     ", "     ",            FACE, 32 }, /* 523 */
  {226, "Oh^6  ", "-F 4c 2 3       ", "F m -3 c                       ", "F 4/m -3 2/c       ", "Fm-3c     ", "     ",            FACE, 32 }, /* 524 */
  {227, "Oh^7  ", "F 4d 2 3 -1d    ", "F d -3 m                       ", "F 4_1/d -3 2/m     ", "Fd-3m     ", "1    ",            FACE, 32 }, /* 525 */
  {227, "Oh^7  ", "-F 4vw 2vw 3    ", "F d -3 m                       ", "F 4_1/d -3 2/m     ", "Fd-3m     ", "2    ",            FACE, 32 }, /* 526 */
  {228, "Oh^8  ", "F 4d 2 3 -1cd   ", "F d -3 c                       ", "F 4_1/d -3 2/c     ", "Fd-3c     ", "1    ",            FACE, 32 }, /* 527 */
  {228, "Oh^8  ", "-F 4cvw 2vw 3   ", "F d -3 c                       ", "F 4_1/d -3 2/c     ", "Fd-3c     ", "2    ",            FACE, 32 }, /* 528 */
  {229, "Oh^9  ", "-I 4 2 3        ", "I m -3 m                       ", "I 4/m -3 2/m       ", "Im-3m     ", "     ",            BODY, 32 }, /* 529 */
  {230, "Oh^10 ", "-I 4bd 2c 3     ", "I a -3 d                       ", "I 4_1/a -3 2/d     ", "Ia-3d     ", "     ",            BODY, 32 }, /* 530 */
};

static const int symmetry_operations[] = { 
  0       ,  /* dummy */
  16484   ,  /*    1 (  1) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*    2 (  2) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*    3 (  2) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*    4 (  3) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*    5 (  3) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*    6 (  4) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*    7 (  4) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*    8 (  5) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*    9 (  5) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*   10 (  6) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /*   11 (  6) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*   12 (  7) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*   13 (  7) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*   14 (  8) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022432,  /*   15 (  8) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  16484   ,  /*   16 (  9) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*   17 (  9) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /*   18 (  9) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /*   19 (  9) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*   20 ( 10) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*   21 ( 10) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*   22 ( 10) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /*   23 ( 10) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*   24 ( 11) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*   25 ( 11) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /*   26 ( 11) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /*   27 ( 11) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*   28 ( 12) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*   29 ( 12) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /*   30 ( 12) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /*   31 ( 12) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*   32 ( 13) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*   33 ( 13) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /*   34 ( 13) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /*   35 ( 13) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*   36 ( 14) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*   37 ( 14) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /*   38 ( 14) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /*   39 ( 14) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*   40 ( 15) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*   41 ( 15) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17140694,  /*   42 ( 15) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /*   43 ( 15) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*   44 ( 16) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*   45 ( 16) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /*   46 ( 16) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /*   47 ( 16) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*   48 ( 17) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*   49 ( 17) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /*   50 ( 17) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /*   51 ( 17) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*   52 ( 18) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*   53 ( 18) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*   54 ( 19) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*   55 ( 19) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*   56 ( 20) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*   57 ( 20) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*   58 ( 21) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134420  ,  /*   59 ( 21) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*   60 ( 22) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140532,  /*   61 ( 22) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*   62 ( 23) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022434,  /*   63 ( 23) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*   64 ( 24) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*   65 ( 24) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16484   ,  /*   66 ( 25) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439770,  /*   67 ( 25) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*   68 ( 26) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*   69 ( 26) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*   70 ( 27) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420538 ,  /*   71 ( 27) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /*   72 ( 28) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538636 ,  /*   73 ( 28) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*   74 ( 29) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*   75 ( 29) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*   76 ( 30) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*   77 ( 30) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*   78 ( 30) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /*   79 ( 30) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*   80 ( 31) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*   81 ( 31) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /*   82 ( 31) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /*   83 ( 31) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*   84 ( 32) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*   85 ( 32) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /*   86 ( 32) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /*   87 ( 32) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*   88 ( 33) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*   89 ( 33) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*   90 ( 33) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*   91 ( 33) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*   92 ( 34) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*   93 ( 34) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140694,  /*   94 ( 34) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*   95 ( 34) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*   96 ( 35) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*   97 ( 35) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /*   98 ( 35) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /*   99 ( 35) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  100 ( 36) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  101 ( 36) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /*  102 ( 36) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127572,  /*  103 ( 36) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  104 ( 37) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  105 ( 37) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*  106 ( 37) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426650,  /*  107 ( 37) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  108 ( 38) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  109 ( 38) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /*  110 ( 38) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544748,  /*  111 ( 38) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  112 ( 39) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134420  ,  /*  113 ( 39) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /*  114 ( 39) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18557708,  /*  115 ( 39) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  116 ( 40) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140532,  /*  117 ( 40) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /*  118 ( 40) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18439610,  /*  119 ( 40) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  120 ( 41) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022434,  /*  121 ( 41) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18557870,  /*  122 ( 41) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1551596 ,  /*  123 ( 41) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  124 ( 42) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022434,  /*  125 ( 42) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /*  126 ( 42) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18557708,  /*  127 ( 42) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  128 ( 43) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551596 ,  /*  129 ( 43) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /*  130 ( 43) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17140532,  /*  131 ( 43) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  132 ( 44) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134420  ,  /*  133 ( 44) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /*  134 ( 44) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18439610,  /*  135 ( 44) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  136 ( 45) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  137 ( 45) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  1551758 ,  /*  138 ( 45) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18557868,  /*  139 ( 45) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  140 ( 46) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /*  141 ( 46) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17140694,  /*  142 ( 46) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18439770,  /*  143 ( 46) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  144 ( 47) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  145 ( 47) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  18557870,  /*  146 ( 47) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17140692,  /*  147 ( 47) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  148 ( 48) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  149 ( 48) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17140694,  /*  150 ( 48) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18557868,  /*  151 ( 48) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  152 ( 49) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140692,  /*  153 ( 49) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  1551758 ,  /*  154 ( 49) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18439770,  /*  155 ( 49) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  156 ( 50) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  157 ( 50) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18557870,  /*  158 ( 50) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1551756 ,  /*  159 ( 50) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  160 ( 51) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420538 ,  /*  161 ( 51) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  17140694,  /*  162 ( 51) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544748,  /*  163 ( 51) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  164 ( 52) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538636 ,  /*  165 ( 52) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /*  166 ( 52) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17127572,  /*  167 ( 52) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  168 ( 53) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*  169 ( 53) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /*  170 ( 53) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18426650,  /*  171 ( 53) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  172 ( 54) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*  173 ( 54) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /*  174 ( 54) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544748,  /*  175 ( 54) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  176 ( 55) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538636 ,  /*  177 ( 55) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /*  178 ( 55) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18426650,  /*  179 ( 55) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  180 ( 56) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420538 ,  /*  181 ( 56) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18557870,  /*  182 ( 56) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127572,  /*  183 ( 56) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  184 ( 57) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  185 ( 57) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  186 ( 57) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  187 ( 57) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  188 ( 58) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  189 ( 58) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /*  190 ( 58) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  191 ( 58) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  192 ( 59) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  193 ( 59) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  194 ( 59) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  195 ( 59) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  196 ( 60) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  197 ( 60) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420536 ,  /*  198 ( 60) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /*  199 ( 60) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /*  200 ( 61) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  201 ( 61) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /*  202 ( 61) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /*  203 ( 61) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16484   ,  /*  204 ( 62) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  205 ( 62) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17022432,  /*  206 ( 62) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /*  207 ( 62) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  208 ( 63) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  209 ( 63) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  210 ( 63) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  211 ( 63) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*  212 ( 63) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /*  213 ( 63) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /*  214 ( 63) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /*  215 ( 63) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  216 ( 64) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  217 ( 64) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  218 ( 64) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  219 ( 64) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /*  220 ( 64) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /*  221 ( 64) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /*  222 ( 64) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /*  223 ( 64) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  224 ( 65) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  225 ( 65) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  226 ( 65) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  227 ( 65) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /*  228 ( 65) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  229 ( 65) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544746,  /*  230 ( 65) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /*  231 ( 65) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  232 ( 66) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  233 ( 66) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /*  234 ( 66) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  235 ( 66) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*  236 ( 66) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /*  237 ( 66) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /*  238 ( 66) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  239 ( 66) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  240 ( 67) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  241 ( 67) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /*  242 ( 67) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  243 ( 67) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140694,  /*  244 ( 67) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /*  245 ( 67) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /*  246 ( 67) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  247 ( 67) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  248 ( 68) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  249 ( 68) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /*  250 ( 68) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  251 ( 68) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /*  252 ( 68) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  253 ( 68) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /*  254 ( 68) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /*  255 ( 68) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  256 ( 69) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  257 ( 69) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  258 ( 69) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  259 ( 69) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /*  260 ( 69) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /*  261 ( 69) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /*  262 ( 69) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /*  263 ( 69) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  264 ( 70) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  265 ( 70) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  266 ( 70) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  267 ( 70) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*  268 ( 70) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /*  269 ( 70) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /*  270 ( 70) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /*  271 ( 70) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  272 ( 71) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  273 ( 71) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  274 ( 71) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  275 ( 71) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /*  276 ( 71) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  277 ( 71) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /*  278 ( 71) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /*  279 ( 71) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  280 ( 72) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  281 ( 72) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121458  ,  /*  282 ( 72) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /*  283 ( 72) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*  284 ( 73) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  285 ( 73) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127570,  /*  286 ( 73) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /*  287 ( 73) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  288 ( 74) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  289 ( 74) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009472,  /*  290 ( 74) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /*  291 ( 74) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  292 ( 75) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  293 ( 75) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /*  294 ( 75) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /*  295 ( 75) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16484   ,  /*  296 ( 76) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  297 ( 76) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /*  298 ( 76) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /*  299 ( 76) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  300 ( 77) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  301 ( 77) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /*  302 ( 77) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /*  303 ( 77) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  304 ( 78) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  305 ( 78) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433496 ,  /*  306 ( 78) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  307 ( 78) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /*  308 ( 79) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  309 ( 79) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1551594 ,  /*  310 ( 79) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /*  311 ( 79) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  312 ( 80) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  313 ( 80) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /*  314 ( 80) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /*  315 ( 80) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*  316 ( 81) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  317 ( 81) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538634 ,  /*  318 ( 81) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /*  319 ( 81) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  320 ( 82) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  321 ( 82) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544746,  /*  322 ( 82) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /*  323 ( 82) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  324 ( 83) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  325 ( 83) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426648,  /*  326 ( 83) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /*  327 ( 83) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  328 ( 84) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  329 ( 84) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /*  330 ( 84) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  331 ( 84) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  332 ( 85) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  333 ( 85) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /*  334 ( 85) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /*  335 ( 85) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  336 ( 86) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  337 ( 86) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /*  338 ( 86) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  339 ( 86) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  340 ( 87) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  341 ( 87) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18439608,  /*  342 ( 87) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /*  343 ( 87) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  344 ( 88) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  345 ( 88) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18557706,  /*  346 ( 88) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /*  347 ( 88) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  348 ( 89) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  349 ( 89) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17140530,  /*  350 ( 89) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /*  351 ( 89) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  352 ( 90) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  353 ( 90) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121458  ,  /*  354 ( 90) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /*  355 ( 90) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /*  356 ( 90) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /*  357 ( 90) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18544746,  /*  358 ( 90) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /*  359 ( 90) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  360 ( 91) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  361 ( 91) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127570,  /*  362 ( 91) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /*  363 ( 91) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /*  364 ( 91) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /*  365 ( 91) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18426648,  /*  366 ( 91) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /*  367 ( 91) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  368 ( 92) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  369 ( 92) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009472,  /*  370 ( 92) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /*  371 ( 92) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18557870,  /*  372 ( 92) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  373 ( 92) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1538634 ,  /*  374 ( 92) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /*  375 ( 92) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  376 ( 93) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  377 ( 93) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009472,  /*  378 ( 93) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /*  379 ( 93) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /*  380 ( 93) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /*  381 ( 93) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18544746,  /*  382 ( 93) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /*  383 ( 93) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  384 ( 94) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  385 ( 94) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538634 ,  /*  386 ( 94) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /*  387 ( 94) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /*  388 ( 94) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /*  389 ( 94) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17127570,  /*  390 ( 94) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /*  391 ( 94) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  392 ( 95) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  393 ( 95) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121458  ,  /*  394 ( 95) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /*  395 ( 95) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /*  396 ( 95) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  397 ( 95) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18426648,  /*  398 ( 95) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /*  399 ( 95) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  400 ( 96) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  401 ( 96) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /*  402 ( 96) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /*  403 ( 96) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  1551758 ,  /*  404 ( 96) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /*  405 ( 96) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18544586,  /*  406 ( 96) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /*  407 ( 96) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  408 ( 97) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  409 ( 97) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /*  410 ( 97) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  411 ( 97) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17140694,  /*  412 ( 97) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /*  413 ( 97) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /*  414 ( 97) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /*  415 ( 97) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  416 ( 98) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  417 ( 98) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /*  418 ( 98) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /*  419 ( 98) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  18557870,  /*  420 ( 98) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  421 ( 98) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17127410,  /*  422 ( 98) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  423 ( 98) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  424 ( 99) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  425 ( 99) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /*  426 ( 99) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /*  427 ( 99) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17140694,  /*  428 ( 99) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /*  429 ( 99) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18544586,  /*  430 ( 99) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /*  431 ( 99) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  432 (100) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  433 (100) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /*  434 (100) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  435 (100) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  1551758 ,  /*  436 (100) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /*  437 (100) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18426488,  /*  438 (100) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /*  439 (100) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  440 (101) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  441 (101) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /*  442 (101) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /*  443 (101) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18557870,  /*  444 (101) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  445 (101) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1538474 ,  /*  446 (101) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  447 (101) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  448 (102) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  449 (102) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433496 ,  /*  450 (102) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  451 (102) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  17140694,  /*  452 (102) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /*  453 (102) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18557706,  /*  454 (102) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /*  455 (102) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  456 (103) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  457 (103) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1551594 ,  /*  458 (103) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /*  459 (103) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /*  460 (103) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /*  461 (103) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17140530,  /*  462 (103) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /*  463 (103) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  464 (104) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  465 (104) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /*  466 (104) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /*  467 (104) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /*  468 (104) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  469 (104) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /*  470 (104) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /*  471 (104) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  472 (105) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  473 (105) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /*  474 (105) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /*  475 (105) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /*  476 (105) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /*  477 (105) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18557706,  /*  478 (105) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /*  479 (105) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  480 (106) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  481 (106) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1551594 ,  /*  482 (106) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /*  483 (106) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /*  484 (106) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /*  485 (106) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /*  486 (106) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /*  487 (106) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  488 (107) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /*  489 (107) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433496 ,  /*  490 (107) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  491 (107) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18557870,  /*  492 (107) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /*  493 (107) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17140530,  /*  494 (107) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /*  495 (107) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  496 (108) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  497 (108) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*  498 (108) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  499 (108) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  500 (109) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  501 (109) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16320   ,  /*  502 (109) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121458  ,  /*  503 (109) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16484   ,  /*  504 (110) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17009312,  /*  505 (110) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022432,  /*  506 (110) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  3360    ,  /*  507 (110) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  508 (111) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  509 (111) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1433496 ,  /*  510 (111) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420536 ,  /*  511 (111) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  512 (112) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  513 (112) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439608,  /*  514 (112) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /*  515 (112) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  516 (113) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538474 ,  /*  517 (113) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16320   ,  /*  518 (113) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538634 ,  /*  519 (113) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  520 (114) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17127410,  /*  521 (114) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /*  522 (114) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  3360    ,  /*  523 (114) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  524 (115) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17127410,  /*  525 (115) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439608,  /*  526 (115) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  1538634 ,  /*  527 (115) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  528 (116) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  529 (116) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16320   ,  /*  530 (116) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121458  ,  /*  531 (116) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  18439772,  /*  532 (116) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544586,  /*  533 (116) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18439608,  /*  534 (116) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18544746,  /*  535 (116) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  536 (117) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17009312,  /*  537 (117) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022432,  /*  538 (117) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  3360    ,  /*  539 (117) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*  540 (117) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18544586,  /*  541 (117) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /*  542 (117) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1538634 ,  /*  543 (117) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  544 (118) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  545 (118) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1433496 ,  /*  546 (118) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420536 ,  /*  547 (118) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17140694,  /*  548 (118) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /*  549 (118) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18557706,  /*  550 (118) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544746,  /*  551 (118) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  552 (119) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  553 (119) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*  554 (119) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  555 (119) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /*  556 (119) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /*  557 (119) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /*  558 (119) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /*  559 (119) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  560 (120) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  561 (120) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*  562 (120) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  563 (120) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*  564 (120) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /*  565 (120) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /*  566 (120) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /*  567 (120) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  568 (121) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  569 (121) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*  570 (121) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  571 (121) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140694,  /*  572 (121) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /*  573 (121) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /*  574 (121) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /*  575 (121) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  576 (122) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  577 (122) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*  578 (122) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  579 (122) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*  580 (122) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /*  581 (122) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /*  582 (122) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /*  583 (122) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17140694,  /*  584 (122) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /*  585 (122) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /*  586 (122) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /*  587 (122) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439772,  /*  588 (122) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /*  589 (122) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /*  590 (122) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /*  591 (122) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  592 (123) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  593 (123) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /*  594 (123) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /*  595 (123) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /*  596 (123) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /*  597 (123) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /*  598 (123) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544746,  /*  599 (123) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  600 (124) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420376 ,  /*  601 (124) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  134418  ,  /*  602 (124) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  1538634 ,  /*  603 (124) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18557870,  /*  604 (124) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127410,  /*  605 (124) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439608,  /*  606 (124) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17009472,  /*  607 (124) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16484   ,  /*  608 (125) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  609 (125) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  610 (125) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*  611 (125) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  612 (126) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  613 (126) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  614 (126) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  615 (126) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  616 (127) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  617 (127) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  618 (127) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*  619 (127) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  620 (128) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  621 (128) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  3362    ,  /*  622 (128) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134420  ,  /*  623 (128) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*  624 (129) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  625 (129) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  121460  ,  /*  626 (129) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16322   ,  /*  627 (129) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  628 (130) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  629 (130) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022432,  /*  630 (130) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  16322   ,  /*  631 (130) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  632 (131) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  633 (131) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17022432,  /*  634 (131) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /*  635 (131) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  636 (132) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  637 (132) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1420538 ,  /*  638 (132) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1420536 ,  /*  639 (132) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  640 (133) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  641 (133) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  3362    ,  /*  642 (133) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /*  643 (133) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  644 (134) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  645 (134) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*  646 (134) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /*  647 (134) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*  648 (135) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  649 (135) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16320   ,  /*  650 (135) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17022434,  /*  651 (135) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  652 (136) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  653 (136) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  654 (136) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /*  655 (136) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  656 (137) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  657 (137) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17009474,  /*  658 (137) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17022434,  /*  659 (137) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  660 (138) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  661 (138) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1420538 ,  /*  662 (138) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1433498 ,  /*  663 (138) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /*  664 (139) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  665 (139) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /*  666 (139) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433498 ,  /*  667 (139) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /*  668 (140) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /*  669 (140) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /*  670 (140) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134420  ,  /*  671 (140) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*  672 (141) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /*  673 (141) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /*  674 (141) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /*  675 (141) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  676 (142) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  677 (142) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /*  678 (142) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  3360    ,  /*  679 (142) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  680 (143) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  681 (143) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17127572,  /*  682 (143) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17022434,  /*  683 (143) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  684 (144) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  685 (144) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  1420538 ,  /*  686 (144) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1551596 ,  /*  687 (144) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  688 (145) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  689 (145) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17022432,  /*  690 (145) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  18439610,  /*  691 (145) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  692 (146) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140692,  /*  693 (146) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17022432,  /*  694 (146) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  134420  ,  /*  695 (146) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /*  696 (147) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /*  697 (147) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  121460  ,  /*  698 (147) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1420536 ,  /*  699 (147) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  700 (148) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  701 (148) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18426650,  /*  702 (148) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  1420536 ,  /*  703 (148) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  704 (149) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  705 (149) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1538636 ,  /*  706 (149) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /*  707 (149) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  708 (150) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  709 (150) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17127572,  /*  710 (150) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /*  711 (150) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  712 (151) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140692,  /*  713 (151) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16320   ,  /*  714 (151) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17140532,  /*  715 (151) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  716 (152) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439770,  /*  717 (152) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16320   ,  /*  718 (152) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18439610,  /*  719 (152) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  720 (153) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439770,  /*  721 (153) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /*  722 (153) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  3360    ,  /*  723 (153) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  724 (154) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /*  725 (154) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /*  726 (154) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  3360    ,  /*  727 (154) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  728 (155) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17127410,  /*  729 (155) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  3362    ,  /*  730 (155) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140532,  /*  731 (155) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  732 (156) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538474 ,  /*  733 (156) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1538636 ,  /*  734 (156) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  16322   ,  /*  735 (156) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  736 (157) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439770,  /*  737 (157) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /*  738 (157) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  16322   ,  /*  739 (157) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /*  740 (158) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  741 (158) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140530,  /*  742 (158) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /*  743 (158) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  744 (159) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  745 (159) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538636 ,  /*  746 (159) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /*  747 (159) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  748 (160) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18439770,  /*  749 (160) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  3362    ,  /*  750 (160) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18426648,  /*  751 (160) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  752 (161) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  753 (161) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18426650,  /*  754 (161) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /*  755 (161) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  756 (162) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /*  757 (162) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16320   ,  /*  758 (162) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1551596 ,  /*  759 (162) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  760 (163) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140692,  /*  761 (163) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /*  762 (163) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  3360    ,  /*  763 (163) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  764 (164) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  765 (164) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18544748,  /*  766 (164) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18439610,  /*  767 (164) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  768 (165) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  769 (165) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18426650,  /*  770 (165) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18557708,  /*  771 (165) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  772 (166) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /*  773 (166) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17022432,  /*  774 (166) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  18557708,  /*  775 (166) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  776 (167) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18557868,  /*  777 (167) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17022432,  /*  778 (167) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  1551596 ,  /*  779 (167) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  780 (168) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18557868,  /*  781 (168) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17127572,  /*  782 (168) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  1420536 ,  /*  783 (168) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  784 (169) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140692,  /*  785 (169) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18544748,  /*  786 (169) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1420536 ,  /*  787 (169) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /*  788 (170) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  789 (170) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18544748,  /*  790 (170) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /*  791 (170) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  792 (171) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18557868,  /*  793 (171) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16320   ,  /*  794 (171) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18557708,  /*  795 (171) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  796 (172) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18557868,  /*  797 (172) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /*  798 (172) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  3360    ,  /*  799 (172) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /*  800 (173) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  801 (173) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  802 (173) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*  803 (173) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*  804 (173) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /*  805 (173) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18426650,  /*  806 (173) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /*  807 (173) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  808 (174) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  809 (174) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  810 (174) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  811 (174) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /*  812 (174) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  813 (174) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /*  814 (174) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /*  815 (174) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  816 (175) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  817 (175) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  818 (175) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*  819 (175) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140694,  /*  820 (175) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  821 (175) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /*  822 (175) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /*  823 (175) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  824 (176) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  825 (176) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  3362    ,  /*  826 (176) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134420  ,  /*  827 (176) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /*  828 (176) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544586,  /*  829 (176) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18426650,  /*  830 (176) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18557708,  /*  831 (176) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  832 (177) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /*  833 (177) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  121460  ,  /*  834 (177) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16322   ,  /*  835 (177) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*  836 (177) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544586,  /*  837 (177) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18544748,  /*  838 (177) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18439610,  /*  839 (177) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  840 (178) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  841 (178) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022432,  /*  842 (178) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  16322   ,  /*  843 (178) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /*  844 (178) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18557868,  /*  845 (178) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /*  846 (178) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1551596 ,  /*  847 (178) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  848 (179) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  849 (179) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17022432,  /*  850 (179) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /*  851 (179) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /*  852 (179) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  853 (179) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18557706,  /*  854 (179) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /*  855 (179) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  856 (180) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  857 (180) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1420538 ,  /*  858 (180) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1420536 ,  /*  859 (180) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17140694,  /*  860 (180) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  861 (180) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18544748,  /*  862 (180) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /*  863 (180) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  864 (181) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  865 (181) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  3362    ,  /*  866 (181) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /*  867 (181) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17140694,  /*  868 (181) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18557868,  /*  869 (181) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17127572,  /*  870 (181) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544746,  /*  871 (181) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /*  872 (182) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  873 (182) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*  874 (182) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /*  875 (182) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /*  876 (182) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /*  877 (182) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18544748,  /*  878 (182) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /*  879 (182) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  880 (183) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /*  881 (183) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16320   ,  /*  882 (183) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17022434,  /*  883 (183) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /*  884 (183) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18557868,  /*  885 (183) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  1551594 ,  /*  886 (183) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18557708,  /*  887 (183) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /*  888 (184) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  889 (184) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  890 (184) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /*  891 (184) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140694,  /*  892 (184) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18557868,  /*  893 (184) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /*  894 (184) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /*  895 (184) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /*  896 (185) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  897 (185) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  898 (185) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*  899 (185) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /*  900 (185) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /*  901 (185) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1538636 ,  /*  902 (185) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /*  903 (185) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /*  904 (186) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  905 (186) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /*  906 (186) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /*  907 (186) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /*  908 (186) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /*  909 (186) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17127572,  /*  910 (186) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /*  911 (186) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  912 (187) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  913 (187) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  914 (187) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  915 (187) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /*  916 (187) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /*  917 (187) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /*  918 (187) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /*  919 (187) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /*  920 (188) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  921 (188) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /*  922 (188) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /*  923 (188) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /*  924 (188) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /*  925 (188) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /*  926 (188) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /*  927 (188) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /*  928 (189) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  929 (189) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  930 (189) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*  931 (189) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /*  932 (189) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /*  933 (189) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /*  934 (189) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /*  935 (189) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  936 (190) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /*  937 (190) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /*  938 (190) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /*  939 (190) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*  940 (190) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /*  941 (190) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /*  942 (190) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /*  943 (190) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  944 (191) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  945 (191) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*  946 (191) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /*  947 (191) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  1551758 ,  /*  948 (191) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /*  949 (191) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1420538 ,  /*  950 (191) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1433498 ,  /*  951 (191) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /*  952 (192) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  953 (192) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /*  954 (192) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /*  955 (192) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17140694,  /*  956 (192) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /*  957 (192) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17009474,  /*  958 (192) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17022434,  /*  959 (192) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  960 (193) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /*  961 (193) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /*  962 (193) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134420  ,  /*  963 (193) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17140694,  /*  964 (193) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17022594,  /*  965 (193) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17140530,  /*  966 (193) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17022434,  /*  967 (193) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  968 (194) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  969 (194) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /*  970 (194) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433498 ,  /*  971 (194) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18439772,  /*  972 (194) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17022594,  /*  973 (194) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18439608,  /*  974 (194) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17022434,  /*  975 (194) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /*  976 (195) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /*  977 (195) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  978 (195) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /*  979 (195) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /*  980 (195) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17022594,  /*  981 (195) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /*  982 (195) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  18426648,  /*  983 (195) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /*  984 (196) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /*  985 (196) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /*  986 (196) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /*  987 (196) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /*  988 (196) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1433658 ,  /*  989 (196) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /*  990 (196) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1538634 ,  /*  991 (196) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /*  992 (197) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /*  993 (197) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17009474,  /*  994 (197) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17022434,  /*  995 (197) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /*  996 (197) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /*  997 (197) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18544748,  /*  998 (197) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /*  999 (197) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1000 (198) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1001 (198) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1420538 ,  /* 1002 (198) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1433498 ,  /* 1003 (198) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  17140694,  /* 1004 (198) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 1005 (198) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18544748,  /* 1006 (198) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 1007 (198) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1008 (199) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /* 1009 (199) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /* 1010 (199) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433498 ,  /* 1011 (199) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  17140694,  /* 1012 (199) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18557868,  /* 1013 (199) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17140530,  /* 1014 (199) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18557708,  /* 1015 (199) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1016 (200) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /* 1017 (200) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 1018 (200) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134420  ,  /* 1019 (200) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /* 1020 (200) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18557868,  /* 1021 (200) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 1022 (200) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18557708,  /* 1023 (200) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1024 (201) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /* 1025 (201) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1026 (201) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 1027 (201) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /* 1028 (201) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18557868,  /* 1029 (201) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1030 (201) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 1031 (201) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /* 1032 (202) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /* 1033 (202) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1034 (202) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  3360    ,  /* 1035 (202) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /* 1036 (202) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18557868,  /* 1037 (202) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1038 (202) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1538634 ,  /* 1039 (202) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /* 1040 (203) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1041 (203) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17127572,  /* 1042 (203) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 1043 (203) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /* 1044 (203) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 1045 (203) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18426650,  /* 1046 (203) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 1047 (203) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1048 (204) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1049 (204) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1538636 ,  /* 1050 (204) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 1051 (204) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 1052 (204) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 1053 (204) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18426650,  /* 1054 (204) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 1055 (204) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1056 (205) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /* 1057 (205) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16320   ,  /* 1058 (205) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1551596 ,  /* 1059 (205) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 1060 (205) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18439770,  /* 1061 (205) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17140530,  /* 1062 (205) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18439610,  /* 1063 (205) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1064 (206) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /* 1065 (206) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16320   ,  /* 1066 (206) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1551596 ,  /* 1067 (206) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /* 1068 (206) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17140692,  /* 1069 (206) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 1070 (206) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17140532,  /* 1071 (206) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1072 (207) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1551756 ,  /* 1073 (207) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1074 (207) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  3360    ,  /* 1075 (207) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439772,  /* 1076 (207) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17140692,  /* 1077 (207) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1078 (207) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18426648,  /* 1079 (207) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /* 1080 (208) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140692,  /* 1081 (208) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1082 (208) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  3360    ,  /* 1083 (208) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /* 1084 (208) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18439770,  /* 1085 (208) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1086 (208) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  1538634 ,  /* 1087 (208) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16484   ,  /* 1088 (209) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1089 (209) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /* 1090 (209) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /* 1091 (209) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 1092 (209) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 1093 (209) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1538636 ,  /* 1094 (209) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 1095 (209) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 1096 (209) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 1097 (209) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17127572,  /* 1098 (209) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 1099 (209) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439772,  /* 1100 (209) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 1101 (209) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18426650,  /* 1102 (209) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 1103 (209) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1104 (210) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1105 (210) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 1106 (210) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1107 (210) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 1108 (210) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1109 (210) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 1110 (210) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1111 (210) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 1112 (210) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1113 (210) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 1114 (210) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1115 (210) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439772,  /* 1116 (210) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1117 (210) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 1118 (210) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1119 (210) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1120 (211) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1121 (211) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1122 (211) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 1123 (211) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551758 ,  /* 1124 (211) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1125 (211) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1126 (211) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 1127 (211) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17140694,  /* 1128 (211) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1129 (211) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1130 (211) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 1131 (211) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439772,  /* 1132 (211) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1133 (211) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1134 (211) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 1135 (211) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1551758 ,  /* 1136 (212) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 1137 (212) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  9274055 ,  /* 1138 (212) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 3, 3] */
  9287015 ,  /* 1139 (212) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 3] */
  16484   ,  /* 1140 (212) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1141 (212) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  10809329,  /* 1142 (212) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 9, 9] */
  10822289,  /* 1143 (212) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 9] */
  18439772,  /* 1144 (212) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 1145 (212) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  26398265,  /* 1146 (212) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 3, 9] */
  26411225,  /* 1147 (212) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 9] */
  17140694,  /* 1148 (212) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 1149 (212) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  27697343,  /* 1150 (212) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 9, 3] */
  27710303,  /* 1151 (212) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 3] */
  18439772,  /* 1152 (213) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  9287175 ,  /* 1153 (213) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 3] */
  18439608,  /* 1154 (213) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  10822289,  /* 1155 (213) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 9] */
  17140694,  /* 1156 (213) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  10822449,  /* 1157 (213) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 9] */
  17140530,  /* 1158 (213) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  9287015 ,  /* 1159 (213) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 3] */
  1551758 ,  /* 1160 (213) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  26411385,  /* 1161 (213) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 9] */
  1551594 ,  /* 1162 (213) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  27710303,  /* 1163 (213) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 3] */
  16484   ,  /* 1164 (213) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  27710463,  /* 1165 (213) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 3] */
  16320   ,  /* 1166 (213) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  26411225,  /* 1167 (213) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 9] */
  17140694,  /* 1168 (214) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  10822449,  /* 1169 (214) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 9] */
  27697343,  /* 1170 (214) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 9, 3] */
  1538634 ,  /* 1171 (214) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18439772,  /* 1172 (214) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  9287175 ,  /* 1173 (214) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 3] */
  26398265,  /* 1174 (214) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 3, 9] */
  3360    ,  /* 1175 (214) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 1176 (214) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  27710463,  /* 1177 (214) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 3] */
  10809329,  /* 1178 (214) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 9, 9] */
  18426648,  /* 1179 (214) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1551758 ,  /* 1180 (214) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  26411385,  /* 1181 (214) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 9] */
  9274055 ,  /* 1182 (214) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 3, 3] */
  17127570,  /* 1183 (214) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /* 1184 (215) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1185 (215) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /* 1186 (215) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /* 1187 (215) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 1188 (215) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 1189 (215) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18544748,  /* 1190 (215) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 1191 (215) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1192 (216) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1193 (216) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 1194 (216) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1195 (216) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 1196 (216) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1197 (216) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 1198 (216) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1199 (216) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1200 (217) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1201 (217) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1202 (217) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 1203 (217) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 1204 (217) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1205 (217) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1206 (217) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 1207 (217) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 1208 (218) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1209 (218) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /* 1210 (218) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /* 1211 (218) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /* 1212 (218) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 1213 (218) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18426650,  /* 1214 (218) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 1215 (218) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1216 (219) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /* 1217 (219) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16320   ,  /* 1218 (219) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17022434,  /* 1219 (219) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18557870,  /* 1220 (219) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1551756 ,  /* 1221 (219) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18557706,  /* 1222 (219) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1551596 ,  /* 1223 (219) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1224 (220) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /* 1225 (220) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1226 (220) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /* 1227 (220) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 1228 (220) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17140692,  /* 1229 (220) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1230 (220) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544746,  /* 1231 (220) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 1232 (221) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1233 (221) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17009474,  /* 1234 (221) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17022434,  /* 1235 (221) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18557870,  /* 1236 (221) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 1237 (221) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  1538636 ,  /* 1238 (221) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 1239 (221) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1240 (222) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1241 (222) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1420538 ,  /* 1242 (222) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1433498 ,  /* 1243 (222) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18557870,  /* 1244 (222) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 1245 (222) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17127572,  /* 1246 (222) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 1247 (222) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1248 (223) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1433658 ,  /* 1249 (223) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /* 1250 (223) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1433498 ,  /* 1251 (223) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18557870,  /* 1252 (223) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17140692,  /* 1253 (223) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18557706,  /* 1254 (223) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17140532,  /* 1255 (223) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1256 (224) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /* 1257 (224) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 1258 (224) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  134420  ,  /* 1259 (224) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /* 1260 (224) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18439770,  /* 1261 (224) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18557706,  /* 1262 (224) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18439610,  /* 1263 (224) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1264 (225) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /* 1265 (225) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1266 (225) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 1267 (225) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 1268 (225) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18439770,  /* 1269 (225) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1270 (225) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544746,  /* 1271 (225) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 1272 (226) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17022594,  /* 1273 (226) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1274 (226) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  3360    ,  /* 1275 (226) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 1276 (226) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1551756 ,  /* 1277 (226) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1278 (226) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18544746,  /* 1279 (226) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 1280 (227) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1281 (227) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1282 (227) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1283 (227) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 1284 (227) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1285 (227) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 1286 (227) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1287 (227) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1288 (228) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1289 (228) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 1290 (228) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 1291 (228) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18544584,  /* 1292 (228) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18557868,  /* 1293 (228) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1294 (228) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 1295 (228) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1296 (229) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1297 (229) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1298 (229) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1299 (229) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1551594 ,  /* 1300 (229) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1301 (229) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  17127570,  /* 1302 (229) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1303 (229) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1304 (230) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1305 (230) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1306 (230) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1307 (230) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 1308 (230) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1309 (230) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 1310 (230) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1311 (230) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1312 (231) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1313 (231) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1314 (231) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1315 (231) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16320   ,  /* 1316 (231) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1317 (231) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17009472,  /* 1318 (231) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1319 (231) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1320 (232) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1321 (232) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1322 (232) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1323 (232) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433496 ,  /* 1324 (232) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1325 (232) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /* 1326 (232) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1327 (232) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1328 (233) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1329 (233) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 1330 (233) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 1331 (233) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18426486,  /* 1332 (233) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 1333 (233) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1334 (233) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 1335 (233) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1336 (234) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1337 (234) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1338 (234) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1339 (234) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1433496 ,  /* 1340 (234) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1341 (234) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  17009472,  /* 1342 (234) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1343 (234) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1344 (235) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1345 (235) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 1346 (235) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 1347 (235) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538472 ,  /* 1348 (235) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551756 ,  /* 1349 (235) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1350 (235) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 1351 (235) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1352 (236) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1353 (236) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1354 (236) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1355 (236) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 1356 (236) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1357 (236) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  121458  ,  /* 1358 (236) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1359 (236) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1360 (237) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1361 (237) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 1362 (237) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 1363 (237) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17127408,  /* 1364 (237) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140692,  /* 1365 (237) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1366 (237) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 1367 (237) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1368 (238) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1369 (238) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1370 (238) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1371 (238) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  134418  ,  /* 1372 (238) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1373 (238) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  17127570,  /* 1374 (238) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1375 (238) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1376 (239) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1377 (239) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1378 (239) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1379 (239) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022432,  /* 1380 (239) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1381 (239) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  3360    ,  /* 1382 (239) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1383 (239) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1384 (240) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1385 (240) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1386 (240) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1387 (240) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /* 1388 (240) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1389 (240) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /* 1390 (240) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1391 (240) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1392 (241) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1393 (241) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1394 (241) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1395 (241) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1433496 ,  /* 1396 (241) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1397 (241) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1420536 ,  /* 1398 (241) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1399 (241) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1400 (242) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1401 (242) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1402 (242) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1403 (242) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134418  ,  /* 1404 (242) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1405 (242) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 1406 (242) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1407 (242) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1408 (243) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1409 (243) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1410 (243) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1411 (243) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 1412 (243) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1413 (243) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121458  ,  /* 1414 (243) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1415 (243) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1416 (244) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1417 (244) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1418 (244) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1419 (244) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17022432,  /* 1420 (244) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1421 (244) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17009472,  /* 1422 (244) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1423 (244) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1424 (245) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1425 (245) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1426 (245) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1427 (245) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  1551594 ,  /* 1428 (245) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1429 (245) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18544746,  /* 1430 (245) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1431 (245) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1432 (246) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1433 (246) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1434 (246) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1435 (246) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  18557706,  /* 1436 (246) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1437 (246) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /* 1438 (246) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1439 (246) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1440 (247) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1441 (247) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1442 (247) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1443 (247) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  1433496 ,  /* 1444 (247) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1445 (247) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  17127570,  /* 1446 (247) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1447 (247) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1448 (248) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1449 (248) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1450 (248) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1451 (248) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  134418  ,  /* 1452 (248) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1453 (248) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18544746,  /* 1454 (248) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1455 (248) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1456 (249) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1457 (249) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1458 (249) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1459 (249) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18557706,  /* 1460 (249) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1461 (249) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  121458  ,  /* 1462 (249) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1463 (249) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1464 (250) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1465 (250) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1466 (250) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1467 (250) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  1551594 ,  /* 1468 (250) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1469 (250) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  17009472,  /* 1470 (250) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1471 (250) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1472 (251) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1473 (251) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1474 (251) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1475 (251) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  16320   ,  /* 1476 (251) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1477 (251) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17127570,  /* 1478 (251) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1479 (251) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1480 (252) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1481 (252) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1482 (252) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1483 (252) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 1484 (252) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1485 (252) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  3360    ,  /* 1486 (252) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1487 (252) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1488 (253) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1489 (253) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1490 (253) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1491 (253) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 1492 (253) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1493 (253) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  3360    ,  /* 1494 (253) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1495 (253) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1496 (254) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1497 (254) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1498 (254) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1499 (254) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140530,  /* 1500 (254) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1501 (254) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 1502 (254) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1503 (254) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1504 (255) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1505 (255) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1506 (255) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1507 (255) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551594 ,  /* 1508 (255) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1509 (255) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 1510 (255) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1511 (255) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1512 (256) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1513 (256) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1514 (256) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1515 (256) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  16320   ,  /* 1516 (256) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1517 (256) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18426648,  /* 1518 (256) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1519 (256) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1520 (257) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1521 (257) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1522 (257) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1523 (257) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17140530,  /* 1524 (257) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1525 (257) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  121458  ,  /* 1526 (257) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1527 (257) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1528 (258) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1529 (258) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1530 (258) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1531 (258) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  134418  ,  /* 1532 (258) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1533 (258) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1538634 ,  /* 1534 (258) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1535 (258) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1536 (259) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1537 (259) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1538 (259) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1539 (259) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  1433496 ,  /* 1540 (259) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1541 (259) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18426648,  /* 1542 (259) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1543 (259) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1544 (260) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1545 (260) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1546 (260) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1547 (260) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  134418  ,  /* 1548 (260) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1549 (260) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  17009472,  /* 1550 (260) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1551 (260) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1552 (261) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1553 (261) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1554 (261) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1555 (261) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1433496 ,  /* 1556 (261) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1557 (261) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  121458  ,  /* 1558 (261) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1559 (261) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1560 (262) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1561 (262) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1562 (262) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1563 (262) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  18439608,  /* 1564 (262) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1565 (262) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17009472,  /* 1566 (262) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1567 (262) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1568 (263) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1569 (263) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1570 (263) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1571 (263) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18439608,  /* 1572 (263) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1573 (263) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 1574 (263) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1575 (263) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1576 (264) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1577 (264) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1578 (264) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1579 (264) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16320   ,  /* 1580 (264) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1581 (264) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538634 ,  /* 1582 (264) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1583 (264) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1584 (265) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1585 (265) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1586 (265) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1587 (265) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 1588 (265) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1589 (265) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  3360    ,  /* 1590 (265) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1591 (265) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1592 (266) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1593 (266) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1594 (266) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1595 (266) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17140530,  /* 1596 (266) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1597 (266) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  1538634 ,  /* 1598 (266) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1599 (266) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1600 (267) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1601 (267) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1602 (267) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1603 (267) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  1551594 ,  /* 1604 (267) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1605 (267) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18426648,  /* 1606 (267) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1607 (267) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1608 (268) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1609 (268) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1610 (268) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1611 (268) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18439608,  /* 1612 (268) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1613 (268) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17127570,  /* 1614 (268) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1615 (268) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1616 (269) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1617 (269) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1618 (269) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1619 (269) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  1433496 ,  /* 1620 (269) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1621 (269) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1538634 ,  /* 1622 (269) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1623 (269) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1624 (270) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1625 (270) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1626 (270) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1627 (270) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  17140530,  /* 1628 (270) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1629 (270) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17009472,  /* 1630 (270) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1631 (270) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1632 (271) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1633 (271) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1634 (271) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1635 (271) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17022432,  /* 1636 (271) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1637 (271) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  121458  ,  /* 1638 (271) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1639 (271) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1640 (272) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1641 (272) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1642 (272) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1643 (272) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17022432,  /* 1644 (272) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1645 (272) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  18426648,  /* 1646 (272) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1647 (272) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1648 (273) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1649 (273) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1650 (273) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1651 (273) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18439608,  /* 1652 (273) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1653 (273) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  1420536 ,  /* 1654 (273) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1655 (273) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1656 (274) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1657 (274) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1658 (274) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1659 (274) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  134418  ,  /* 1660 (274) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1661 (274) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1420536 ,  /* 1662 (274) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1663 (274) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1664 (275) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1665 (275) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1666 (275) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1667 (275) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18557706,  /* 1668 (275) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1669 (275) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 1670 (275) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1671 (275) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1672 (276) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1673 (276) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1674 (276) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1675 (276) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  16320   ,  /* 1676 (276) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1677 (276) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18544746,  /* 1678 (276) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1679 (276) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1680 (277) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1681 (277) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1682 (277) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1683 (277) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 1684 (277) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1685 (277) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  3360    ,  /* 1686 (277) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1687 (277) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1688 (278) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 1689 (278) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439608,  /* 1690 (278) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 1691 (278) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426486,  /* 1692 (278) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 1693 (278) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  3362    ,  /* 1694 (278) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /* 1695 (278) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1696 (279) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1697 (279) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 1698 (279) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1699 (279) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17022432,  /* 1700 (279) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1701 (279) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  1420536 ,  /* 1702 (279) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1703 (279) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1704 (280) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538474 ,  /* 1705 (280) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16320   ,  /* 1706 (280) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538634 ,  /* 1707 (280) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538472 ,  /* 1708 (280) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  16482   ,  /* 1709 (280) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538636 ,  /* 1710 (280) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  16322   ,  /* 1711 (280) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 1712 (281) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1713 (281) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1714 (281) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1715 (281) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  1551594 ,  /* 1716 (281) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1717 (281) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1420536 ,  /* 1718 (281) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1719 (281) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1720 (282) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17127410,  /* 1721 (282) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 1722 (282) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  3360    ,  /* 1723 (282) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17127408,  /* 1724 (282) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  16482   ,  /* 1725 (282) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1726 (282) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17140532,  /* 1727 (282) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1728 (283) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1729 (283) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1730 (283) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1731 (283) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  17022432,  /* 1732 (283) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1733 (283) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17127570,  /* 1734 (283) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1735 (283) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1736 (284) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1737 (284) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1738 (284) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1739 (284) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 1740 (284) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1741 (284) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  121458  ,  /* 1742 (284) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1743 (284) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 1744 (285) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1745 (285) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1746 (285) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1747 (285) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  134418  ,  /* 1748 (285) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1749 (285) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18426648,  /* 1750 (285) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1751 (285) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1752 (286) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1753 (286) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1754 (286) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1755 (286) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18557706,  /* 1756 (286) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1757 (286) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1538634 ,  /* 1758 (286) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1759 (286) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1760 (287) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1761 (287) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1762 (287) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1763 (287) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18557706,  /* 1764 (287) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1765 (287) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17009472,  /* 1766 (287) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1767 (287) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 1768 (288) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1769 (288) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1770 (288) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1771 (288) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  1433496 ,  /* 1772 (288) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1773 (288) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18544746,  /* 1774 (288) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1775 (288) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1776 (289) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1777 (289) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1778 (289) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1779 (289) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  17140530,  /* 1780 (289) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1781 (289) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544746,  /* 1782 (289) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1783 (289) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1784 (290) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1785 (290) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1786 (290) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1787 (290) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 1788 (290) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1789 (290) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  1538634 ,  /* 1790 (290) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1791 (290) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1792 (291) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1793 (291) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1794 (291) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1795 (291) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17140530,  /* 1796 (291) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1797 (291) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18426648,  /* 1798 (291) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1799 (291) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1800 (292) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1801 (292) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1802 (292) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1803 (292) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18557706,  /* 1804 (292) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1805 (292) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1420536 ,  /* 1806 (292) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1807 (292) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1808 (293) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1809 (293) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1810 (293) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1811 (293) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17022432,  /* 1812 (293) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1813 (293) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  18544746,  /* 1814 (293) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1815 (293) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1816 (294) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1817 (294) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1818 (294) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1819 (294) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  18439608,  /* 1820 (294) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1821 (294) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544746,  /* 1822 (294) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1823 (294) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1824 (295) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1825 (295) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1826 (295) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1827 (295) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17140530,  /* 1828 (295) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1829 (295) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  1420536 ,  /* 1830 (295) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1831 (295) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 1832 (296) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1833 (296) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18544586,  /* 1834 (296) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1835 (296) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17022432,  /* 1836 (296) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1837 (296) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  1538634 ,  /* 1838 (296) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1839 (296) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1840 (297) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1841 (297) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1842 (297) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1843 (297) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  18557706,  /* 1844 (297) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1845 (297) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 1846 (297) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1847 (297) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1848 (298) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1849 (298) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1850 (298) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1851 (298) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 1852 (298) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1853 (298) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121458  ,  /* 1854 (298) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 1855 (298) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /* 1856 (298) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 1857 (298) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18544586,  /* 1858 (298) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1859 (298) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 1860 (298) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1861 (298) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544746,  /* 1862 (298) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1863 (298) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1864 (299) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1865 (299) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 1866 (299) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 1867 (299) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134418  ,  /* 1868 (299) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 1869 (299) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 1870 (299) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1871 (299) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /* 1872 (299) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 1873 (299) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18544586,  /* 1874 (299) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1875 (299) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 1876 (299) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1877 (299) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 1878 (299) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1879 (299) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1880 (300) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1881 (300) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 1882 (300) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 1883 (300) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022432,  /* 1884 (300) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1885 (300) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  3360    ,  /* 1886 (300) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1887 (300) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 1888 (300) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 1889 (300) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18544586,  /* 1890 (300) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1891 (300) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 1892 (300) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1893 (300) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1538634 ,  /* 1894 (300) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1895 (300) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1896 (301) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1897 (301) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1898 (301) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1899 (301) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17022432,  /* 1900 (301) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 1901 (301) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17009472,  /* 1902 (301) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 1903 (301) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /* 1904 (301) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 1905 (301) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 1906 (301) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1907 (301) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18557706,  /* 1908 (301) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1909 (301) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 1910 (301) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1911 (301) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1912 (302) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1913 (302) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1914 (302) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1915 (302) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1433496 ,  /* 1916 (302) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 1917 (302) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1420536 ,  /* 1918 (302) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1919 (302) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  17140694,  /* 1920 (302) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 1921 (302) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 1922 (302) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1923 (302) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18557706,  /* 1924 (302) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 1925 (302) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 1926 (302) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1927 (302) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1928 (303) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1929 (303) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 1930 (303) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 1931 (303) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /* 1932 (303) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1933 (303) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /* 1934 (303) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 1935 (303) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  17140694,  /* 1936 (303) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 1937 (303) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18544586,  /* 1938 (303) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 1939 (303) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17140530,  /* 1940 (303) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1941 (303) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544746,  /* 1942 (303) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 1943 (303) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 1944 (304) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1945 (304) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1946 (304) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1947 (304) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16320   ,  /* 1948 (304) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 1949 (304) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538634 ,  /* 1950 (304) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1951 (304) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /* 1952 (304) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 1953 (304) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17127410,  /* 1954 (304) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1955 (304) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 1956 (304) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1957 (304) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17127570,  /* 1958 (304) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1959 (304) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 1960 (305) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1961 (305) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 1962 (305) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 1963 (305) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 1964 (305) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 1965 (305) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  3360    ,  /* 1966 (305) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1967 (305) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /* 1968 (305) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 1969 (305) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17127410,  /* 1970 (305) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1971 (305) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 1972 (305) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1973 (305) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18426648,  /* 1974 (305) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 1975 (305) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 1976 (306) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1977 (306) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 1978 (306) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 1979 (306) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 1980 (306) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1981 (306) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  3360    ,  /* 1982 (306) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 1983 (306) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 1984 (306) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 1985 (306) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18426488,  /* 1986 (306) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 1987 (306) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 1988 (306) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 1989 (306) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  1538634 ,  /* 1990 (306) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 1991 (306) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 1992 (307) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 1993 (307) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 1994 (307) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 1995 (307) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17140530,  /* 1996 (307) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 1997 (307) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 1998 (307) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 1999 (307) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /* 2000 (307) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2001 (307) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 2002 (307) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2003 (307) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18439608,  /* 2004 (307) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2005 (307) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 2006 (307) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2007 (307) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2008 (308) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2009 (308) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2010 (308) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2011 (308) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1551594 ,  /* 2012 (308) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2013 (308) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 2014 (308) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2015 (308) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 2016 (308) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2017 (308) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 2018 (308) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2019 (308) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 2020 (308) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2021 (308) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 2022 (308) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2023 (308) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2024 (309) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2025 (309) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 2026 (309) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2027 (309) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  16320   ,  /* 2028 (309) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2029 (309) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1538634 ,  /* 2030 (309) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2031 (309) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 2032 (309) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2033 (309) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /* 2034 (309) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2035 (309) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17140530,  /* 2036 (309) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2037 (309) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18426648,  /* 2038 (309) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2039 (309) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2040 (310) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2041 (310) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2042 (310) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2043 (310) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 2044 (310) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2045 (310) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 2046 (310) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2047 (310) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /* 2048 (310) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2049 (310) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426488,  /* 2050 (310) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2051 (310) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 2052 (310) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2053 (310) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 2054 (310) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2055 (310) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2056 (311) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2057 (311) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2058 (311) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2059 (311) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 2060 (311) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2061 (311) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 2062 (311) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2063 (311) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 2064 (311) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2065 (311) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 2066 (311) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2067 (311) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 2068 (311) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2069 (311) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 2070 (311) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2071 (311) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 2072 (312) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2073 (312) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2074 (312) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2075 (312) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 2076 (312) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2077 (312) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 2078 (312) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2079 (312) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /* 2080 (312) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2081 (312) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 2082 (312) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2083 (312) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 2084 (312) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2085 (312) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 2086 (312) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2087 (312) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2088 (313) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2089 (313) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2090 (313) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2091 (313) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 2092 (313) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2093 (313) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 2094 (313) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2095 (313) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /* 2096 (313) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2097 (313) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426488,  /* 2098 (313) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2099 (313) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18557706,  /* 2100 (313) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2101 (313) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 2102 (313) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2103 (313) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2104 (314) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2105 (314) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 2106 (314) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2107 (314) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16320   ,  /* 2108 (314) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2109 (314) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17009472,  /* 2110 (314) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2111 (314) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /* 2112 (314) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2113 (314) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18544586,  /* 2114 (314) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2115 (314) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  1551594 ,  /* 2116 (314) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2117 (314) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18544746,  /* 2118 (314) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2119 (314) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2120 (315) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2121 (315) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2122 (315) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2123 (315) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433496 ,  /* 2124 (315) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2125 (315) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /* 2126 (315) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2127 (315) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  17140694,  /* 2128 (315) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2129 (315) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18544586,  /* 2130 (315) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2131 (315) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 2132 (315) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2133 (315) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /* 2134 (315) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2135 (315) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2136 (316) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2137 (316) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2138 (316) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2139 (316) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /* 2140 (316) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2141 (316) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /* 2142 (316) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 2143 (316) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18439772,  /* 2144 (316) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2145 (316) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17009312,  /* 2146 (316) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2147 (316) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18439608,  /* 2148 (316) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2149 (316) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17009472,  /* 2150 (316) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2151 (316) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 2152 (317) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2153 (317) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2154 (317) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2155 (317) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433496 ,  /* 2156 (317) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2157 (317) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /* 2158 (317) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2159 (317) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18439772,  /* 2160 (317) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2161 (317) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17009312,  /* 2162 (317) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2163 (317) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022432,  /* 2164 (317) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 2165 (317) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  18426648,  /* 2166 (317) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2167 (317) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2168 (318) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2169 (318) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 2170 (318) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 2171 (318) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134418  ,  /* 2172 (318) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2173 (318) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 2174 (318) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2175 (318) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 2176 (318) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2177 (318) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1420376 ,  /* 2178 (318) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2179 (318) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433496 ,  /* 2180 (318) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2181 (318) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1538634 ,  /* 2182 (318) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2183 (318) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 2184 (319) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2185 (319) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2186 (319) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2187 (319) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 2188 (319) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2189 (319) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 2190 (319) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2191 (319) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  1551758 ,  /* 2192 (319) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2193 (319) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 2194 (319) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2195 (319) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1433496 ,  /* 2196 (319) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2197 (319) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1420536 ,  /* 2198 (319) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 2199 (319) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  16484   ,  /* 2200 (320) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2201 (320) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2202 (320) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2203 (320) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 2204 (320) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2205 (320) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 2206 (320) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2207 (320) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17140694,  /* 2208 (320) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2209 (320) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 2210 (320) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2211 (320) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17022432,  /* 2212 (320) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 2213 (320) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17009472,  /* 2214 (320) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2215 (320) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 2216 (321) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2217 (321) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 2218 (321) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 2219 (321) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 2220 (321) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2221 (321) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121458  ,  /* 2222 (321) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2223 (321) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17140694,  /* 2224 (321) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2225 (321) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17009312,  /* 2226 (321) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2227 (321) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17140530,  /* 2228 (321) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2229 (321) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17009472,  /* 2230 (321) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2231 (321) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 2232 (322) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2233 (322) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2234 (322) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2235 (322) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538472 ,  /* 2236 (322) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551756 ,  /* 2237 (322) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2238 (322) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 2239 (322) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /* 2240 (322) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 2241 (322) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 2242 (322) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 2243 (322) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17127408,  /* 2244 (322) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140692,  /* 2245 (322) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2246 (322) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 2247 (322) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2248 (323) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2249 (323) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2250 (323) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2251 (323) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 2252 (323) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2253 (323) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  121458  ,  /* 2254 (323) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2255 (323) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18439772,  /* 2256 (323) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2257 (323) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17009312,  /* 2258 (323) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2259 (323) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17140530,  /* 2260 (323) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2261 (323) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544746,  /* 2262 (323) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2263 (323) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2264 (324) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2265 (324) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2266 (324) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2267 (324) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538472 ,  /* 2268 (324) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551756 ,  /* 2269 (324) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2270 (324) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 2271 (324) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /* 2272 (324) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 2273 (324) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 2274 (324) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 2275 (324) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17127408,  /* 2276 (324) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140692,  /* 2277 (324) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2278 (324) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 2279 (324) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2280 (325) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2281 (325) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2282 (325) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2283 (325) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  134418  ,  /* 2284 (325) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2285 (325) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1538634 ,  /* 2286 (325) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2287 (325) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18439772,  /* 2288 (325) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2289 (325) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17009312,  /* 2290 (325) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2291 (325) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18557706,  /* 2292 (325) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2293 (325) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /* 2294 (325) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2295 (325) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2296 (326) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2297 (326) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2298 (326) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2299 (326) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17127408,  /* 2300 (326) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140692,  /* 2301 (326) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2302 (326) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 2303 (326) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /* 2304 (326) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 2305 (326) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /* 2306 (326) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /* 2307 (326) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18426486,  /* 2308 (326) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 2309 (326) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2310 (326) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 2311 (326) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2312 (327) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2313 (327) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 2314 (327) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2315 (327) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  134418  ,  /* 2316 (327) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2317 (327) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  17127570,  /* 2318 (327) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2319 (327) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /* 2320 (327) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2321 (327) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18544586,  /* 2322 (327) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2323 (327) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  1433496 ,  /* 2324 (327) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2325 (327) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18426648,  /* 2326 (327) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2327 (327) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2328 (328) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2329 (328) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2330 (328) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2331 (328) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17127408,  /* 2332 (328) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17140692,  /* 2333 (328) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2334 (328) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17140532,  /* 2335 (328) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  1551758 ,  /* 2336 (328) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 2337 (328) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /* 2338 (328) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /* 2339 (328) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18426486,  /* 2340 (328) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 2341 (328) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2342 (328) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 2343 (328) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2344 (329) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2345 (329) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 2346 (329) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2347 (329) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  134418  ,  /* 2348 (329) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2349 (329) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  17009472,  /* 2350 (329) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2351 (329) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  1551758 ,  /* 2352 (329) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2353 (329) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  18426488,  /* 2354 (329) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2355 (329) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1433496 ,  /* 2356 (329) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2357 (329) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18544746,  /* 2358 (329) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2359 (329) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2360 (330) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2361 (330) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2362 (330) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2363 (330) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538472 ,  /* 2364 (330) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551756 ,  /* 2365 (330) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2366 (330) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 2367 (330) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 2368 (330) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 2369 (330) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 2370 (330) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /* 2371 (330) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18426486,  /* 2372 (330) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 2373 (330) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2374 (330) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 2375 (330) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2376 (331) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2377 (331) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1538474 ,  /* 2378 (331) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2379 (331) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1433496 ,  /* 2380 (331) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2381 (331) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  121458  ,  /* 2382 (331) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2383 (331) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17140694,  /* 2384 (331) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2385 (331) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /* 2386 (331) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2387 (331) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18557706,  /* 2388 (331) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2389 (331) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17009472,  /* 2390 (331) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2391 (331) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 2392 (332) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2393 (332) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2394 (332) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2395 (332) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1538472 ,  /* 2396 (332) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1551756 ,  /* 2397 (332) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2398 (332) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1551596 ,  /* 2399 (332) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 2400 (332) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 2401 (332) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 2402 (332) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /* 2403 (332) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18426486,  /* 2404 (332) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 2405 (332) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2406 (332) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 2407 (332) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2408 (333) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2409 (333) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2410 (333) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2411 (333) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 2412 (333) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2413 (333) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  121458  ,  /* 2414 (333) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2415 (333) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  17140694,  /* 2416 (333) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2417 (333) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  18544586,  /* 2418 (333) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2419 (333) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 2420 (333) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2421 (333) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17009472,  /* 2422 (333) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2423 (333) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 2424 (334) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2425 (334) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2426 (334) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2427 (334) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 2428 (334) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2429 (334) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 2430 (334) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2431 (334) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1551758 ,  /* 2432 (334) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2433 (334) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 2434 (334) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2435 (334) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 2436 (334) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2437 (334) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 2438 (334) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2439 (334) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  17140694,  /* 2440 (334) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2441 (334) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 2442 (334) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2443 (334) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 2444 (334) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2445 (334) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 2446 (334) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2447 (334) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439772,  /* 2448 (334) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2449 (334) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426488,  /* 2450 (334) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2451 (334) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 2452 (334) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2453 (334) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 2454 (334) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2455 (334) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2456 (335) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2457 (335) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2458 (335) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2459 (335) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  9273891 ,  /* 2460 (335) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 3, 3] */
  9287175 ,  /* 2461 (335) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 3] */
  9274055 ,  /* 2462 (335) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 3, 3] */
  9287015 ,  /* 2463 (335) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 3] */
  1551758 ,  /* 2464 (335) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 2465 (335) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /* 2466 (335) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /* 2467 (335) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  10809165,  /* 2468 (335) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 9, 9] */
  10822449,  /* 2469 (335) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 9] */
  10809329,  /* 2470 (335) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 9, 9] */
  10822289,  /* 2471 (335) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 9] */
  17140694,  /* 2472 (335) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 2473 (335) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 2474 (335) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /* 2475 (335) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  26398101,  /* 2476 (335) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 3, 9] */
  26411385,  /* 2477 (335) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 9] */
  26398265,  /* 2478 (335) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 3, 9] */
  26411225,  /* 2479 (335) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 9] */
  18439772,  /* 2480 (335) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 2481 (335) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 2482 (335) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 2483 (335) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  27697179,  /* 2484 (335) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 9, 3] */
  27710463,  /* 2485 (335) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 3] */
  27697343,  /* 2486 (335) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 9, 3] */
  27710303,  /* 2487 (335) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 3] */
  16484   ,  /* 2488 (336) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2489 (336) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  9214844 ,  /* 2490 (336) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 0] */
  9228126 ,  /* 2491 (336) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 0] */
  783957  ,  /* 2492 (336) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 3, 3] */
  770999  ,  /* 2493 (336) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 3, 3] */
  8565465 ,  /* 2494 (336) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 0, 3] */
  8578427 ,  /* 2495 (336) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 0, 3] */
  1551758 ,  /* 2496 (336) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 2497 (336) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  10750118,  /* 2498 (336) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 6] */
  10763400,  /* 2499 (336) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 6] */
  2319231 ,  /* 2500 (336) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 9, 9] */
  2306273 ,  /* 2501 (336) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 9, 9] */
  10100739,  /* 2502 (336) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 6, 9] */
  10113701,  /* 2503 (336) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 6, 9] */
  17140694,  /* 2504 (336) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 2505 (336) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  26339054,  /* 2506 (336) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 6] */
  26352336,  /* 2507 (336) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 6] */
  17908167,  /* 2508 (336) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 3, 9] */
  17895209,  /* 2509 (336) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 3, 9] */
  25689675,  /* 2510 (336) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 0, 9] */
  25702637,  /* 2511 (336) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 0, 9] */
  18439772,  /* 2512 (336) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2513 (336) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  27638132,  /* 2514 (336) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 0] */
  27651414,  /* 2515 (336) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 0] */
  19207245,  /* 2516 (336) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 9, 3] */
  19194287,  /* 2517 (336) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 9, 3] */
  26988753,  /* 2518 (336) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 6, 3] */
  27001715,  /* 2519 (336) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 6, 3] */
  16484   ,  /* 2520 (337) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2521 (337) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2522 (337) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2523 (337) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 2524 (337) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2525 (337) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 2526 (337) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2527 (337) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 2528 (337) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2529 (337) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 2530 (337) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2531 (337) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 2532 (337) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2533 (337) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 2534 (337) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2535 (337) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2536 (338) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2537 (338) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2538 (338) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2539 (338) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 2540 (338) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2541 (338) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 2542 (338) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2543 (338) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /* 2544 (338) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2545 (338) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 2546 (338) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2547 (338) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 2548 (338) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2549 (338) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 2550 (338) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2551 (338) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2552 (339) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2553 (339) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 2554 (339) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2555 (339) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  16320   ,  /* 2556 (339) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2557 (339) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17009472,  /* 2558 (339) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2559 (339) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18557870,  /* 2560 (339) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2561 (339) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1538474 ,  /* 2562 (339) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2563 (339) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  18557706,  /* 2564 (339) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2565 (339) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  1538634 ,  /* 2566 (339) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2567 (339) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 2568 (340) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2569 (340) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2570 (340) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2571 (340) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433496 ,  /* 2572 (340) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2573 (340) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  3360    ,  /* 2574 (340) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2575 (340) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 2576 (340) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2577 (340) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17127410,  /* 2578 (340) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2579 (340) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 2580 (340) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2581 (340) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  18544746,  /* 2582 (340) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2583 (340) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2584 (341) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2585 (341) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2586 (341) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2587 (341) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  134418  ,  /* 2588 (341) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2589 (341) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1538634 ,  /* 2590 (341) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2591 (341) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18557870,  /* 2592 (341) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2593 (341) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17127410,  /* 2594 (341) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2595 (341) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 2596 (341) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2597 (341) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17009472,  /* 2598 (341) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2599 (341) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  16484   ,  /* 2600 (342) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2601 (342) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 2602 (342) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2603 (342) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  1433496 ,  /* 2604 (342) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2605 (342) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18426648,  /* 2606 (342) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2607 (342) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18557870,  /* 2608 (342) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2609 (342) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1538474 ,  /* 2610 (342) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2611 (342) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  17140530,  /* 2612 (342) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2613 (342) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  121458  ,  /* 2614 (342) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2615 (342) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 2616 (343) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2617 (343) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 2618 (343) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2619 (343) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  16320   ,  /* 2620 (343) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2621 (343) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /* 2622 (343) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 2623 (343) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18557870,  /* 2624 (343) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2625 (343) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17127410,  /* 2626 (343) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2627 (343) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18557706,  /* 2628 (343) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2629 (343) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /* 2630 (343) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2631 (343) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2632 (344) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2633 (344) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17009312,  /* 2634 (344) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  17022594,  /* 2635 (344) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022432,  /* 2636 (344) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 2637 (344) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  3360    ,  /* 2638 (344) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2639 (344) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 2640 (344) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2641 (344) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  1538474 ,  /* 2642 (344) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 2643 (344) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 2644 (344) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2645 (344) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18544746,  /* 2646 (344) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2647 (344) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2648 (345) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2649 (345) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 2650 (345) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 2651 (345) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134418  ,  /* 2652 (345) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 2653 (345) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 2654 (345) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 2655 (345) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 2656 (345) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2657 (345) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18426488,  /* 2658 (345) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2659 (345) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 2660 (345) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 2661 (345) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18544746,  /* 2662 (345) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 2663 (345) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2664 (346) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2665 (346) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2666 (346) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2667 (346) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  1433496 ,  /* 2668 (346) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 2669 (346) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  1420536 ,  /* 2670 (346) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 2671 (346) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18557870,  /* 2672 (346) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2673 (346) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 2674 (346) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2675 (346) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  17140530,  /* 2676 (346) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 2677 (346) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 2678 (346) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 2679 (346) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  16484   ,  /* 2680 (347) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2681 (347) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2682 (347) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2683 (347) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  17022432,  /* 2684 (347) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 2685 (347) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  17009472,  /* 2686 (347) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 2687 (347) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18557870,  /* 2688 (347) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2689 (347) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 2690 (347) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2691 (347) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  1551594 ,  /* 2692 (347) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 2693 (347) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 2694 (347) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 2695 (347) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  16484   ,  /* 2696 (348) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2697 (348) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 2698 (348) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 2699 (348) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 2700 (348) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 2701 (348) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  121458  ,  /* 2702 (348) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 2703 (348) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /* 2704 (348) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2705 (348) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18426488,  /* 2706 (348) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2707 (348) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18557706,  /* 2708 (348) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 2709 (348) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 2710 (348) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 2711 (348) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2712 (349) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2713 (349) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2714 (349) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2715 (349) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 2716 (350) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  66947   ,  /* 2717 (350) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 3] */
  121298  ,  /* 2718 (350) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  188933  ,  /* 2719 (350) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 9] */
  16484   ,  /* 2720 (351) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125996  ,  /* 2721 (351) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3200    ,  /* 2722 (351) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  129884  ,  /* 2723 (351) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 2724 (352) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  185045  ,  /* 2725 (352) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 9] */
  121298  ,  /* 2726 (352) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  70835   ,  /* 2727 (352) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 3] */
  16484   ,  /* 2728 (353) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2729 (353) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2730 (353) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2731 (353) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 2732 (353) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18549284,  /* 2733 (353) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 2734 (353) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553172,  /* 2735 (353) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2736 (354) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 2737 (354) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 2738 (354) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 2739 (354) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  18557870,  /* 2740 (354) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 2741 (354) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 2742 (354) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 2743 (354) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  16484   ,  /* 2744 (355) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 2745 (355) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2746 (355) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 2747 (355) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2748 (356) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 2749 (356) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2750 (356) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 2751 (356) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 2752 (356) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 2753 (356) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 2754 (356) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 2755 (356) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 2756 (357) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2757 (357) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 2758 (357) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 2759 (357) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2760 (357) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2761 (357) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 2762 (357) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 2763 (357) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2764 (358) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2765 (358) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  125996  ,  /* 2766 (358) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  129882  ,  /* 2767 (358) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3200    ,  /* 2768 (358) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2769 (358) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  129884  ,  /* 2770 (358) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  125994  ,  /* 2771 (358) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  16484   ,  /* 2772 (359) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18431186,  /* 2773 (359) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  3200    ,  /* 2774 (359) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18435074,  /* 2775 (359) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 2776 (359) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  11784   ,  /* 2777 (359) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18439770,  /* 2778 (359) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  7896    ,  /* 2779 (359) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2780 (360) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2781 (360) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17014010,  /* 2782 (360) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 2783 (360) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  18426488,  /* 2784 (360) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2785 (360) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1428962 ,  /* 2786 (360) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 2787 (360) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  16484   ,  /* 2788 (361) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 2789 (361) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 2790 (361) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 2791 (361) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2792 (361) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  11784   ,  /* 2793 (361) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557868,  /* 2794 (361) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  7896    ,  /* 2795 (361) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2796 (362) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2797 (362) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1543172 ,  /* 2798 (362) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1547058 ,  /* 2799 (362) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  18426488,  /* 2800 (362) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 2801 (362) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17135996,  /* 2802 (362) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17132106,  /* 2803 (362) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  16484   ,  /* 2804 (363) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2805 (363) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 2806 (363) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 2807 (363) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 2808 (363) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 2809 (363) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 2810 (363) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 2811 (363) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 2812 (363) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2813 (363) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18549284,  /* 2814 (363) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 2815 (363) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 2816 (363) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 2817 (363) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18553172,  /* 2818 (363) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 2819 (363) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 2820 (364) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 2821 (364) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 2822 (364) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 2823 (364) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  1479423 ,  /* 2824 (364) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 3] */
  11784   ,  /* 2825 (364) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  17199741,  /* 2826 (364) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 9] */
  18549282,  /* 2827 (364) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557870,  /* 2828 (364) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 2829 (364) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 2830 (364) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 2831 (364) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  17186457,  /* 2832 (364) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 9] */
  18553170,  /* 2833 (364) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  1492707 ,  /* 2834 (364) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 3] */
  7896    ,  /* 2835 (364) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2836 (365) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 2837 (365) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  26284703,  /* 2838 (365) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 3] */
  26288589,  /* 2839 (365) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 3] */
  17127410,  /* 2840 (365) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 2841 (365) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  27823865,  /* 2842 (365) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 9] */
  27819975,  /* 2843 (365) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 9] */
  18557870,  /* 2844 (365) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 2845 (365) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  10813865,  /* 2846 (365) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 9] */
  10817751,  /* 2847 (365) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 9] */
  1420376 ,  /* 2848 (365) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 2849 (365) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  9282479 ,  /* 2850 (365) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 3] */
  9278589 ,  /* 2851 (365) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 3] */
  16484   ,  /* 2852 (366) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2853 (366) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2854 (366) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2855 (366) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2856 (366) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 2857 (366) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2858 (366) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 2859 (366) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2860 (367) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18431186,  /* 2861 (367) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  3200    ,  /* 2862 (367) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18435074,  /* 2863 (367) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 2864 (367) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  7410    ,  /* 2865 (367) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18426648,  /* 2866 (367) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  12270   ,  /* 2867 (367) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2868 (368) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  66947   ,  /* 2869 (368) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 3] */
  121298  ,  /* 2870 (368) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  188933  ,  /* 2871 (368) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 9] */
  134418  ,  /* 2872 (368) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  66459   ,  /* 2873 (368) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 3] */
  3360    ,  /* 2874 (368) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  189417  ,  /* 2875 (368) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 9] */
  16484   ,  /* 2876 (369) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18490235,  /* 2877 (369) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 3] */
  121298  ,  /* 2878 (369) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18612221,  /* 2879 (369) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 9] */
  18616755,  /* 2880 (369) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 9] */
  125508  ,  /* 2881 (369) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18485697,  /* 2882 (369) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 3] */
  12270   ,  /* 2883 (369) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2884 (370) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125996  ,  /* 2885 (370) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3200    ,  /* 2886 (370) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  129884  ,  /* 2887 (370) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16320   ,  /* 2888 (370) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  125508  ,  /* 2889 (370) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3360    ,  /* 2890 (370) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  130368  ,  /* 2891 (370) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  16484   ,  /* 2892 (371) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 2893 (371) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 2894 (371) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 2895 (371) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /* 2896 (371) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  7410    ,  /* 2897 (371) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544746,  /* 2898 (371) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  12270   ,  /* 2899 (371) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2900 (372) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  185045  ,  /* 2901 (372) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 9] */
  121298  ,  /* 2902 (372) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  70835   ,  /* 2903 (372) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 3] */
  134418  ,  /* 2904 (372) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  184557  ,  /* 2905 (372) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 9] */
  3360    ,  /* 2906 (372) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  71319   ,  /* 2907 (372) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 3] */
  16484   ,  /* 2908 (373) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18608333,  /* 2909 (373) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 9] */
  121298  ,  /* 2910 (373) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18494123,  /* 2911 (373) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 3] */
  18498657,  /* 2912 (373) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 3] */
  125508  ,  /* 2913 (373) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18603795,  /* 2914 (373) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 9] */
  12270   ,  /* 2915 (373) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2916 (374) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2917 (374) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2918 (374) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2919 (374) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 2920 (374) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 2921 (374) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 2922 (374) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 2923 (374) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 2924 (374) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18549284,  /* 2925 (374) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 2926 (374) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553172,  /* 2927 (374) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /* 2928 (374) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18548796,  /* 2929 (374) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544746,  /* 2930 (374) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18553656,  /* 2931 (374) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 2932 (375) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 2933 (375) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 2934 (375) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 2935 (375) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  1492545 ,  /* 2936 (375) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 3] */
  7410    ,  /* 2937 (375) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  17186619,  /* 2938 (375) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 9] */
  18553656,  /* 2939 (375) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557870,  /* 2940 (375) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 2941 (375) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 2942 (375) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 2943 (375) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  17199579,  /* 2944 (375) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 9] */
  18548796,  /* 2945 (375) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  1479585 ,  /* 2946 (375) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 3] */
  12270   ,  /* 2947 (375) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 2948 (376) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2949 (376) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2950 (376) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2951 (376) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /* 2952 (376) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  12272   ,  /* 2953 (376) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /* 2954 (376) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7412    ,  /* 2955 (376) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 2956 (377) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2957 (377) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2958 (377) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2959 (377) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18426650,  /* 2960 (377) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18435560,  /* 2961 (377) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 2962 (377) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18430700,  /* 2963 (377) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 2964 (378) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125996  ,  /* 2965 (378) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3200    ,  /* 2966 (378) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  129884  ,  /* 2967 (378) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  121460  ,  /* 2968 (378) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  12272   ,  /* 2969 (378) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  134420  ,  /* 2970 (378) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  7412    ,  /* 2971 (378) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 2972 (379) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 2973 (379) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 2974 (379) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 2975 (379) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544748,  /* 2976 (379) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  12272   ,  /* 2977 (379) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557708,  /* 2978 (379) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  7412    ,  /* 2979 (379) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 2980 (380) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2981 (380) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2982 (380) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2983 (380) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /* 2984 (380) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  130370  ,  /* 2985 (380) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /* 2986 (380) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  125510  ,  /* 2987 (380) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 2988 (381) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 2989 (381) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 2990 (381) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 2991 (381) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18544748,  /* 2992 (381) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553658,  /* 2993 (381) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 2994 (381) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18548798,  /* 2995 (381) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 2996 (382) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125996  ,  /* 2997 (382) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3200    ,  /* 2998 (382) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  129884  ,  /* 2999 (382) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3362    ,  /* 3000 (382) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  130370  ,  /* 3001 (382) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16322   ,  /* 3002 (382) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  125510  ,  /* 3003 (382) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3004 (383) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125996  ,  /* 3005 (383) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3200    ,  /* 3006 (383) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  129884  ,  /* 3007 (383) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  18426650,  /* 3008 (383) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18553658,  /* 3009 (383) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18439610,  /* 3010 (383) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18548798,  /* 3011 (383) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3012 (384) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 3013 (384) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 3014 (384) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 3015 (384) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3362    ,  /* 3016 (384) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  12272   ,  /* 3017 (384) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16322   ,  /* 3018 (384) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7412    ,  /* 3019 (384) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 3020 (384) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18549284,  /* 3021 (384) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 3022 (384) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553172,  /* 3023 (384) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544748,  /* 3024 (384) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553658,  /* 3025 (384) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 3026 (384) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18548798,  /* 3027 (384) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3028 (385) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 3029 (385) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 3030 (385) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 3031 (385) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  121460  ,  /* 3032 (385) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  130370  ,  /* 3033 (385) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  134420  ,  /* 3034 (385) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  125510  ,  /* 3035 (385) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /* 3036 (385) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18549284,  /* 3037 (385) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 3038 (385) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553172,  /* 3039 (385) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18426650,  /* 3040 (385) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18435560,  /* 3041 (385) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 3042 (385) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18430700,  /* 3043 (385) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3044 (386) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 3045 (386) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 3046 (386) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 3047 (386) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3362    ,  /* 3048 (386) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1488497 ,  /* 3049 (386) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18557708,  /* 3050 (386) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17190671,  /* 3051 (386) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  18557870,  /* 3052 (386) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 3053 (386) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 3054 (386) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 3055 (386) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544748,  /* 3056 (386) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17195531,  /* 3057 (386) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  16322   ,  /* 3058 (386) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1483637 ,  /* 3059 (386) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  16484   ,  /* 3060 (387) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 3061 (387) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 3062 (387) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 3063 (387) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  121460  ,  /* 3064 (387) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1606595 ,  /* 3065 (387) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 9] */
  18439610,  /* 3066 (387) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  17072573,  /* 3067 (387) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 3] */
  18557870,  /* 3068 (387) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 3069 (387) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 3070 (387) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 3071 (387) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18426650,  /* 3072 (387) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17077433,  /* 3073 (387) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 3] */
  134420  ,  /* 3074 (387) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  1601735 ,  /* 3075 (387) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 9] */
  16484   ,  /* 3076 (388) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3077 (388) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3078 (388) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3079 (388) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 3080 (388) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3081 (388) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 3082 (388) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3083 (388) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3084 (389) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3085 (389) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3086 (389) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3087 (389) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 3088 (389) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3089 (389) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 3090 (389) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 3091 (389) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3092 (390) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3093 (390) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3094 (390) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3095 (390) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18439608,  /* 3096 (390) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 3097 (390) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 3098 (390) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 3099 (390) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3100 (391) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3101 (391) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3102 (391) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3103 (391) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557706,  /* 3104 (391) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3105 (391) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 3106 (391) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3107 (391) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3108 (392) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3109 (392) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3110 (392) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3111 (392) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3112 (392) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 3113 (392) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 3114 (392) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 3115 (392) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 3116 (393) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3117 (393) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3118 (393) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3119 (393) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  121460  ,  /* 3120 (393) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  125508  ,  /* 3121 (393) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 3122 (393) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  130368  ,  /* 3123 (393) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  16484   ,  /* 3124 (394) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3125 (394) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3126 (394) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3127 (394) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18426650,  /* 3128 (394) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18430698,  /* 3129 (394) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 3130 (394) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18435558,  /* 3131 (394) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /* 3132 (395) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3133 (395) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3134 (395) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3135 (395) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544748,  /* 3136 (395) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18548796,  /* 3137 (395) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 3138 (395) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553656,  /* 3139 (395) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 3140 (396) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3141 (396) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3142 (396) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3143 (396) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3144 (396) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 3145 (396) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 3146 (396) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 3147 (396) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557870,  /* 3148 (396) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3149 (396) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 3150 (396) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3151 (396) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 3152 (396) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18548796,  /* 3153 (396) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 3154 (396) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553656,  /* 3155 (396) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16484   ,  /* 3156 (397) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3157 (397) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3158 (397) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3159 (397) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  121460  ,  /* 3160 (397) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  125508  ,  /* 3161 (397) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 3162 (397) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  130368  ,  /* 3163 (397) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18557870,  /* 3164 (397) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3165 (397) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 3166 (397) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3167 (397) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18426650,  /* 3168 (397) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18430698,  /* 3169 (397) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 3170 (397) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18435558,  /* 3171 (397) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  16484   ,  /* 3172 (398) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3173 (398) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3174 (398) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3175 (398) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 3176 (398) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3177 (398) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 3178 (398) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3179 (398) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 3180 (398) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3181 (398) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 3182 (398) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3183 (398) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 3184 (398) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3185 (398) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 3186 (398) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3187 (398) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3188 (399) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3189 (399) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3190 (399) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3191 (399) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  1492545 ,  /* 3192 (399) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 3] */
  1488497 ,  /* 3193 (399) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  1479585 ,  /* 3194 (399) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 3] */
  1483637 ,  /* 3195 (399) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18557870,  /* 3196 (399) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3197 (399) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 3198 (399) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3199 (399) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  17199579,  /* 3200 (399) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 9] */
  17195531,  /* 3201 (399) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  17186619,  /* 3202 (399) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 9] */
  17190671,  /* 3203 (399) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  16484   ,  /* 3204 (400) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3205 (400) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 3206 (400) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3207 (400) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3208 (400) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3209 (400) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 3210 (400) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3211 (400) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 3212 (400) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3213 (400) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 3214 (400) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3215 (400) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 3216 (400) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 3217 (400) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 3218 (400) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3219 (400) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3220 (401) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3221 (401) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 3222 (401) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3223 (401) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3224 (401) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3225 (401) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 3226 (401) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3227 (401) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 3228 (401) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 3229 (401) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  125508  ,  /* 3230 (401) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3231 (401) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 3232 (401) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 3233 (401) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  130368  ,  /* 3234 (401) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 3235 (401) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3236 (402) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 3237 (402) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 3238 (402) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 3239 (402) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 3240 (402) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 3241 (402) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 3242 (402) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 3243 (402) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18426486,  /* 3244 (402) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18435072,  /* 3245 (402) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18439770,  /* 3246 (402) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18431184,  /* 3247 (402) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 3248 (402) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18435560,  /* 3249 (402) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18439610,  /* 3250 (402) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18430700,  /* 3251 (402) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3252 (403) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3253 (403) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17014010,  /* 3254 (403) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 3255 (403) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  18426488,  /* 3256 (403) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3257 (403) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1428962 ,  /* 3258 (403) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 3259 (403) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1433496 ,  /* 3260 (403) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 3261 (403) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18430698,  /* 3262 (403) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 3263 (403) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  17009472,  /* 3264 (403) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 3265 (403) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  12270   ,  /* 3266 (403) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3267 (403) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3268 (404) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 3269 (404) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 3270 (404) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 3271 (404) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 3272 (404) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 3273 (404) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 3274 (404) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 3275 (404) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544584,  /* 3276 (404) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18553170,  /* 3277 (404) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557868,  /* 3278 (404) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18549282,  /* 3279 (404) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 3280 (404) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553658,  /* 3281 (404) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 3282 (404) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18548798,  /* 3283 (404) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3284 (405) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3285 (405) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17014010,  /* 3286 (405) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 3287 (405) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  18426488,  /* 3288 (405) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3289 (405) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1428962 ,  /* 3290 (405) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 3291 (405) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 3292 (405) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 3293 (405) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18548796,  /* 3294 (405) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3295 (405) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /* 3296 (405) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 3297 (405) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  130368  ,  /* 3298 (405) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 3299 (405) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3300 (406) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3301 (406) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 3302 (406) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3303 (406) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3304 (406) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3305 (406) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 3306 (406) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3307 (406) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18439608,  /* 3308 (406) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 3309 (406) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18430698,  /* 3310 (406) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 3311 (406) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 3312 (406) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 3313 (406) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18435558,  /* 3314 (406) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 3315 (406) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3316 (407) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3317 (407) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 3318 (407) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3319 (407) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3320 (407) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3321 (407) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 3322 (407) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3323 (407) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557706,  /* 3324 (407) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 3325 (407) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18548796,  /* 3326 (407) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3327 (407) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 3328 (407) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 3329 (407) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553656,  /* 3330 (407) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3331 (407) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3332 (408) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18431186,  /* 3333 (408) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  3200    ,  /* 3334 (408) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18435074,  /* 3335 (408) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 3336 (408) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  7410    ,  /* 3337 (408) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18426648,  /* 3338 (408) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  12270   ,  /* 3339 (408) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18426486,  /* 3340 (408) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  11784   ,  /* 3341 (408) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18439770,  /* 3342 (408) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  7896    ,  /* 3343 (408) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3344 (408) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18435560,  /* 3345 (408) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16322   ,  /* 3346 (408) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18430700,  /* 3347 (408) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3348 (409) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3349 (409) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17014010,  /* 3350 (409) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 3351 (409) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  18426488,  /* 3352 (409) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3353 (409) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1428962 ,  /* 3354 (409) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 3355 (409) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  17022432,  /* 3356 (409) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 3357 (409) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  7410    ,  /* 3358 (409) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3359 (409) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  1420536 ,  /* 3360 (409) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 3361 (409) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18435558,  /* 3362 (409) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 3363 (409) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3364 (410) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18431186,  /* 3365 (410) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  3200    ,  /* 3366 (410) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18435074,  /* 3367 (410) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18557706,  /* 3368 (410) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  125508  ,  /* 3369 (410) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18544746,  /* 3370 (410) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  130368  ,  /* 3371 (410) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18426486,  /* 3372 (410) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  11784   ,  /* 3373 (410) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18439770,  /* 3374 (410) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  7896    ,  /* 3375 (410) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  121460  ,  /* 3376 (410) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18553658,  /* 3377 (410) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  134420  ,  /* 3378 (410) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18548798,  /* 3379 (410) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3380 (411) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3381 (411) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17014010,  /* 3382 (411) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 3383 (411) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  18426488,  /* 3384 (411) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3385 (411) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1428962 ,  /* 3386 (411) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 3387 (411) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  17140530,  /* 3388 (411) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 3389 (411) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  125508  ,  /* 3390 (411) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3391 (411) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  1538634 ,  /* 3392 (411) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 3393 (411) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18553656,  /* 3394 (411) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3395 (411) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3396 (412) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3397 (412) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  125996  ,  /* 3398 (412) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  129882  ,  /* 3399 (412) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3200    ,  /* 3400 (412) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3401 (412) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  129884  ,  /* 3402 (412) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  125994  ,  /* 3403 (412) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 3404 (412) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3405 (412) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125508  ,  /* 3406 (412) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3407 (412) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 3408 (412) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 3409 (412) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  130368  ,  /* 3410 (412) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 3411 (412) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3412 (413) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3413 (413) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  125996  ,  /* 3414 (413) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  129882  ,  /* 3415 (413) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3200    ,  /* 3416 (413) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3417 (413) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  129884  ,  /* 3418 (413) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  125994  ,  /* 3419 (413) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  134418  ,  /* 3420 (413) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 3421 (413) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  7410    ,  /* 3422 (413) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3423 (413) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  121458  ,  /* 3424 (413) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 3425 (413) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  12270   ,  /* 3426 (413) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3427 (413) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3428 (414) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 3429 (414) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 3430 (414) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 3431 (414) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  134418  ,  /* 3432 (414) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  18430698,  /* 3433 (414) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  121458  ,  /* 3434 (414) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  18435558,  /* 3435 (414) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18544584,  /* 3436 (414) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  11784   ,  /* 3437 (414) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557868,  /* 3438 (414) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  7896    ,  /* 3439 (414) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18426650,  /* 3440 (414) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  130370  ,  /* 3441 (414) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  18439610,  /* 3442 (414) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  125510  ,  /* 3443 (414) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3444 (415) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3445 (415) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17132108,  /* 3446 (415) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135994,  /* 3447 (415) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /* 3448 (415) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3449 (415) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1547060 ,  /* 3450 (415) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1543170 ,  /* 3451 (415) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1433496 ,  /* 3452 (415) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 0] */
  1420538 ,  /* 3453 (415) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 0] */
  18548796,  /* 3454 (415) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3455 (415) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  17009472,  /* 3456 (415) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 3457 (415) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  130368  ,  /* 3458 (415) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 3459 (415) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3460 (416) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 3461 (416) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 3462 (416) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 3463 (416) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16320   ,  /* 3464 (416) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18548796,  /* 3465 (416) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3360    ,  /* 3466 (416) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18553656,  /* 3467 (416) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544584,  /* 3468 (416) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  11784   ,  /* 3469 (416) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557868,  /* 3470 (416) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  7896    ,  /* 3471 (416) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544748,  /* 3472 (416) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  12272   ,  /* 3473 (416) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557708,  /* 3474 (416) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  7412    ,  /* 3475 (416) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3476 (417) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3477 (417) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17132108,  /* 3478 (417) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135994,  /* 3479 (417) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /* 3480 (417) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3481 (417) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1547060 ,  /* 3482 (417) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1543170 ,  /* 3483 (417) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 3484 (417) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 3485 (417) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18430698,  /* 3486 (417) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 3487 (417) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  17127570,  /* 3488 (417) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 3489 (417) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  12270   ,  /* 3490 (417) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3491 (417) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3492 (418) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3493 (418) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  125996  ,  /* 3494 (418) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  129882  ,  /* 3495 (418) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3200    ,  /* 3496 (418) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3497 (418) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  129884  ,  /* 3498 (418) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  125994  ,  /* 3499 (418) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18439608,  /* 3500 (418) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 3501 (418) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18548796,  /* 3502 (418) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3503 (418) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 3504 (418) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 3505 (418) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18553656,  /* 3506 (418) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3507 (418) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3508 (419) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3509 (419) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18549284,  /* 3510 (419) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3511 (419) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3200    ,  /* 3512 (419) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3513 (419) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18553172,  /* 3514 (419) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3515 (419) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 3516 (419) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 3517 (419) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  7410    ,  /* 3518 (419) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3519 (419) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18544746,  /* 3520 (419) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 3521 (419) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  12270   ,  /* 3522 (419) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3523 (419) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3524 (420) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 3525 (420) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 3526 (420) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 3527 (420) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /* 3528 (420) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  7410    ,  /* 3529 (420) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544746,  /* 3530 (420) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  12270   ,  /* 3531 (420) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544584,  /* 3532 (420) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  11784   ,  /* 3533 (420) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557868,  /* 3534 (420) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  7896    ,  /* 3535 (420) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3536 (420) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18553658,  /* 3537 (420) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16322   ,  /* 3538 (420) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18548798,  /* 3539 (420) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3540 (421) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3541 (421) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17132108,  /* 3542 (421) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135994,  /* 3543 (421) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /* 3544 (421) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3545 (421) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1547060 ,  /* 3546 (421) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1543170 ,  /* 3547 (421) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  17022432,  /* 3548 (421) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 0] */
  17009474,  /* 3549 (421) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 0] */
  125508  ,  /* 3550 (421) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3551 (421) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  1420536 ,  /* 3552 (421) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 3553 (421) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  18553656,  /* 3554 (421) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3555 (421) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3556 (422) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 3557 (422) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 3558 (422) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 3559 (422) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18439608,  /* 3560 (422) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  125508  ,  /* 3561 (422) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18426648,  /* 3562 (422) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  130368  ,  /* 3563 (422) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  18544584,  /* 3564 (422) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  11784   ,  /* 3565 (422) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557868,  /* 3566 (422) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  7896    ,  /* 3567 (422) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  121460  ,  /* 3568 (422) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  18435560,  /* 3569 (422) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  134420  ,  /* 3570 (422) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  18430700,  /* 3571 (422) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3572 (423) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3573 (423) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17132108,  /* 3574 (423) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135994,  /* 3575 (423) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  18426488,  /* 3576 (423) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 3577 (423) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1547060 ,  /* 3578 (423) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1543170 ,  /* 3579 (423) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  17140530,  /* 3580 (423) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 3581 (423) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  7410    ,  /* 3582 (423) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3583 (423) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  1538634 ,  /* 3584 (423) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 3585 (423) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  18435558,  /* 3586 (423) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 3587 (423) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3588 (424) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3589 (424) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 3590 (424) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3591 (424) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3592 (424) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3593 (424) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 3594 (424) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3595 (424) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 3596 (424) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3597 (424) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 3598 (424) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3599 (424) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 3600 (424) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 3601 (424) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 3602 (424) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3603 (424) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557870,  /* 3604 (424) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 3605 (424) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18549284,  /* 3606 (424) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3607 (424) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 3608 (424) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 3609 (424) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18553172,  /* 3610 (424) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3611 (424) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 3612 (424) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 3613 (424) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18548796,  /* 3614 (424) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 3615 (424) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 3616 (424) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 3617 (424) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553656,  /* 3618 (424) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 3619 (424) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16484   ,  /* 3620 (425) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3621 (425) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 3622 (425) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 3623 (425) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 3624 (425) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 3625 (425) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 3626 (425) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 3627 (425) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  134418  ,  /* 3628 (425) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 3629 (425) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  125508  ,  /* 3630 (425) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3631 (425) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  121458  ,  /* 3632 (425) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  134420  ,  /* 3633 (425) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  130368  ,  /* 3634 (425) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 3635 (425) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  18557870,  /* 3636 (425) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 3637 (425) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18549284,  /* 3638 (425) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 3639 (425) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 3640 (425) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 3641 (425) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18553172,  /* 3642 (425) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 3643 (425) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 3644 (425) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 3645 (425) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18430698,  /* 3646 (425) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 3647 (425) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 3648 (425) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 3649 (425) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18435558,  /* 3650 (425) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 3651 (425) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  16484   ,  /* 3652 (426) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 3653 (426) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 3654 (426) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 3655 (426) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  1492545 ,  /* 3656 (426) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 3] */
  7410    ,  /* 3657 (426) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  17186619,  /* 3658 (426) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 9] */
  18553656,  /* 3659 (426) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  1479423 ,  /* 3660 (426) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 3] */
  11784   ,  /* 3661 (426) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  17199741,  /* 3662 (426) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 9] */
  18549282,  /* 3663 (426) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3362    ,  /* 3664 (426) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1488497 ,  /* 3665 (426) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18557708,  /* 3666 (426) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17190671,  /* 3667 (426) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  18557870,  /* 3668 (426) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 3669 (426) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 3670 (426) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 3671 (426) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  17199579,  /* 3672 (426) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 9] */
  18548796,  /* 3673 (426) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  1479585 ,  /* 3674 (426) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 3] */
  12270   ,  /* 3675 (426) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  17186457,  /* 3676 (426) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 9] */
  18553170,  /* 3677 (426) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  1492707 ,  /* 3678 (426) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 3] */
  7896    ,  /* 3679 (426) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544748,  /* 3680 (426) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17195531,  /* 3681 (426) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  16322   ,  /* 3682 (426) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1483637 ,  /* 3683 (426) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  16484   ,  /* 3684 (427) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3685 (427) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10695767,  /* 3686 (427) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 3] */
  10699653,  /* 3687 (427) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 3] */
  17127410,  /* 3688 (427) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 3689 (427) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  9400577 ,  /* 3690 (427) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 9] */
  9396687 ,  /* 3691 (427) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 9] */
  16320   ,  /* 3692 (427) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 3693 (427) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  9396201 ,  /* 3694 (427) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 9] */
  9401063 ,  /* 3695 (427) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 9] */
  17127570,  /* 3696 (427) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 3697 (427) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  10700139,  /* 3698 (427) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 3] */
  10695281,  /* 3699 (427) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 3] */
  18557870,  /* 3700 (427) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 3701 (427) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  26402801,  /* 3702 (427) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  26406687,  /* 3703 (427) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  1420376 ,  /* 3704 (427) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 3705 (427) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  27705767,  /* 3706 (427) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  27701877,  /* 3707 (427) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  18557706,  /* 3708 (427) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 3709 (427) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  27701391,  /* 3710 (427) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 3] */
  27706253,  /* 3711 (427) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 3] */
  1420536 ,  /* 3712 (427) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  1433498 ,  /* 3713 (427) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  26407173,  /* 3714 (427) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 9] */
  26402315,  /* 3715 (427) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 9] */
  16484   ,  /* 3716 (428) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1484123 ,  /* 3717 (428) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 3] */
  18544586,  /* 3718 (428) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  17195045,  /* 3719 (428) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 9] */
  17081481,  /* 3720 (428) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 3] */
  18430698,  /* 3721 (428) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  1597683 ,  /* 3722 (428) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 9] */
  130368  ,  /* 3723 (428) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  1479423 ,  /* 3724 (428) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 3] */
  11784   ,  /* 3725 (428) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  17199741,  /* 3726 (428) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 9] */
  18549282,  /* 3727 (428) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18426650,  /* 3728 (428) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17077433,  /* 3729 (428) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 3] */
  134420  ,  /* 3730 (428) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  1601735 ,  /* 3731 (428) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 9] */
  18557870,  /* 3732 (428) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17191157,  /* 3733 (428) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 9] */
  3200    ,  /* 3734 (428) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  1488011 ,  /* 3735 (428) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 3] */
  1610643 ,  /* 3736 (428) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 9] */
  125508  ,  /* 3737 (428) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  17068521,  /* 3738 (428) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 3] */
  18435558,  /* 3739 (428) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  17186457,  /* 3740 (428) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 9] */
  18553170,  /* 3741 (428) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  1492707 ,  /* 3742 (428) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 3] */
  7896    ,  /* 3743 (428) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  121460  ,  /* 3744 (428) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1606595 ,  /* 3745 (428) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 9] */
  18439610,  /* 3746 (428) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  17072573,  /* 3747 (428) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 3] */
  16484   ,  /* 3748 (429) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3749 (429) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10695767,  /* 3750 (429) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 3] */
  10699653,  /* 3751 (429) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 3] */
  17127410,  /* 3752 (429) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 3753 (429) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  9400577 ,  /* 3754 (429) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 9] */
  9396687 ,  /* 3755 (429) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 9] */
  134418  ,  /* 3756 (429) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 3757 (429) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  9278103 ,  /* 3758 (429) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  9282965 ,  /* 3759 (429) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 3] */
  17009472,  /* 3760 (429) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 3761 (429) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  10818237,  /* 3762 (429) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  10813379,  /* 3763 (429) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 9] */
  18557870,  /* 3764 (429) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 3765 (429) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  26402801,  /* 3766 (429) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  26406687,  /* 3767 (429) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  1420376 ,  /* 3768 (429) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 3769 (429) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  27705767,  /* 3770 (429) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  27701877,  /* 3771 (429) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  18439608,  /* 3772 (429) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 3773 (429) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  27819489,  /* 3774 (429) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 9] */
  27824351,  /* 3775 (429) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 9] */
  1538634 ,  /* 3776 (429) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 3777 (429) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  26289075,  /* 3778 (429) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 3] */
  26284217,  /* 3779 (429) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 3] */
  16484   ,  /* 3780 (430) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3781 (430) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3782 (430) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3783 (431) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  86549   ,  /* 3784 (431) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  162689  ,  /* 3785 (431) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  16484   ,  /* 3786 (432) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  165281  ,  /* 3787 (432) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  83957   ,  /* 3788 (432) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  16484   ,  /* 3789 (433) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3790 (433) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3791 (433) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  23714816,  /* 3792 (433) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23706149,  /* 3793 (433) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23703557,  /* 3794 (433) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  13400924,  /* 3795 (433) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13392257,  /* 3796 (433) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13389665,  /* 3797 (433) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  16484   ,  /* 3798 (434) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 3799 (434) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 3800 (434) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16484   ,  /* 3801 (435) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3802 (435) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 3803 (435) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 3804 (435) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 3805 (435) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 3806 (435) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 3807 (436) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3808 (436) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 3809 (436) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 3810 (436) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 3811 (436) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 3812 (436) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  23714816,  /* 3813 (436) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23701530,  /* 3814 (436) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  8, 4, 4] */
  23706149,  /* 3815 (436) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23710197,  /* 3816 (436) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  8, 4, 4] */
  23703557,  /* 3817 (436) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23712789,  /* 3818 (436) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  8, 4, 4] */
  13400924,  /* 3819 (436) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13387638,  /* 3820 (436) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  4, 8, 8] */
  13392257,  /* 3821 (436) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13396305,  /* 3822 (436) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  4, 8, 8] */
  13389665,  /* 3823 (436) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13398897,  /* 3824 (436) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  4, 8, 8] */
  16484   ,  /* 3825 (437) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3826 (437) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 3827 (437) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 3828 (437) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  12064   ,  /* 3829 (437) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 3830 (437) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16484   ,  /* 3831 (438) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3832 (438) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3833 (438) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 3834 (438) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  5547    ,  /* 3835 (438) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16563   ,  /* 3836 (438) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 3837 (439) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3838 (439) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3839 (439) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 3840 (439) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  14133   ,  /* 3841 (439) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3117    ,  /* 3842 (439) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 3843 (440) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  86549   ,  /* 3844 (440) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  162689  ,  /* 3845 (440) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  164874  ,  /* 3846 (440) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 8] */
  84279   ,  /* 3847 (440) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 4] */
  16563   ,  /* 3848 (440) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 3849 (441) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  86549   ,  /* 3850 (441) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  162689  ,  /* 3851 (441) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  12270   ,  /* 3852 (441) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  171597  ,  /* 3853 (441) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 8] */
  81849   ,  /* 3854 (441) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 4] */
  16484   ,  /* 3855 (442) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  165281  ,  /* 3856 (442) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  83957   ,  /* 3857 (442) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  86142   ,  /* 3858 (442) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 4] */
  163011  ,  /* 3859 (442) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 8] */
  16563   ,  /* 3860 (442) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 3861 (443) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  165281  ,  /* 3862 (443) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  83957   ,  /* 3863 (443) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  12270   ,  /* 3864 (443) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  92865   ,  /* 3865 (443) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 4] */
  160581  ,  /* 3866 (443) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 8] */
  16484   ,  /* 3867 (444) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3868 (444) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3869 (444) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 3870 (444) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  14133   ,  /* 3871 (444) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3117    ,  /* 3872 (444) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  23714816,  /* 3873 (444) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23706149,  /* 3874 (444) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23703557,  /* 3875 (444) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23710602,  /* 3876 (444) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  8, 4, 4] */
  23712465,  /* 3877 (444) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  8, 4, 4] */
  23701449,  /* 3878 (444) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  8, 4, 4] */
  13400924,  /* 3879 (444) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13392257,  /* 3880 (444) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13389665,  /* 3881 (444) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13396710,  /* 3882 (444) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  4, 8, 8] */
  13398573,  /* 3883 (444) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  4, 8, 8] */
  13387557,  /* 3884 (444) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  4, 8, 8] */
  16484   ,  /* 3885 (445) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 3886 (445) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 3887 (445) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7410    ,  /* 3888 (445) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3250    ,  /* 3889 (445) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  9022    ,  /* 3890 (445) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  16484   ,  /* 3891 (446) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3892 (446) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3893 (446) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7412    ,  /* 3894 (446) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  5549    ,  /* 3895 (446) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16565   ,  /* 3896 (446) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3897 (447) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3898 (447) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3899 (447) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  12272   ,  /* 3900 (447) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14135   ,  /* 3901 (447) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  3119    ,  /* 3902 (447) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3903 (448) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3904 (448) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3905 (448) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  125510  ,  /* 3906 (448) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  123647  ,  /* 3907 (448) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134663  ,  /* 3908 (448) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3909 (449) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3910 (449) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3911 (449) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  130370  ,  /* 3912 (449) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  132233  ,  /* 3913 (449) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  121217  ,  /* 3914 (449) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3915 (450) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3916 (450) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3917 (450) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7412    ,  /* 3918 (450) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  5549    ,  /* 3919 (450) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16565   ,  /* 3920 (450) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  23714816,  /* 3921 (450) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23706149,  /* 3922 (450) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23703557,  /* 3923 (450) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23705744,  /* 3924 (450) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23703881,  /* 3925 (450) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23714897,  /* 3926 (450) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  13400924,  /* 3927 (450) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13392257,  /* 3928 (450) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13389665,  /* 3929 (450) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13391852,  /* 3930 (450) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13389989,  /* 3931 (450) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13401005,  /* 3932 (450) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  16484   ,  /* 3933 (451) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 3934 (451) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 3935 (451) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  12272   ,  /* 3936 (451) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16432   ,  /* 3937 (451) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  10660   ,  /* 3938 (451) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  16484   ,  /* 3939 (452) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 3940 (452) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 3941 (452) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  125510  ,  /* 3942 (452) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  123647  ,  /* 3943 (452) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  134663  ,  /* 3944 (452) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  23714816,  /* 3945 (452) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23706149,  /* 3946 (452) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23703557,  /* 3947 (452) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23823842,  /* 3948 (452) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  8, 4,10] */
  23821979,  /* 3949 (452) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  8, 4,10] */
  23832995,  /* 3950 (452) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  8, 4,10] */
  13400924,  /* 3951 (452) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13392257,  /* 3952 (452) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13389665,  /* 3953 (452) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13273754,  /* 3954 (452) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 2] */
  13271891,  /* 3955 (452) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  4, 8, 2] */
  13282907,  /* 3956 (452) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  4, 8, 2] */
  16484   ,  /* 3957 (453) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 3958 (453) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 3959 (453) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  18553658,  /* 3960 (453) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557818,  /* 3961 (453) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18552046,  /* 3962 (453) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  16484   ,  /* 3963 (454) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3964 (454) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 3965 (454) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 3966 (454) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 3967 (454) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 3968 (454) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 3969 (454) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 3970 (454) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  5547    ,  /* 3971 (454) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  14135   ,  /* 3972 (454) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16563   ,  /* 3973 (454) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  3119    ,  /* 3974 (454) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3975 (455) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3976 (455) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 3977 (455) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 3978 (455) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 3979 (455) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 3980 (455) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  125508  ,  /* 3981 (455) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 3982 (455) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  123645  ,  /* 3983 (455) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  132233  ,  /* 3984 (455) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134661  ,  /* 3985 (455) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  121217  ,  /* 3986 (455) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 3987 (456) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 3988 (456) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 3989 (456) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 3990 (456) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 3991 (456) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 3992 (456) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 3993 (456) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 3994 (456) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14133   ,  /* 3995 (456) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  5549    ,  /* 3996 (456) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3117    ,  /* 3997 (456) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  16565   ,  /* 3998 (456) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 3999 (457) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4000 (457) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4001 (457) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4002 (457) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4003 (457) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4004 (457) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  130368  ,  /* 4005 (457) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 4006 (457) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  132231  ,  /* 4007 (457) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  123647  ,  /* 4008 (457) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121215  ,  /* 4009 (457) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 6] */
  134663  ,  /* 4010 (457) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4011 (458) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4012 (458) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4013 (458) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4014 (458) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4015 (458) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4016 (458) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 4017 (458) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 4018 (458) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14133   ,  /* 4019 (458) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  5549    ,  /* 4020 (458) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3117    ,  /* 4021 (458) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  16565   ,  /* 4022 (458) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  23714816,  /* 4023 (458) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23701530,  /* 4024 (458) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  8, 4, 4] */
  23706149,  /* 4025 (458) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23710197,  /* 4026 (458) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  8, 4, 4] */
  23703557,  /* 4027 (458) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23712789,  /* 4028 (458) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  8, 4, 4] */
  23710602,  /* 4029 (458) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  8, 4, 4] */
  23705744,  /* 4030 (458) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23712465,  /* 4031 (458) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  8, 4, 4] */
  23703881,  /* 4032 (458) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23701449,  /* 4033 (458) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  8, 4, 4] */
  23714897,  /* 4034 (458) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  13400924,  /* 4035 (458) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13387638,  /* 4036 (458) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  4, 8, 8] */
  13392257,  /* 4037 (458) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13396305,  /* 4038 (458) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  4, 8, 8] */
  13389665,  /* 4039 (458) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13398897,  /* 4040 (458) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  4, 8, 8] */
  13396710,  /* 4041 (458) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  4, 8, 8] */
  13391852,  /* 4042 (458) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13398573,  /* 4043 (458) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  4, 8, 8] */
  13389989,  /* 4044 (458) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13387557,  /* 4045 (458) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  4, 8, 8] */
  13401005,  /* 4046 (458) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  16484   ,  /* 4047 (459) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4048 (459) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4049 (459) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4050 (459) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  12064   ,  /* 4051 (459) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4052 (459) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7410    ,  /* 4053 (459) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 4054 (459) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3250    ,  /* 4055 (459) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  16432   ,  /* 4056 (459) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  9022    ,  /* 4057 (459) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10660   ,  /* 4058 (459) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  16484   ,  /* 4059 (460) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4060 (460) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4061 (460) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4062 (460) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4063 (460) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4064 (460) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  130368  ,  /* 4065 (460) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 4066 (460) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  132231  ,  /* 4067 (460) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  123647  ,  /* 4068 (460) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  121215  ,  /* 4069 (460) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 6] */
  134663  ,  /* 4070 (460) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  23714816,  /* 4071 (460) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  8, 4, 4] */
  23701530,  /* 4072 (460) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  8, 4, 4] */
  23706149,  /* 4073 (460) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  8, 4, 4] */
  23710197,  /* 4074 (460) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  8, 4, 4] */
  23703557,  /* 4075 (460) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  8, 4, 4] */
  23712789,  /* 4076 (460) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  8, 4, 4] */
  23828700,  /* 4077 (460) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  8, 4,10] */
  23823842,  /* 4078 (460) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  8, 4,10] */
  23830563,  /* 4079 (460) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  8, 4,10] */
  23821979,  /* 4080 (460) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  8, 4,10] */
  23819547,  /* 4081 (460) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  8, 4,10] */
  23832995,  /* 4082 (460) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  8, 4,10] */
  13400924,  /* 4083 (460) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  4, 8, 8] */
  13387638,  /* 4084 (460) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  4, 8, 8] */
  13392257,  /* 4085 (460) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  4, 8, 8] */
  13396305,  /* 4086 (460) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  4, 8, 8] */
  13389665,  /* 4087 (460) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 8] */
  13398897,  /* 4088 (460) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  4, 8, 8] */
  13278612,  /* 4089 (460) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  4, 8, 2] */
  13273754,  /* 4090 (460) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  4, 8, 2] */
  13280475,  /* 4091 (460) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  4, 8, 2] */
  13271891,  /* 4092 (460) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  4, 8, 2] */
  13269459,  /* 4093 (460) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  4, 8, 2] */
  13282907,  /* 4094 (460) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  4, 8, 2] */
  16484   ,  /* 4095 (461) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4096 (461) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4097 (461) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4098 (461) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  12064   ,  /* 4099 (461) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4100 (461) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18548796,  /* 4101 (461) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 4102 (461) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544636,  /* 4103 (461) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18557818,  /* 4104 (461) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18550408,  /* 4105 (461) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18552046,  /* 4106 (461) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  16484   ,  /* 4107 (462) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  14459   ,  /* 4108 (462) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 4109 (462) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4110 (462) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 4111 (462) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11867   ,  /* 4112 (462) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 4113 (463) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  53825   ,  /* 4114 (463) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 2] */
  86549   ,  /* 4115 (463) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  121298  ,  /* 4116 (463) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  162689  ,  /* 4117 (463) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  208697  ,  /* 4118 (463) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0,10] */
  16484   ,  /* 4119 (464) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  211289  ,  /* 4120 (464) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0,10] */
  165281  ,  /* 4121 (464) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  121298  ,  /* 4122 (464) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  83957   ,  /* 4123 (464) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  51233   ,  /* 4124 (464) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 2] */
  16484   ,  /* 4125 (465) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  93191   ,  /* 4126 (465) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 4] */
  165281  ,  /* 4127 (465) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  3200    ,  /* 4128 (465) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  83957   ,  /* 4129 (465) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  169331  ,  /* 4130 (465) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 8] */
  16484   ,  /* 4131 (466) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  171923  ,  /* 4132 (466) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 8] */
  86549   ,  /* 4133 (466) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  3200    ,  /* 4134 (466) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  162689  ,  /* 4135 (466) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  90599   ,  /* 4136 (466) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 4] */
  16484   ,  /* 4137 (467) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  132557  ,  /* 4138 (467) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  7817    ,  /* 4139 (467) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /* 4140 (467) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  5225    ,  /* 4141 (467) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  129965  ,  /* 4142 (467) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4143 (468) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  5223    ,  /* 4144 (468) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4145 (468) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4146 (468) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4147 (468) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7815    ,  /* 4148 (468) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4149 (469) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4150 (469) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  14459   ,  /* 4151 (469) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  5223    ,  /* 4152 (469) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4153 (469) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4154 (469) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 4155 (469) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4156 (469) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4157 (469) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4158 (469) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  11867   ,  /* 4159 (469) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  7815    ,  /* 4160 (469) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4161 (470) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4162 (470) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  132557  ,  /* 4163 (470) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  123321  ,  /* 4164 (470) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  7817    ,  /* 4165 (470) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4166 (470) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 4167 (470) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 4168 (470) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  5225    ,  /* 4169 (470) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4170 (470) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  129965  ,  /* 4171 (470) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  125913  ,  /* 4172 (470) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  16484   ,  /* 4173 (471) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  14459   ,  /* 4174 (471) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 4175 (471) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4176 (471) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 4177 (471) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11867   ,  /* 4178 (471) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 4179 (471) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3117    ,  /* 4180 (471) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  5547    ,  /* 4181 (471) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 4182 (471) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16563   ,  /* 4183 (471) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  14133   ,  /* 4184 (471) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4185 (472) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  53825   ,  /* 4186 (472) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 2] */
  86549   ,  /* 4187 (472) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  121298  ,  /* 4188 (472) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  162689  ,  /* 4189 (472) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  208697  ,  /* 4190 (472) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0,10] */
  204240  ,  /* 4191 (472) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0,10] */
  160581  ,  /* 4192 (472) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 8] */
  123645  ,  /* 4193 (472) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  91002   ,  /* 4194 (472) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 4] */
  55929   ,  /* 4195 (472) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 2] */
  14133   ,  /* 4196 (472) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4197 (473) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  211289  ,  /* 4198 (473) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0,10] */
  165281  ,  /* 4199 (473) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  121298  ,  /* 4200 (473) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  83957   ,  /* 4201 (473) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  51233   ,  /* 4202 (473) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 2] */
  46776   ,  /* 4203 (473) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 2] */
  81849   ,  /* 4204 (473) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 4] */
  123645  ,  /* 4205 (473) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  169734  ,  /* 4206 (473) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 8] */
  213393  ,  /* 4207 (473) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0,10] */
  14133   ,  /* 4208 (473) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4209 (474) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  93191   ,  /* 4210 (474) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 4] */
  165281  ,  /* 4211 (474) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 8] */
  3200    ,  /* 4212 (474) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  83957   ,  /* 4213 (474) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 4] */
  169331  ,  /* 4214 (474) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 8] */
  164874  ,  /* 4215 (474) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 8] */
  81849   ,  /* 4216 (474) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 4] */
  5547    ,  /* 4217 (474) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  169734  ,  /* 4218 (474) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 8] */
  95295   ,  /* 4219 (474) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 4] */
  14133   ,  /* 4220 (474) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4221 (475) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  171923  ,  /* 4222 (475) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 8] */
  86549   ,  /* 4223 (475) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 4] */
  3200    ,  /* 4224 (475) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  162689  ,  /* 4225 (475) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 8] */
  90599   ,  /* 4226 (475) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 4] */
  86142   ,  /* 4227 (475) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 4] */
  160581  ,  /* 4228 (475) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 8] */
  5547    ,  /* 4229 (475) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  91002   ,  /* 4230 (475) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 4] */
  174027  ,  /* 4231 (475) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 8] */
  14133   ,  /* 4232 (475) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4233 (476) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  132557  ,  /* 4234 (476) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  7817    ,  /* 4235 (476) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /* 4236 (476) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  5225    ,  /* 4237 (476) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  129965  ,  /* 4238 (476) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  125508  ,  /* 4239 (476) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3117    ,  /* 4240 (476) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  123645  ,  /* 4241 (476) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  12270   ,  /* 4242 (476) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  134661  ,  /* 4243 (476) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  14133   ,  /* 4244 (476) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4245 (477) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  14459   ,  /* 4246 (477) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 4247 (477) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4248 (477) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 4249 (477) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11867   ,  /* 4250 (477) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  12272   ,  /* 4251 (477) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16565   ,  /* 4252 (477) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  14135   ,  /* 4253 (477) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7412    ,  /* 4254 (477) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3119    ,  /* 4255 (477) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  5549    ,  /* 4256 (477) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 4257 (478) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  14459   ,  /* 4258 (478) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7817    ,  /* 4259 (478) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4260 (478) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  5225    ,  /* 4261 (478) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11867   ,  /* 4262 (478) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  130370  ,  /* 4263 (478) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  134663  ,  /* 4264 (478) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  132233  ,  /* 4265 (478) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  125510  ,  /* 4266 (478) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  121217  ,  /* 4267 (478) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  123647  ,  /* 4268 (478) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4269 (479) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  132557  ,  /* 4270 (479) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  7817    ,  /* 4271 (479) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /* 4272 (479) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  5225    ,  /* 4273 (479) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  129965  ,  /* 4274 (479) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  12272   ,  /* 4275 (479) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  134663  ,  /* 4276 (479) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  14135   ,  /* 4277 (479) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  125510  ,  /* 4278 (479) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3119    ,  /* 4279 (479) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  123647  ,  /* 4280 (479) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4281 (480) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  132557  ,  /* 4282 (480) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  7817    ,  /* 4283 (480) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  121298  ,  /* 4284 (480) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  5225    ,  /* 4285 (480) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  129965  ,  /* 4286 (480) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  130370  ,  /* 4287 (480) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16565   ,  /* 4288 (480) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  132233  ,  /* 4289 (480) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  7412    ,  /* 4290 (480) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  121217  ,  /* 4291 (480) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  5549    ,  /* 4292 (480) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 4293 (481) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  5223    ,  /* 4294 (481) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4295 (481) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4296 (481) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4297 (481) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7815    ,  /* 4298 (481) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 4299 (481) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16565   ,  /* 4300 (481) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5547    ,  /* 4301 (481) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 4302 (481) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16563   ,  /* 4303 (481) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  5549    ,  /* 4304 (481) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 4305 (482) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  123321  ,  /* 4306 (482) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  7817    ,  /* 4307 (482) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /* 4308 (482) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  5225    ,  /* 4309 (482) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  125913  ,  /* 4310 (482) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  7410    ,  /* 4311 (482) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  134663  ,  /* 4312 (482) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  5547    ,  /* 4313 (482) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  125510  ,  /* 4314 (482) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16563   ,  /* 4315 (482) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  123647  ,  /* 4316 (482) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4317 (483) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  5223    ,  /* 4318 (483) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4319 (483) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4320 (483) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4321 (483) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7815    ,  /* 4322 (483) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 4323 (483) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3117    ,  /* 4324 (483) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  14135   ,  /* 4325 (483) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 4326 (483) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3119    ,  /* 4327 (483) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  14133   ,  /* 4328 (483) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4329 (484) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  123321  ,  /* 4330 (484) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  7817    ,  /* 4331 (484) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  134580  ,  /* 4332 (484) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  5225    ,  /* 4333 (484) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  125913  ,  /* 4334 (484) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 4335 (484) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3117    ,  /* 4336 (484) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  132233  ,  /* 4337 (484) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  12270   ,  /* 4338 (484) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  121217  ,  /* 4339 (484) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  14133   ,  /* 4340 (484) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  16484   ,  /* 4341 (485) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4342 (485) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  14459   ,  /* 4343 (485) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  5223    ,  /* 4344 (485) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4345 (485) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4346 (485) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 4347 (485) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4348 (485) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4349 (485) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4350 (485) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  11867   ,  /* 4351 (485) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  7815    ,  /* 4352 (485) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 4353 (485) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 4354 (485) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3117    ,  /* 4355 (485) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  16565   ,  /* 4356 (485) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  5547    ,  /* 4357 (485) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  14135   ,  /* 4358 (485) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 4359 (485) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 4360 (485) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16563   ,  /* 4361 (485) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  3119    ,  /* 4362 (485) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  14133   ,  /* 4363 (485) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  5549    ,  /* 4364 (485) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 4365 (486) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4366 (486) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  14459   ,  /* 4367 (486) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  5223    ,  /* 4368 (486) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7817    ,  /* 4369 (486) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4370 (486) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 4371 (486) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4372 (486) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  5225    ,  /* 4373 (486) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4374 (486) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  11867   ,  /* 4375 (486) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  7815    ,  /* 4376 (486) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  125508  ,  /* 4377 (486) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 4378 (486) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  121215  ,  /* 4379 (486) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 6] */
  134663  ,  /* 4380 (486) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  123645  ,  /* 4381 (486) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  132233  ,  /* 4382 (486) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  130368  ,  /* 4383 (486) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 4384 (486) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  134661  ,  /* 4385 (486) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  121217  ,  /* 4386 (486) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  132231  ,  /* 4387 (486) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  123647  ,  /* 4388 (486) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4389 (487) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4390 (487) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  132557  ,  /* 4391 (487) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  123321  ,  /* 4392 (487) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  7817    ,  /* 4393 (487) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4394 (487) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 4395 (487) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 4396 (487) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  5225    ,  /* 4397 (487) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4398 (487) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  129965  ,  /* 4399 (487) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  125913  ,  /* 4400 (487) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  7410    ,  /* 4401 (487) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 4402 (487) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  121215  ,  /* 4403 (487) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 6] */
  134663  ,  /* 4404 (487) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 6] */
  5547    ,  /* 4405 (487) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  14135   ,  /* 4406 (487) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  130368  ,  /* 4407 (487) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 4408 (487) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  16563   ,  /* 4409 (487) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 0] */
  3119    ,  /* 4410 (487) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 0] */
  132231  ,  /* 4411 (487) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  123647  ,  /* 4412 (487) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  16484   ,  /* 4413 (488) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4414 (488) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  132557  ,  /* 4415 (488) [  1,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  123321  ,  /* 4416 (488) [ -1, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  7817    ,  /* 4417 (488) [  0,-1, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  11865   ,  /* 4418 (488) [  0, 1, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  121298  ,  /* 4419 (488) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  134580  ,  /* 4420 (488) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  5225    ,  /* 4421 (488) [ -1, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  14457   ,  /* 4422 (488) [  1,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  129965  ,  /* 4423 (488) [  0, 1, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  125913  ,  /* 4424 (488) [  0,-1, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  125508  ,  /* 4425 (488) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 4426 (488) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3117    ,  /* 4427 (488) [ -1, 0, 0,-1, 1, 0, 0, 0,-1,  0, 0, 0] */
  16565   ,  /* 4428 (488) [  1, 0, 0, 1,-1, 0, 0, 0, 1,  0, 0, 0] */
  123645  ,  /* 4429 (488) [ -1, 1, 0, 0, 1, 0, 0, 0,-1,  0, 0, 6] */
  132233  ,  /* 4430 (488) [  1,-1, 0, 0,-1, 0, 0, 0, 1,  0, 0, 6] */
  12270   ,  /* 4431 (488) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 4432 (488) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  134661  ,  /* 4433 (488) [  1, 0, 0, 1,-1, 0, 0, 0,-1,  0, 0, 6] */
  121217  ,  /* 4434 (488) [ -1, 0, 0,-1, 1, 0, 0, 0, 1,  0, 0, 6] */
  14133   ,  /* 4435 (488) [  1,-1, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  5549    ,  /* 4436 (488) [ -1, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  16484   ,  /* 4437 (489) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4438 (489) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 4439 (489) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 4440 (489) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4441 (489) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10324   ,  /* 4442 (489) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9352    ,  /* 4443 (489) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  8872    ,  /* 4444 (489) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 4445 (489) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4446 (489) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4447 (489) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4448 (489) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16484   ,  /* 4449 (490) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4450 (490) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 4451 (490) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 4452 (490) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4453 (490) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10324   ,  /* 4454 (490) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9352    ,  /* 4455 (490) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  8872    ,  /* 4456 (490) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 4457 (490) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4458 (490) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4459 (490) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4460 (490) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  1551758 ,  /* 4461 (490) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 4462 (490) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /* 4463 (490) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /* 4464 (490) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1546090 ,  /* 4465 (490) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1545598 ,  /* 4466 (490) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544626 ,  /* 4467 (490) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544146 ,  /* 4468 (490) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1547338 ,  /* 4469 (490) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542946 ,  /* 4470 (490) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1542910 ,  /* 4471 (490) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1547266 ,  /* 4472 (490) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  17140694,  /* 4473 (490) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 4474 (490) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 4475 (490) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /* 4476 (490) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17135026,  /* 4477 (490) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17134534,  /* 4478 (490) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133562,  /* 4479 (490) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133082,  /* 4480 (490) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17136274,  /* 4481 (490) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131882,  /* 4482 (490) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17131846,  /* 4483 (490) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136202,  /* 4484 (490) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  18439772,  /* 4485 (490) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 4486 (490) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 4487 (490) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 4488 (490) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18434104,  /* 4489 (490) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18433612,  /* 4490 (490) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432640,  /* 4491 (490) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432160,  /* 4492 (490) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18435352,  /* 4493 (490) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430960,  /* 4494 (490) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18430924,  /* 4495 (490) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435280,  /* 4496 (490) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  16484   ,  /* 4497 (491) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4498 (491) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 4499 (491) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 4500 (491) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4501 (491) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10324   ,  /* 4502 (491) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9352    ,  /* 4503 (491) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  8872    ,  /* 4504 (491) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 4505 (491) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4506 (491) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4507 (491) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4508 (491) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18557870,  /* 4509 (491) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 4510 (491) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /* 4511 (491) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544746,  /* 4512 (491) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18552202,  /* 4513 (491) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18551710,  /* 4514 (491) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550738,  /* 4515 (491) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550258,  /* 4516 (491) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18553450,  /* 4517 (491) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18549058,  /* 4518 (491) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18549022,  /* 4519 (491) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18553378,  /* 4520 (491) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  16484   ,  /* 4521 (492) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  17127410,  /* 4522 (492) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439608,  /* 4523 (492) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  1538634 ,  /* 4524 (492) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  10816   ,  /* 4525 (492) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18433612,  /* 4526 (492) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1544626 ,  /* 4527 (492) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  17133082,  /* 4528 (492) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  12064   ,  /* 4529 (492) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  1542946 ,  /* 4530 (492) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  17131846,  /* 4531 (492) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  18435280,  /* 4532 (492) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  16484   ,  /* 4533 (493) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  1420376 ,  /* 4534 (493) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  134418  ,  /* 4535 (493) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  1538634 ,  /* 4536 (493) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  10816   ,  /* 4537 (493) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  128422  ,  /* 4538 (493) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 6] */
  17015464,  /* 4539 (493) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 0] */
  17133082,  /* 4540 (493) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  12064   ,  /* 4541 (493) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  17013784,  /* 4542 (493) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 0] */
  1424812 ,  /* 4543 (493) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 0] */
  18435280,  /* 4544 (493) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18557870,  /* 4545 (493) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  17127410,  /* 4546 (493) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18439608,  /* 4547 (493) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  17009472,  /* 4548 (493) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  18552202,  /* 4549 (493) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18433612,  /* 4550 (493) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1544626 ,  /* 4551 (493) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1426048 ,  /* 4552 (493) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 0] */
  18553450,  /* 4553 (493) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  1542946 ,  /* 4554 (493) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  17131846,  /* 4555 (493) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  130090  ,  /* 4556 (493) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 6] */
  16484   ,  /* 4557 (494) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4558 (494) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 4559 (494) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4560 (494) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 4561 (494) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 4562 (494) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 4563 (494) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 4564 (494) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 4565 (494) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4566 (494) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10324   ,  /* 4567 (494) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 4568 (494) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9352    ,  /* 4569 (494) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 4570 (494) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8872    ,  /* 4571 (494) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 4572 (494) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  12064   ,  /* 4573 (494) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4574 (494) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4575 (494) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 4576 (494) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4577 (494) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 4578 (494) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4579 (494) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 4580 (494) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16484   ,  /* 4581 (495) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4582 (495) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 4583 (495) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 4584 (495) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4585 (495) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10324   ,  /* 4586 (495) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9352    ,  /* 4587 (495) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  8872    ,  /* 4588 (495) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 4589 (495) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4590 (495) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4591 (495) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4592 (495) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18544584,  /* 4593 (495) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18557868,  /* 4594 (495) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 4595 (495) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 4596 (495) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18550252,  /* 4597 (495) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550744,  /* 4598 (495) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18551716,  /* 4599 (495) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18552196,  /* 4600 (495) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18549004,  /* 4601 (495) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18553396,  /* 4602 (495) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18553432,  /* 4603 (495) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18549076,  /* 4604 (495) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  16484   ,  /* 4605 (496) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4606 (496) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18426488,  /* 4607 (496) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 4608 (496) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1551594 ,  /* 4609 (496) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 4610 (496) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  17127570,  /* 4611 (496) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 4612 (496) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  10816   ,  /* 4613 (496) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4614 (496) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  1545598 ,  /* 4615 (496) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544632 ,  /* 4616 (496) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  17133562,  /* 4617 (496) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134540,  /* 4618 (496) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  18432160,  /* 4619 (496) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18434098,  /* 4620 (496) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  12064   ,  /* 4621 (496) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4622 (496) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  17131882,  /* 4623 (496) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136220,  /* 4624 (496) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  18430924,  /* 4625 (496) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435334,  /* 4626 (496) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  1547266 ,  /* 4627 (496) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542964 ,  /* 4628 (496) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  16484   ,  /* 4629 (497) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4630 (497) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 4631 (497) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4632 (497) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 4633 (497) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 4634 (497) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 4635 (497) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 4636 (497) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 4637 (497) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4638 (497) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10324   ,  /* 4639 (497) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 4640 (497) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9352    ,  /* 4641 (497) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 4642 (497) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8872    ,  /* 4643 (497) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 4644 (497) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  12064   ,  /* 4645 (497) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4646 (497) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4647 (497) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 4648 (497) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4649 (497) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 4650 (497) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4651 (497) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 4652 (497) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  1551758 ,  /* 4653 (497) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 4654 (497) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 4655 (497) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 4656 (497) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 4657 (497) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 4658 (497) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 4659 (497) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 4660 (497) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1546090 ,  /* 4661 (497) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544140 ,  /* 4662 (497) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545598 ,  /* 4663 (497) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544632 ,  /* 4664 (497) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544626 ,  /* 4665 (497) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545604 ,  /* 4666 (497) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544146 ,  /* 4667 (497) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1546084 ,  /* 4668 (497) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1547338 ,  /* 4669 (497) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542892 ,  /* 4670 (497) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542946 ,  /* 4671 (497) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547284 ,  /* 4672 (497) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1542910 ,  /* 4673 (497) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1547320 ,  /* 4674 (497) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547266 ,  /* 4675 (497) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542964 ,  /* 4676 (497) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  17140694,  /* 4677 (497) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 4678 (497) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 4679 (497) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 4680 (497) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 4681 (497) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 4682 (497) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 4683 (497) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 4684 (497) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17135026,  /* 4685 (497) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133076,  /* 4686 (497) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134534,  /* 4687 (497) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133568,  /* 4688 (497) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133562,  /* 4689 (497) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134540,  /* 4690 (497) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133082,  /* 4691 (497) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17135020,  /* 4692 (497) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17136274,  /* 4693 (497) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131828,  /* 4694 (497) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17131882,  /* 4695 (497) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136220,  /* 4696 (497) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17131846,  /* 4697 (497) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136256,  /* 4698 (497) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136202,  /* 4699 (497) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17131900,  /* 4700 (497) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  18439772,  /* 4701 (497) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 4702 (497) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426488,  /* 4703 (497) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 4704 (497) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 4705 (497) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 4706 (497) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 4707 (497) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 4708 (497) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18434104,  /* 4709 (497) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432154,  /* 4710 (497) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433612,  /* 4711 (497) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432646,  /* 4712 (497) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432640,  /* 4713 (497) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433618,  /* 4714 (497) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432160,  /* 4715 (497) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18434098,  /* 4716 (497) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18435352,  /* 4717 (497) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430906,  /* 4718 (497) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430960,  /* 4719 (497) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18435298,  /* 4720 (497) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18430924,  /* 4721 (497) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435334,  /* 4722 (497) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18435280,  /* 4723 (497) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430978,  /* 4724 (497) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  16484   ,  /* 4725 (498) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 4726 (498) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 4727 (498) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 4728 (498) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 4729 (498) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10324   ,  /* 4730 (498) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9352    ,  /* 4731 (498) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  8872    ,  /* 4732 (498) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  12064   ,  /* 4733 (498) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4734 (498) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4735 (498) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4736 (498) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  9273891 ,  /* 4737 (498) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 3, 3] */
  9287175 ,  /* 4738 (498) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 3] */
  9274055 ,  /* 4739 (498) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 3, 3] */
  9287015 ,  /* 4740 (498) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 3] */
  9279559 ,  /* 4741 (498) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  3, 3, 3] */
  9280051 ,  /* 4742 (498) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  3, 3, 3] */
  9281023 ,  /* 4743 (498) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 3, 3] */
  9281503 ,  /* 4744 (498) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 3, 3] */
  9278311 ,  /* 4745 (498) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  3, 3, 3] */
  9282703 ,  /* 4746 (498) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 3, 3] */
  9282739 ,  /* 4747 (498) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 3, 3] */
  9278383 ,  /* 4748 (498) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  3, 3, 3] */
  1551758 ,  /* 4749 (498) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 4750 (498) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /* 4751 (498) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /* 4752 (498) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1546090 ,  /* 4753 (498) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1545598 ,  /* 4754 (498) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544626 ,  /* 4755 (498) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544146 ,  /* 4756 (498) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1547338 ,  /* 4757 (498) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542946 ,  /* 4758 (498) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1542910 ,  /* 4759 (498) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1547266 ,  /* 4760 (498) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  10809165,  /* 4761 (498) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 9, 9] */
  10822449,  /* 4762 (498) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 9] */
  10809329,  /* 4763 (498) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 9, 9] */
  10822289,  /* 4764 (498) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 9] */
  10814833,  /* 4765 (498) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  3, 9, 9] */
  10815325,  /* 4766 (498) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  3, 9, 9] */
  10816297,  /* 4767 (498) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 9, 9] */
  10816777,  /* 4768 (498) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 9, 9] */
  10813585,  /* 4769 (498) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  3, 9, 9] */
  10817977,  /* 4770 (498) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 9, 9] */
  10818013,  /* 4771 (498) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 9, 9] */
  10813657,  /* 4772 (498) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  3, 9, 9] */
  17140694,  /* 4773 (498) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 4774 (498) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 4775 (498) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /* 4776 (498) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17135026,  /* 4777 (498) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17134534,  /* 4778 (498) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133562,  /* 4779 (498) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133082,  /* 4780 (498) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17136274,  /* 4781 (498) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131882,  /* 4782 (498) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17131846,  /* 4783 (498) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136202,  /* 4784 (498) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  26398101,  /* 4785 (498) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 3, 9] */
  26411385,  /* 4786 (498) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 9] */
  26398265,  /* 4787 (498) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 3, 9] */
  26411225,  /* 4788 (498) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 9] */
  26403769,  /* 4789 (498) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  9, 3, 9] */
  26404261,  /* 4790 (498) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  9, 3, 9] */
  26405233,  /* 4791 (498) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 3, 9] */
  26405713,  /* 4792 (498) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 3, 9] */
  26402521,  /* 4793 (498) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  9, 3, 9] */
  26406913,  /* 4794 (498) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 3, 9] */
  26406949,  /* 4795 (498) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 3, 9] */
  26402593,  /* 4796 (498) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  9, 3, 9] */
  18439772,  /* 4797 (498) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 4798 (498) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 4799 (498) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 4800 (498) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18434104,  /* 4801 (498) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18433612,  /* 4802 (498) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432640,  /* 4803 (498) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432160,  /* 4804 (498) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18435352,  /* 4805 (498) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430960,  /* 4806 (498) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18430924,  /* 4807 (498) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435280,  /* 4808 (498) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  27697179,  /* 4809 (498) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 9, 3] */
  27710463,  /* 4810 (498) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 3] */
  27697343,  /* 4811 (498) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 9, 3] */
  27710303,  /* 4812 (498) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 3] */
  27702847,  /* 4813 (498) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  9, 9, 3] */
  27703339,  /* 4814 (498) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  9, 9, 3] */
  27704311,  /* 4815 (498) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 9, 3] */
  27704791,  /* 4816 (498) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 9, 3] */
  27701599,  /* 4817 (498) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  9, 9, 3] */
  27705991,  /* 4818 (498) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 9, 3] */
  27706027,  /* 4819 (498) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 9, 3] */
  27701671,  /* 4820 (498) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  9, 9, 3] */
  16484   ,  /* 4821 (499) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4822 (499) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  9214844 ,  /* 4823 (499) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 0] */
  9228126 ,  /* 4824 (499) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 0] */
  783957  ,  /* 4825 (499) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 3, 3] */
  770999  ,  /* 4826 (499) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 3, 3] */
  8565465 ,  /* 4827 (499) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 0, 3] */
  8578427 ,  /* 4828 (499) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 0, 3] */
  10816   ,  /* 4829 (499) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4830 (499) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  777961  ,  /* 4831 (499) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 3, 3] */
  776995  ,  /* 4832 (499) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 3, 3] */
  8571457 ,  /* 4833 (499) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  3, 0, 3] */
  8572435 ,  /* 4834 (499) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 0, 3] */
  9220516 ,  /* 4835 (499) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  3, 3, 0] */
  9222454 ,  /* 4836 (499) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 3, 0] */
  12064   ,  /* 4837 (499) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4838 (499) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  8569777 ,  /* 4839 (499) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  3, 0, 3] */
  8574115 ,  /* 4840 (499) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 0, 3] */
  9219280 ,  /* 4841 (499) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  3, 3, 0] */
  9223690 ,  /* 4842 (499) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 3, 0] */
  779629  ,  /* 4843 (499) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 3, 3] */
  775327  ,  /* 4844 (499) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 3, 3] */
  1551758 ,  /* 4845 (499) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 4846 (499) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  10750118,  /* 4847 (499) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 6] */
  10763400,  /* 4848 (499) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 6] */
  2319231 ,  /* 4849 (499) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 9, 9] */
  2306273 ,  /* 4850 (499) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 9, 9] */
  10100739,  /* 4851 (499) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 6, 9] */
  10113701,  /* 4852 (499) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 6, 9] */
  1546090 ,  /* 4853 (499) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544140 ,  /* 4854 (499) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  2313235 ,  /* 4855 (499) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 9, 9] */
  2312269 ,  /* 4856 (499) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 9, 9] */
  10106731,  /* 4857 (499) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  3, 6, 9] */
  10107709,  /* 4858 (499) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 6, 9] */
  10755790,  /* 4859 (499) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  3, 9, 6] */
  10757728,  /* 4860 (499) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 9, 6] */
  1547338 ,  /* 4861 (499) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542892 ,  /* 4862 (499) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  10105051,  /* 4863 (499) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  3, 6, 9] */
  10109389,  /* 4864 (499) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 6, 9] */
  10754554,  /* 4865 (499) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  3, 9, 6] */
  10758964,  /* 4866 (499) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 9, 6] */
  2314903 ,  /* 4867 (499) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 9, 9] */
  2310601 ,  /* 4868 (499) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 9, 9] */
  17140694,  /* 4869 (499) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 4870 (499) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  26339054,  /* 4871 (499) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 6] */
  26352336,  /* 4872 (499) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 6] */
  17908167,  /* 4873 (499) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 3, 9] */
  17895209,  /* 4874 (499) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 3, 9] */
  25689675,  /* 4875 (499) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 0, 9] */
  25702637,  /* 4876 (499) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 0, 9] */
  17135026,  /* 4877 (499) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133076,  /* 4878 (499) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17902171,  /* 4879 (499) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 3, 9] */
  17901205,  /* 4880 (499) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 3, 9] */
  25695667,  /* 4881 (499) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  9, 0, 9] */
  25696645,  /* 4882 (499) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 0, 9] */
  26344726,  /* 4883 (499) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  9, 3, 6] */
  26346664,  /* 4884 (499) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 3, 6] */
  17136274,  /* 4885 (499) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131828,  /* 4886 (499) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  25693987,  /* 4887 (499) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  9, 0, 9] */
  25698325,  /* 4888 (499) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 0, 9] */
  26343490,  /* 4889 (499) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  9, 3, 6] */
  26347900,  /* 4890 (499) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 3, 6] */
  17903839,  /* 4891 (499) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 3, 9] */
  17899537,  /* 4892 (499) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 3, 9] */
  18439772,  /* 4893 (499) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 4894 (499) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  27638132,  /* 4895 (499) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 0] */
  27651414,  /* 4896 (499) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 0] */
  19207245,  /* 4897 (499) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 9, 3] */
  19194287,  /* 4898 (499) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 9, 3] */
  26988753,  /* 4899 (499) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 6, 3] */
  27001715,  /* 4900 (499) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 6, 3] */
  18434104,  /* 4901 (499) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432154,  /* 4902 (499) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  19201249,  /* 4903 (499) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 9, 3] */
  19200283,  /* 4904 (499) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 9, 3] */
  26994745,  /* 4905 (499) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  9, 6, 3] */
  26995723,  /* 4906 (499) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 6, 3] */
  27643804,  /* 4907 (499) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  9, 9, 0] */
  27645742,  /* 4908 (499) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 9, 0] */
  18435352,  /* 4909 (499) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430906,  /* 4910 (499) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  26993065,  /* 4911 (499) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  9, 6, 3] */
  26997403,  /* 4912 (499) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 6, 3] */
  27642568,  /* 4913 (499) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  9, 9, 0] */
  27646978,  /* 4914 (499) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 9, 0] */
  19202917,  /* 4915 (499) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 9, 3] */
  19198615,  /* 4916 (499) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 9, 3] */
  16484   ,  /* 4917 (500) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4918 (500) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 4919 (500) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 4920 (500) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 4921 (500) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 4922 (500) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 4923 (500) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 4924 (500) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 4925 (500) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4926 (500) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10324   ,  /* 4927 (500) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 4928 (500) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9352    ,  /* 4929 (500) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 4930 (500) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8872    ,  /* 4931 (500) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 4932 (500) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  12064   ,  /* 4933 (500) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4934 (500) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7672    ,  /* 4935 (500) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 4936 (500) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  7636    ,  /* 4937 (500) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 4938 (500) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  11992   ,  /* 4939 (500) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 4940 (500) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  18557870,  /* 4941 (500) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 4942 (500) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 4943 (500) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 4944 (500) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 4945 (500) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 4946 (500) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 4947 (500) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 4948 (500) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18552202,  /* 4949 (500) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550252,  /* 4950 (500) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18551710,  /* 4951 (500) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550744,  /* 4952 (500) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550738,  /* 4953 (500) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18551716,  /* 4954 (500) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550258,  /* 4955 (500) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18552196,  /* 4956 (500) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18553450,  /* 4957 (500) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18549004,  /* 4958 (500) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18549058,  /* 4959 (500) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18553396,  /* 4960 (500) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18549022,  /* 4961 (500) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18553432,  /* 4962 (500) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18553378,  /* 4963 (500) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18549076,  /* 4964 (500) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  16484   ,  /* 4965 (501) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4966 (501) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17127410,  /* 4967 (501) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 4968 (501) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 4969 (501) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 4970 (501) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  1538634 ,  /* 4971 (501) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 4972 (501) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  10816   ,  /* 4973 (501) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4974 (501) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18433612,  /* 4975 (501) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432646,  /* 4976 (501) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  1544626 ,  /* 4977 (501) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545604 ,  /* 4978 (501) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  17133082,  /* 4979 (501) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17135020,  /* 4980 (501) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  12064   ,  /* 4981 (501) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 4982 (501) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  1542946 ,  /* 4983 (501) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547284 ,  /* 4984 (501) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  17131846,  /* 4985 (501) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136256,  /* 4986 (501) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  18435280,  /* 4987 (501) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430978,  /* 4988 (501) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  16484   ,  /* 4989 (502) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 4990 (502) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1420376 ,  /* 4991 (502) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 4992 (502) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  134418  ,  /* 4993 (502) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 4994 (502) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  1538634 ,  /* 4995 (502) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 4996 (502) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  10816   ,  /* 4997 (502) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 4998 (502) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  128422  ,  /* 4999 (502) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 6] */
  127456  ,  /* 5000 (502) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 6] */
  17015464,  /* 5001 (502) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 0] */
  17016442,  /* 5002 (502) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 0] */
  17133082,  /* 5003 (502) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17135020,  /* 5004 (502) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  12064   ,  /* 5005 (502) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 5006 (502) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  17013784,  /* 5007 (502) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 0] */
  17018122,  /* 5008 (502) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 0] */
  1424812 ,  /* 5009 (502) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 0] */
  1429222 ,  /* 5010 (502) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 0] */
  18435280,  /* 5011 (502) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430978,  /* 5012 (502) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18557870,  /* 5013 (502) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 5014 (502) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  17127410,  /* 5015 (502) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 5016 (502) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  18439608,  /* 5017 (502) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 5018 (502) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  17009472,  /* 5019 (502) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 5020 (502) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  18552202,  /* 5021 (502) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550252,  /* 5022 (502) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18433612,  /* 5023 (502) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432646,  /* 5024 (502) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  1544626 ,  /* 5025 (502) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545604 ,  /* 5026 (502) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1426048 ,  /* 5027 (502) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 0] */
  1427986 ,  /* 5028 (502) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 6, 0] */
  18553450,  /* 5029 (502) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18549004,  /* 5030 (502) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  1542946 ,  /* 5031 (502) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547284 ,  /* 5032 (502) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  17131846,  /* 5033 (502) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136256,  /* 5034 (502) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  130090  ,  /* 5035 (502) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 6] */
  125788  ,  /* 5036 (502) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 6] */
  16484   ,  /* 5037 (503) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 5038 (503) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 5039 (503) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 5040 (503) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 5041 (503) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 5042 (503) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 5043 (503) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 5044 (503) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 5045 (503) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10498   ,  /* 5046 (503) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5047 (503) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10642   ,  /* 5048 (503) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5049 (503) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9022    ,  /* 5050 (503) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5051 (503) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9202    ,  /* 5052 (503) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5053 (503) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16426   ,  /* 5054 (503) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  7672    ,  /* 5055 (503) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3310    ,  /* 5056 (503) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  7636    ,  /* 5057 (503) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3250    ,  /* 5058 (503) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  11992   ,  /* 5059 (503) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16378   ,  /* 5060 (503) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  16484   ,  /* 5061 (504) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 5062 (504) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 5063 (504) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 5064 (504) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16320   ,  /* 5065 (504) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18548796,  /* 5066 (504) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3360    ,  /* 5067 (504) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18553656,  /* 5068 (504) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  10816   ,  /* 5069 (504) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18551884,  /* 5070 (504) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  10324   ,  /* 5071 (504) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18552028,  /* 5072 (504) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  9352    ,  /* 5073 (504) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18550408,  /* 5074 (504) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  8872    ,  /* 5075 (504) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18550588,  /* 5076 (504) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  12064   ,  /* 5077 (504) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  18557812,  /* 5078 (504) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  7672    ,  /* 5079 (504) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  18544696,  /* 5080 (504) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  7636    ,  /* 5081 (504) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  18544636,  /* 5082 (504) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  11992   ,  /* 5083 (504) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18557764,  /* 5084 (504) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  16484   ,  /* 5085 (505) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 5086 (505) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 5087 (505) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 5088 (505) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 5089 (505) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 5090 (505) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 5091 (505) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 5092 (505) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 5093 (505) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10498   ,  /* 5094 (505) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5095 (505) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10642   ,  /* 5096 (505) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5097 (505) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9022    ,  /* 5098 (505) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5099 (505) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9202    ,  /* 5100 (505) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5101 (505) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16426   ,  /* 5102 (505) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  7672    ,  /* 5103 (505) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3310    ,  /* 5104 (505) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  7636    ,  /* 5105 (505) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3250    ,  /* 5106 (505) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  11992   ,  /* 5107 (505) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16378   ,  /* 5108 (505) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  1551758 ,  /* 5109 (505) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1543172 ,  /* 5110 (505) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1538474 ,  /* 5111 (505) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1547060 ,  /* 5112 (505) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1551594 ,  /* 5113 (505) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1542684 ,  /* 5114 (505) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1538634 ,  /* 5115 (505) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1547544 ,  /* 5116 (505) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1546090 ,  /* 5117 (505) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1545772 ,  /* 5118 (505) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 6, 6] */
  1545598 ,  /* 5119 (505) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545916 ,  /* 5120 (505) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 6, 6] */
  1544626 ,  /* 5121 (505) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544296 ,  /* 5122 (505) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  1544146 ,  /* 5123 (505) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544476 ,  /* 5124 (505) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  1547338 ,  /* 5125 (505) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1551700 ,  /* 5126 (505) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 6] */
  1542946 ,  /* 5127 (505) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1538584 ,  /* 5128 (505) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  1542910 ,  /* 5129 (505) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1538524 ,  /* 5130 (505) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  1547266 ,  /* 5131 (505) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1551652 ,  /* 5132 (505) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 6] */
  17140694,  /* 5133 (505) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17132108,  /* 5134 (505) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17127410,  /* 5135 (505) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17135996,  /* 5136 (505) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17140530,  /* 5137 (505) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17131620,  /* 5138 (505) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17127570,  /* 5139 (505) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17136480,  /* 5140 (505) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17135026,  /* 5141 (505) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17134708,  /* 5142 (505) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 0, 6] */
  17134534,  /* 5143 (505) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134852,  /* 5144 (505) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 0, 6] */
  17133562,  /* 5145 (505) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133232,  /* 5146 (505) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 0, 6] */
  17133082,  /* 5147 (505) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133412,  /* 5148 (505) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 0, 6] */
  17136274,  /* 5149 (505) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17140636,  /* 5150 (505) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 6] */
  17131882,  /* 5151 (505) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17127520,  /* 5152 (505) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 6] */
  17131846,  /* 5153 (505) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17127460,  /* 5154 (505) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 6] */
  17136202,  /* 5155 (505) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17140588,  /* 5156 (505) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 6] */
  18439772,  /* 5157 (505) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18431186,  /* 5158 (505) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426488,  /* 5159 (505) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18435074,  /* 5160 (505) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18439608,  /* 5161 (505) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18430698,  /* 5162 (505) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18426648,  /* 5163 (505) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18435558,  /* 5164 (505) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18434104,  /* 5165 (505) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18433786,  /* 5166 (505) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 0] */
  18433612,  /* 5167 (505) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433930,  /* 5168 (505) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 0] */
  18432640,  /* 5169 (505) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432310,  /* 5170 (505) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 0] */
  18432160,  /* 5171 (505) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432490,  /* 5172 (505) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 0] */
  18435352,  /* 5173 (505) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18439714,  /* 5174 (505) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 0] */
  18430960,  /* 5175 (505) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18426598,  /* 5176 (505) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 0] */
  18430924,  /* 5177 (505) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18426538,  /* 5178 (505) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 0] */
  18435280,  /* 5179 (505) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18439666,  /* 5180 (505) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 0] */
  16484   ,  /* 5181 (506) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  9278591 ,  /* 5182 (506) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 3] */
  1538474 ,  /* 5183 (506) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  26406689,  /* 5184 (506) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 9] */
  16320   ,  /* 5185 (506) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10813377,  /* 5186 (506) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 9] */
  1538634 ,  /* 5187 (506) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  27706251,  /* 5188 (506) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  10816   ,  /* 5189 (506) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9281191 ,  /* 5190 (506) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 3, 3] */
  17134534,  /* 5191 (506) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  27704623,  /* 5192 (506) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 9, 3] */
  9352    ,  /* 5193 (506) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  26403925,  /* 5194 (506) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 3, 9] */
  17133082,  /* 5195 (506) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  10815169,  /* 5196 (506) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 9, 9] */
  12064   ,  /* 5197 (506) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  9287119 ,  /* 5198 (506) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 3] */
  18430960,  /* 5199 (506) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  10809277,  /* 5200 (506) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 9, 9] */
  7636    ,  /* 5201 (506) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  27697231,  /* 5202 (506) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 9, 3] */
  18435280,  /* 5203 (506) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  26411281,  /* 5204 (506) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 9] */
  1551758 ,  /* 5205 (506) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  10813865,  /* 5206 (506) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 9] */
  3200    ,  /* 5207 (506) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  27705767,  /* 5208 (506) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  1551594 ,  /* 5209 (506) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  9278103 ,  /* 5210 (506) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  3360    ,  /* 5211 (506) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  26407173,  /* 5212 (506) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 9] */
  1546090 ,  /* 5213 (506) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  10816465,  /* 5214 (506) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 9, 9] */
  18433612,  /* 5215 (506) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26405545,  /* 5216 (506) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 9] */
  1544626 ,  /* 5217 (506) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27703003,  /* 5218 (506) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 9, 3] */
  18432160,  /* 5219 (506) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  9279895 ,  /* 5220 (506) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 3, 3] */
  1547338 ,  /* 5221 (506) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  10822393,  /* 5222 (506) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 9] */
  17131882,  /* 5223 (506) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  9274003 ,  /* 5224 (506) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 3, 3] */
  1542910 ,  /* 5225 (506) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  26398153,  /* 5226 (506) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 3, 9] */
  17136202,  /* 5227 (506) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  27710359,  /* 5228 (506) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 3] */
  17140694,  /* 5229 (506) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  26402801,  /* 5230 (506) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  18426488,  /* 5231 (506) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  9282479 ,  /* 5232 (506) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 3] */
  17140530,  /* 5233 (506) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  27701391,  /* 5234 (506) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 3] */
  18426648,  /* 5235 (506) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  10818237,  /* 5236 (506) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  17135026,  /* 5237 (506) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  26405401,  /* 5238 (506) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 3, 9] */
  10324   ,  /* 5239 (506) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10816609,  /* 5240 (506) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 9] */
  17133562,  /* 5241 (506) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  9279715 ,  /* 5242 (506) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 3, 3] */
  8872    ,  /* 5243 (506) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  27703183,  /* 5244 (506) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 9, 3] */
  17136274,  /* 5245 (506) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  26411329,  /* 5246 (506) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 9] */
  1542946 ,  /* 5247 (506) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  27697291,  /* 5248 (506) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 9, 3] */
  17131846,  /* 5249 (506) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  10809217,  /* 5250 (506) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 9, 9] */
  1547266 ,  /* 5251 (506) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  9287071 ,  /* 5252 (506) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 3] */
  18439772,  /* 5253 (506) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  27701879,  /* 5254 (506) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 3] */
  17127410,  /* 5255 (506) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  10817753,  /* 5256 (506) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 9] */
  18439608,  /* 5257 (506) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  26402313,  /* 5258 (506) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  17127570,  /* 5259 (506) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  9282963 ,  /* 5260 (506) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 3] */
  18434104,  /* 5261 (506) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  27704479,  /* 5262 (506) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 9, 3] */
  1545598 ,  /* 5263 (506) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  9281335 ,  /* 5264 (506) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 3, 3] */
  18432640,  /* 5265 (506) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  10814989,  /* 5266 (506) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 9, 9] */
  1544146 ,  /* 5267 (506) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  26404105,  /* 5268 (506) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 3, 9] */
  18435352,  /* 5269 (506) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  27710407,  /* 5270 (506) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 3] */
  7672    ,  /* 5271 (506) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  26398213,  /* 5272 (506) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 3, 9] */
  18430924,  /* 5273 (506) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  9273943 ,  /* 5274 (506) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 3, 3] */
  11992   ,  /* 5275 (506) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  10822345,  /* 5276 (506) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 9] */
  16484   ,  /* 5277 (507) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 5278 (507) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 5279 (507) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 5280 (507) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 5281 (507) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 5282 (507) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 5283 (507) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 5284 (507) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 5285 (507) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10498   ,  /* 5286 (507) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5287 (507) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10642   ,  /* 5288 (507) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5289 (507) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9022    ,  /* 5290 (507) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5291 (507) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9202    ,  /* 5292 (507) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5293 (507) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16426   ,  /* 5294 (507) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  7672    ,  /* 5295 (507) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3310    ,  /* 5296 (507) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  7636    ,  /* 5297 (507) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3250    ,  /* 5298 (507) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  11992   ,  /* 5299 (507) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16378   ,  /* 5300 (507) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  18557870,  /* 5301 (507) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18549284,  /* 5302 (507) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544586,  /* 5303 (507) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553172,  /* 5304 (507) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557706,  /* 5305 (507) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18548796,  /* 5306 (507) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544746,  /* 5307 (507) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18553656,  /* 5308 (507) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18552202,  /* 5309 (507) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18551884,  /* 5310 (507) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18551710,  /* 5311 (507) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18552028,  /* 5312 (507) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18550738,  /* 5313 (507) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550408,  /* 5314 (507) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18550258,  /* 5315 (507) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550588,  /* 5316 (507) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18553450,  /* 5317 (507) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18557812,  /* 5318 (507) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  18549058,  /* 5319 (507) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18544696,  /* 5320 (507) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18549022,  /* 5321 (507) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18544636,  /* 5322 (507) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18553378,  /* 5323 (507) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18557764,  /* 5324 (507) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  16484   ,  /* 5325 (508) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  26402801,  /* 5326 (508) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  17127410,  /* 5327 (508) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  27705767,  /* 5328 (508) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  18439608,  /* 5329 (508) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  9278103 ,  /* 5330 (508) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  1538634 ,  /* 5331 (508) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  10818237,  /* 5332 (508) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  10816   ,  /* 5333 (508) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  27704479,  /* 5334 (508) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 9, 3] */
  18433612,  /* 5335 (508) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  10816609,  /* 5336 (508) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 9] */
  1544626 ,  /* 5337 (508) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  9279715 ,  /* 5338 (508) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 3, 3] */
  17133082,  /* 5339 (508) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  26404105,  /* 5340 (508) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 3, 9] */
  12064   ,  /* 5341 (508) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  10822393,  /* 5342 (508) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 9] */
  1542946 ,  /* 5343 (508) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  26398213,  /* 5344 (508) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 3, 9] */
  17131846,  /* 5345 (508) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  9273943 ,  /* 5346 (508) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 3, 3] */
  18435280,  /* 5347 (508) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  27710359,  /* 5348 (508) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 3] */
  16484   ,  /* 5349 (509) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10695767,  /* 5350 (509) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 3] */
  17127410,  /* 5351 (509) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  9400577 ,  /* 5352 (509) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 9] */
  18439608,  /* 5353 (509) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  27819489,  /* 5354 (509) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 9] */
  1538634 ,  /* 5355 (509) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  26289075,  /* 5356 (509) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 3] */
  10816   ,  /* 5357 (509) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9399289 ,  /* 5358 (509) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 3, 9] */
  18433612,  /* 5359 (509) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26287447,  /* 5360 (509) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 3] */
  1544626 ,  /* 5361 (509) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27821101,  /* 5362 (509) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 9, 9] */
  17133082,  /* 5363 (509) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  10697071,  /* 5364 (509) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 9, 3] */
  12064   ,  /* 5365 (509) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  26293231,  /* 5366 (509) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 3] */
  1542946 ,  /* 5367 (509) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  10691179,  /* 5368 (509) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 9, 3] */
  17131846,  /* 5369 (509) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  27815329,  /* 5370 (509) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 9, 9] */
  18435280,  /* 5371 (509) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  9405169 ,  /* 5372 (509) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 9] */
  16484   ,  /* 5373 (510) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10695767,  /* 5374 (510) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 3] */
  17127410,  /* 5375 (510) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  9400577 ,  /* 5376 (510) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 9] */
  134418  ,  /* 5377 (510) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  9278103 ,  /* 5378 (510) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  17009472,  /* 5379 (510) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  10818237,  /* 5380 (510) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  10816   ,  /* 5381 (510) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9399289 ,  /* 5382 (510) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 3, 9] */
  18433612,  /* 5383 (510) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26287447,  /* 5384 (510) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 3] */
  17015464,  /* 5385 (510) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 0] */
  9279715 ,  /* 5386 (510) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 3, 3] */
  1426048 ,  /* 5387 (510) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 0] */
  26404105,  /* 5388 (510) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 3, 9] */
  12064   ,  /* 5389 (510) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  26293231,  /* 5390 (510) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 3] */
  1542946 ,  /* 5391 (510) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  10691179,  /* 5392 (510) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 9, 3] */
  1424812 ,  /* 5393 (510) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 0] */
  9273943 ,  /* 5394 (510) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 3, 3] */
  130090  ,  /* 5395 (510) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 6] */
  27710359,  /* 5396 (510) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 3] */
  18557870,  /* 5397 (510) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  26402801,  /* 5398 (510) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  1420376 ,  /* 5399 (510) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  27705767,  /* 5400 (510) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  18439608,  /* 5401 (510) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  27819489,  /* 5402 (510) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 9] */
  1538634 ,  /* 5403 (510) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  26289075,  /* 5404 (510) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 3] */
  18552202,  /* 5405 (510) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  27704479,  /* 5406 (510) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 9, 3] */
  128422  ,  /* 5407 (510) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 6] */
  10816609,  /* 5408 (510) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 9] */
  1544626 ,  /* 5409 (510) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27821101,  /* 5410 (510) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 9, 9] */
  17133082,  /* 5411 (510) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  10697071,  /* 5412 (510) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 9, 3] */
  18553450,  /* 5413 (510) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  10822393,  /* 5414 (510) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 9] */
  17013784,  /* 5415 (510) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 0] */
  26398213,  /* 5416 (510) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 3, 9] */
  17131846,  /* 5417 (510) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  27815329,  /* 5418 (510) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 9, 9] */
  18435280,  /* 5419 (510) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  9405169 ,  /* 5420 (510) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 9] */
  16484   ,  /* 5421 (511) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 5422 (511) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 5423 (511) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 5424 (511) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 5425 (511) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 5426 (511) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 5427 (511) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 5428 (511) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 5429 (511) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9184    ,  /* 5430 (511) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5431 (511) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9040    ,  /* 5432 (511) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5433 (511) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10660   ,  /* 5434 (511) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5435 (511) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10480   ,  /* 5436 (511) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5437 (511) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  3256    ,  /* 5438 (511) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  7672    ,  /* 5439 (511) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  16372   ,  /* 5440 (511) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  7636    ,  /* 5441 (511) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  16432   ,  /* 5442 (511) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  11992   ,  /* 5443 (511) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  3304    ,  /* 5444 (511) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  16484   ,  /* 5445 (512) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 5446 (512) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 5447 (512) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 5448 (512) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 5449 (512) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 5450 (512) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 5451 (512) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 5452 (512) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 5453 (512) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9184    ,  /* 5454 (512) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5455 (512) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9040    ,  /* 5456 (512) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5457 (512) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10660   ,  /* 5458 (512) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5459 (512) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10480   ,  /* 5460 (512) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5461 (512) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  3256    ,  /* 5462 (512) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  7672    ,  /* 5463 (512) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  16372   ,  /* 5464 (512) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  7636    ,  /* 5465 (512) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  16432   ,  /* 5466 (512) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  11992   ,  /* 5467 (512) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  3304    ,  /* 5468 (512) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  1551758 ,  /* 5469 (512) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1547058 ,  /* 5470 (512) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 5471 (512) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1543170 ,  /* 5472 (512) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 5473 (512) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1547546 ,  /* 5474 (512) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 5475 (512) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1542686 ,  /* 5476 (512) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1546090 ,  /* 5477 (512) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544458 ,  /* 5478 (512) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 6] */
  1545598 ,  /* 5479 (512) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544314 ,  /* 5480 (512) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 6, 6] */
  1544626 ,  /* 5481 (512) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545934 ,  /* 5482 (512) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  1544146 ,  /* 5483 (512) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1545754 ,  /* 5484 (512) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  1547338 ,  /* 5485 (512) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1538530 ,  /* 5486 (512) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 6] */
  1542946 ,  /* 5487 (512) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1551646 ,  /* 5488 (512) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  1542910 ,  /* 5489 (512) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1551706 ,  /* 5490 (512) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  1547266 ,  /* 5491 (512) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1538578 ,  /* 5492 (512) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 6] */
  17140694,  /* 5493 (512) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17135994,  /* 5494 (512) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 5495 (512) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17132106,  /* 5496 (512) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 5497 (512) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17136482,  /* 5498 (512) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 5499 (512) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17131622,  /* 5500 (512) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135026,  /* 5501 (512) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133394,  /* 5502 (512) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 6] */
  17134534,  /* 5503 (512) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133250,  /* 5504 (512) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 0, 6] */
  17133562,  /* 5505 (512) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134870,  /* 5506 (512) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 6] */
  17133082,  /* 5507 (512) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17134690,  /* 5508 (512) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 6] */
  17136274,  /* 5509 (512) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17127466,  /* 5510 (512) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 6] */
  17131882,  /* 5511 (512) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17140582,  /* 5512 (512) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 6] */
  17131846,  /* 5513 (512) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17140642,  /* 5514 (512) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 6] */
  17136202,  /* 5515 (512) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17127514,  /* 5516 (512) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 6] */
  18439772,  /* 5517 (512) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18435072,  /* 5518 (512) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18426488,  /* 5519 (512) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18431184,  /* 5520 (512) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 5521 (512) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 5522 (512) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 5523 (512) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 5524 (512) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18434104,  /* 5525 (512) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432472,  /* 5526 (512) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 0] */
  18433612,  /* 5527 (512) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432328,  /* 5528 (512) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 0] */
  18432640,  /* 5529 (512) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433948,  /* 5530 (512) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 0] */
  18432160,  /* 5531 (512) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18433768,  /* 5532 (512) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 0] */
  18435352,  /* 5533 (512) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18426544,  /* 5534 (512) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 0] */
  18430960,  /* 5535 (512) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18439660,  /* 5536 (512) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 0] */
  18430924,  /* 5537 (512) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18439720,  /* 5538 (512) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 0] */
  18435280,  /* 5539 (512) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18426592,  /* 5540 (512) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 0] */
  16484   ,  /* 5541 (513) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 5542 (513) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 5543 (513) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 5544 (513) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 5545 (513) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 5546 (513) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 5547 (513) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 5548 (513) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 5549 (513) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9184    ,  /* 5550 (513) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5551 (513) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9040    ,  /* 5552 (513) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5553 (513) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10660   ,  /* 5554 (513) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5555 (513) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10480   ,  /* 5556 (513) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5557 (513) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  3256    ,  /* 5558 (513) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  7672    ,  /* 5559 (513) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  16372   ,  /* 5560 (513) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  7636    ,  /* 5561 (513) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  16432   ,  /* 5562 (513) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  11992   ,  /* 5563 (513) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  3304    ,  /* 5564 (513) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  18557870,  /* 5565 (513) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 5566 (513) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 5567 (513) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 5568 (513) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 5569 (513) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 5570 (513) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 5571 (513) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 5572 (513) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18552202,  /* 5573 (513) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550570,  /* 5574 (513) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18551710,  /* 5575 (513) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550426,  /* 5576 (513) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18550738,  /* 5577 (513) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18552046,  /* 5578 (513) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18550258,  /* 5579 (513) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18551866,  /* 5580 (513) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18553450,  /* 5581 (513) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18544642,  /* 5582 (513) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18549058,  /* 5583 (513) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18557758,  /* 5584 (513) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18549022,  /* 5585 (513) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18557818,  /* 5586 (513) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18553378,  /* 5587 (513) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18544690,  /* 5588 (513) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  16484   ,  /* 5589 (514) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18553170,  /* 5590 (514) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3200    ,  /* 5591 (514) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18549282,  /* 5592 (514) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16320   ,  /* 5593 (514) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18553658,  /* 5594 (514) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3360    ,  /* 5595 (514) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18548798,  /* 5596 (514) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  10816   ,  /* 5597 (514) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18550570,  /* 5598 (514) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  10324   ,  /* 5599 (514) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18550426,  /* 5600 (514) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  9352    ,  /* 5601 (514) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18552046,  /* 5602 (514) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  8872    ,  /* 5603 (514) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18551866,  /* 5604 (514) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  12064   ,  /* 5605 (514) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  18544642,  /* 5606 (514) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  7672    ,  /* 5607 (514) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  18557758,  /* 5608 (514) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  7636    ,  /* 5609 (514) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  18557818,  /* 5610 (514) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  11992   ,  /* 5611 (514) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18544690,  /* 5612 (514) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  16484   ,  /* 5613 (515) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  129882  ,  /* 5614 (515) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3200    ,  /* 5615 (515) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  125994  ,  /* 5616 (515) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 5617 (515) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  130370  ,  /* 5618 (515) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 5619 (515) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  125510  ,  /* 5620 (515) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  10816   ,  /* 5621 (515) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  17015296,  /* 5622 (515) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 0] */
  10324   ,  /* 5623 (515) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  17015152,  /* 5624 (515) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 0, 0] */
  9352    ,  /* 5625 (515) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  17016772,  /* 5626 (515) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 0] */
  8872    ,  /* 5627 (515) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  17016592,  /* 5628 (515) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 0] */
  12064   ,  /* 5629 (515) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  1420432 ,  /* 5630 (515) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 0] */
  7672    ,  /* 5631 (515) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  1433548 ,  /* 5632 (515) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 0] */
  7636    ,  /* 5633 (515) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  1433608 ,  /* 5634 (515) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 0] */
  11992   ,  /* 5635 (515) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  1420480 ,  /* 5636 (515) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 0] */
  1551758 ,  /* 5637 (515) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1428960 ,  /* 5638 (515) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1538474 ,  /* 5639 (515) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1425072 ,  /* 5640 (515) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 5641 (515) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1429448 ,  /* 5642 (515) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1538634 ,  /* 5643 (515) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1424588 ,  /* 5644 (515) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1546090 ,  /* 5645 (515) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  18550570,  /* 5646 (515) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  1545598 ,  /* 5647 (515) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  18550426,  /* 5648 (515) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  1544626 ,  /* 5649 (515) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  18552046,  /* 5650 (515) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  1544146 ,  /* 5651 (515) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  18551866,  /* 5652 (515) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  1547338 ,  /* 5653 (515) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  121354  ,  /* 5654 (515) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 6] */
  1542946 ,  /* 5655 (515) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  134470  ,  /* 5656 (515) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 6] */
  1542910 ,  /* 5657 (515) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  134530  ,  /* 5658 (515) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 6] */
  1547266 ,  /* 5659 (515) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  121402  ,  /* 5660 (515) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 6] */
  17140694,  /* 5661 (515) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17017896,  /* 5662 (515) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17127410,  /* 5663 (515) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17014008,  /* 5664 (515) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17140530,  /* 5665 (515) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17018384,  /* 5666 (515) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17127570,  /* 5667 (515) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17013524,  /* 5668 (515) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17135026,  /* 5669 (515) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  127282  ,  /* 5670 (515) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 6] */
  17134534,  /* 5671 (515) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  127138  ,  /* 5672 (515) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 6] */
  17133562,  /* 5673 (515) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  128758  ,  /* 5674 (515) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 6] */
  17133082,  /* 5675 (515) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  128578  ,  /* 5676 (515) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 6] */
  17136274,  /* 5677 (515) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  18544642,  /* 5678 (515) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  17131882,  /* 5679 (515) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  18557758,  /* 5680 (515) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  17131846,  /* 5681 (515) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  18557818,  /* 5682 (515) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  17136202,  /* 5683 (515) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  18544690,  /* 5684 (515) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  18439772,  /* 5685 (515) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18553170,  /* 5686 (515) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18426488,  /* 5687 (515) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18549282,  /* 5688 (515) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 5689 (515) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18553658,  /* 5690 (515) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 5691 (515) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18548798,  /* 5692 (515) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18434104,  /* 5693 (515) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  1426360 ,  /* 5694 (515) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 0] */
  18433612,  /* 5695 (515) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1426216 ,  /* 5696 (515) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 6, 0] */
  18432640,  /* 5697 (515) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1427836 ,  /* 5698 (515) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 0] */
  18432160,  /* 5699 (515) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  1427656 ,  /* 5700 (515) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 0] */
  18435352,  /* 5701 (515) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  17009368,  /* 5702 (515) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 0] */
  18430960,  /* 5703 (515) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  17022484,  /* 5704 (515) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 0] */
  18430924,  /* 5705 (515) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  17022544,  /* 5706 (515) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 0] */
  18435280,  /* 5707 (515) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  17009416,  /* 5708 (515) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 0] */
  16484   ,  /* 5709 (516) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  10699653,  /* 5710 (516) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 3] */
  1420376 ,  /* 5711 (516) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  27701877,  /* 5712 (516) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  134418  ,  /* 5713 (516) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  9282965 ,  /* 5714 (516) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 3] */
  1538634 ,  /* 5715 (516) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  26284217,  /* 5716 (516) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 3] */
  10816   ,  /* 5717 (516) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9397975 ,  /* 5718 (516) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  3, 3, 9] */
  128422  ,  /* 5719 (516) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 6] */
  10815007,  /* 5720 (516) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  3, 9, 9] */
  17015464,  /* 5721 (516) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 0] */
  9281353 ,  /* 5722 (516) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  3, 3, 3] */
  17133082,  /* 5723 (516) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  10698349,  /* 5724 (516) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  3, 9, 3] */
  12064   ,  /* 5725 (516) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  26280061,  /* 5726 (516) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 3] */
  17013784,  /* 5727 (516) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 0] */
  26411275,  /* 5728 (516) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 3, 9] */
  1424812 ,  /* 5729 (516) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 0] */
  9287125 ,  /* 5730 (516) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 3, 3] */
  18435280,  /* 5731 (516) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  9392095 ,  /* 5732 (516) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 9] */
  18557870,  /* 5733 (516) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  26406687,  /* 5734 (516) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  17127410,  /* 5735 (516) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  9396687 ,  /* 5736 (516) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 9] */
  18439608,  /* 5737 (516) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  27824351,  /* 5738 (516) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 9] */
  17009472,  /* 5739 (516) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  10813379,  /* 5740 (516) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 9] */
  18552202,  /* 5741 (516) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  27703165,  /* 5742 (516) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  9, 9, 3] */
  18433612,  /* 5743 (516) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26285845,  /* 5744 (516) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  9, 3, 3] */
  1544626 ,  /* 5745 (516) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27822739,  /* 5746 (516) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  9, 9, 9] */
  1426048 ,  /* 5747 (516) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 0] */
  26405383,  /* 5748 (516) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  9, 3, 9] */
  18553450,  /* 5749 (516) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  10809223,  /* 5750 (516) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 9] */
  1542946 ,  /* 5751 (516) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  10704241,  /* 5752 (516) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 9, 3] */
  17131846,  /* 5753 (516) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  27828511,  /* 5754 (516) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 9, 9] */
  130090  ,  /* 5755 (516) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 6] */
  27697285,  /* 5756 (516) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 3] */
  16484   ,  /* 5757 (517) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 5758 (517) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 5759 (517) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 5760 (517) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 5761 (517) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 5762 (517) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 5763 (517) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 5764 (517) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 5765 (517) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 5766 (517) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 5767 (517) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 5768 (517) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 5769 (517) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 5770 (517) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 5771 (517) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 5772 (517) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 5773 (517) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 5774 (517) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10498   ,  /* 5775 (517) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9184    ,  /* 5776 (517) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5777 (517) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 5778 (517) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10642   ,  /* 5779 (517) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9040    ,  /* 5780 (517) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5781 (517) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 5782 (517) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9022    ,  /* 5783 (517) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10660   ,  /* 5784 (517) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5785 (517) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 5786 (517) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9202    ,  /* 5787 (517) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  10480   ,  /* 5788 (517) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5789 (517) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 5790 (517) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16426   ,  /* 5791 (517) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  3256    ,  /* 5792 (517) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  7672    ,  /* 5793 (517) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 5794 (517) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3310    ,  /* 5795 (517) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  16372   ,  /* 5796 (517) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  7636    ,  /* 5797 (517) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 5798 (517) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3250    ,  /* 5799 (517) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  16432   ,  /* 5800 (517) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  11992   ,  /* 5801 (517) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 5802 (517) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16378   ,  /* 5803 (517) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  3304    ,  /* 5804 (517) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  16484   ,  /* 5805 (518) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7898    ,  /* 5806 (518) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3200    ,  /* 5807 (518) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  11786   ,  /* 5808 (518) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  16320   ,  /* 5809 (518) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7410    ,  /* 5810 (518) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3360    ,  /* 5811 (518) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  12270   ,  /* 5812 (518) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  10816   ,  /* 5813 (518) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10498   ,  /* 5814 (518) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  10324   ,  /* 5815 (518) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10642   ,  /* 5816 (518) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9352    ,  /* 5817 (518) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9022    ,  /* 5818 (518) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  8872    ,  /* 5819 (518) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9202    ,  /* 5820 (518) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  12064   ,  /* 5821 (518) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16426   ,  /* 5822 (518) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  7672    ,  /* 5823 (518) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3310    ,  /* 5824 (518) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  7636    ,  /* 5825 (518) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3250    ,  /* 5826 (518) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  11992   ,  /* 5827 (518) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16378   ,  /* 5828 (518) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  18544584,  /* 5829 (518) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18553170,  /* 5830 (518) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557868,  /* 5831 (518) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18549282,  /* 5832 (518) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 5833 (518) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18553658,  /* 5834 (518) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18557708,  /* 5835 (518) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18548798,  /* 5836 (518) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18550252,  /* 5837 (518) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550570,  /* 5838 (518) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18550744,  /* 5839 (518) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550426,  /* 5840 (518) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18551716,  /* 5841 (518) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18552046,  /* 5842 (518) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18552196,  /* 5843 (518) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18551866,  /* 5844 (518) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18549004,  /* 5845 (518) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18544642,  /* 5846 (518) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18553396,  /* 5847 (518) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18557758,  /* 5848 (518) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18553432,  /* 5849 (518) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18557818,  /* 5850 (518) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18549076,  /* 5851 (518) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18544690,  /* 5852 (518) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  16484   ,  /* 5853 (519) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 5854 (519) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  17014010,  /* 5855 (519) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 5856 (519) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  18426488,  /* 5857 (519) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 5858 (519) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  1428962 ,  /* 5859 (519) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 5860 (519) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 5861 (519) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 5862 (519) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  18548796,  /* 5863 (519) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 5864 (519) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  17127570,  /* 5865 (519) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 5866 (519) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  130368  ,  /* 5867 (519) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 5868 (519) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  10816   ,  /* 5869 (519) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 5870 (519) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  1427674 ,  /* 5871 (519) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 6, 0] */
  1426360 ,  /* 5872 (519) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 0] */
  1545598 ,  /* 5873 (519) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544632 ,  /* 5874 (519) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  128740  ,  /* 5875 (519) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 6] */
  127138  ,  /* 5876 (519) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 6] */
  17133562,  /* 5877 (519) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134540,  /* 5878 (519) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  18550408,  /* 5879 (519) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18552046,  /* 5880 (519) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18432160,  /* 5881 (519) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18434098,  /* 5882 (519) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  17015314,  /* 5883 (519) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 0, 0] */
  17016592,  /* 5884 (519) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 0] */
  12064   ,  /* 5885 (519) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 5886 (519) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  134524  ,  /* 5887 (519) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 6] */
  121354  ,  /* 5888 (519) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 6] */
  17131882,  /* 5889 (519) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136220,  /* 5890 (519) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17009422,  /* 5891 (519) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 0] */
  17022484,  /* 5892 (519) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 0] */
  18430924,  /* 5893 (519) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435334,  /* 5894 (519) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18544636,  /* 5895 (519) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18557818,  /* 5896 (519) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  1547266 ,  /* 5897 (519) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542964 ,  /* 5898 (519) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1433554 ,  /* 5899 (519) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 0] */
  1420480 ,  /* 5900 (519) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 0] */
  16484   ,  /* 5901 (520) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 5902 (520) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18549284,  /* 5903 (520) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 5904 (520) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3200    ,  /* 5905 (520) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 5906 (520) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18553172,  /* 5907 (520) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 5908 (520) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  16320   ,  /* 5909 (520) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 5910 (520) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18548796,  /* 5911 (520) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 5912 (520) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3360    ,  /* 5913 (520) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 5914 (520) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553656,  /* 5915 (520) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 5916 (520) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  10816   ,  /* 5917 (520) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 5918 (520) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18551884,  /* 5919 (520) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18550570,  /* 5920 (520) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  10324   ,  /* 5921 (520) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 5922 (520) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18552028,  /* 5923 (520) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18550426,  /* 5924 (520) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  9352    ,  /* 5925 (520) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 5926 (520) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18550408,  /* 5927 (520) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18552046,  /* 5928 (520) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  8872    ,  /* 5929 (520) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 5930 (520) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18550588,  /* 5931 (520) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18551866,  /* 5932 (520) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  12064   ,  /* 5933 (520) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 5934 (520) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18557812,  /* 5935 (520) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  18544642,  /* 5936 (520) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  7672    ,  /* 5937 (520) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 5938 (520) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  18544696,  /* 5939 (520) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18557758,  /* 5940 (520) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  7636    ,  /* 5941 (520) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 5942 (520) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  18544636,  /* 5943 (520) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18557818,  /* 5944 (520) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  11992   ,  /* 5945 (520) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 5946 (520) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  18557764,  /* 5947 (520) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18544690,  /* 5948 (520) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  16484   ,  /* 5949 (521) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  18549284,  /* 5950 (521) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  3200    ,  /* 5951 (521) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  18553172,  /* 5952 (521) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  16320   ,  /* 5953 (521) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  18548796,  /* 5954 (521) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  3360    ,  /* 5955 (521) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  18553656,  /* 5956 (521) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  10816   ,  /* 5957 (521) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18551884,  /* 5958 (521) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  10324   ,  /* 5959 (521) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18552028,  /* 5960 (521) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  9352    ,  /* 5961 (521) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  18550408,  /* 5962 (521) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  8872    ,  /* 5963 (521) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  18550588,  /* 5964 (521) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  12064   ,  /* 5965 (521) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  18557812,  /* 5966 (521) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  7672    ,  /* 5967 (521) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  18544696,  /* 5968 (521) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  7636    ,  /* 5969 (521) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  18544636,  /* 5970 (521) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  11992   ,  /* 5971 (521) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18557764,  /* 5972 (521) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18544584,  /* 5973 (521) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  11784   ,  /* 5974 (521) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18557868,  /* 5975 (521) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  7896    ,  /* 5976 (521) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  18544748,  /* 5977 (521) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  12272   ,  /* 5978 (521) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18557708,  /* 5979 (521) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  7412    ,  /* 5980 (521) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  18550252,  /* 5981 (521) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  9184    ,  /* 5982 (521) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  18550744,  /* 5983 (521) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  9040    ,  /* 5984 (521) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  18551716,  /* 5985 (521) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  10660   ,  /* 5986 (521) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  18552196,  /* 5987 (521) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  10480   ,  /* 5988 (521) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  18549004,  /* 5989 (521) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  3256    ,  /* 5990 (521) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  18553396,  /* 5991 (521) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  16372   ,  /* 5992 (521) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  18553432,  /* 5993 (521) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  16432   ,  /* 5994 (521) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  18549076,  /* 5995 (521) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  3304    ,  /* 5996 (521) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  16484   ,  /* 5997 (522) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 5998 (522) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  1543172 ,  /* 5999 (522) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1547058 ,  /* 6000 (522) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  18426488,  /* 6001 (522) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 6002 (522) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  17135996,  /* 6003 (522) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17132106,  /* 6004 (522) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  1551594 ,  /* 6005 (522) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 6006 (522) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  7410    ,  /* 6007 (522) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 6008 (522) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  17127570,  /* 6009 (522) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 6010 (522) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  18435558,  /* 6011 (522) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 6012 (522) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  10816   ,  /* 6013 (522) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 6014 (522) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  17134708,  /* 6015 (522) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 0, 6] */
  17133394,  /* 6016 (522) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 6] */
  1545598 ,  /* 6017 (522) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544632 ,  /* 6018 (522) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  18433930,  /* 6019 (522) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 0] */
  18432328,  /* 6020 (522) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 0] */
  17133562,  /* 6021 (522) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134540,  /* 6022 (522) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  9022    ,  /* 6023 (522) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10660   ,  /* 6024 (522) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  18432160,  /* 6025 (522) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18434098,  /* 6026 (522) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1544476 ,  /* 6027 (522) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  1545754 ,  /* 6028 (522) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  12064   ,  /* 6029 (522) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 6030 (522) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  18439714,  /* 6031 (522) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 0] */
  18426544,  /* 6032 (522) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 0] */
  17131882,  /* 6033 (522) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136220,  /* 6034 (522) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  1538584 ,  /* 6035 (522) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  1551646 ,  /* 6036 (522) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  18430924,  /* 6037 (522) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435334,  /* 6038 (522) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  3250    ,  /* 6039 (522) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  16432   ,  /* 6040 (522) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  1547266 ,  /* 6041 (522) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542964 ,  /* 6042 (522) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  17140588,  /* 6043 (522) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 6] */
  17127514,  /* 6044 (522) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 6] */
  16484   ,  /* 6045 (523) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 6046 (523) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 6047 (523) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 6048 (523) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 6049 (523) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 6050 (523) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 6051 (523) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 6052 (523) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 6053 (523) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 6054 (523) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 6055 (523) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 6056 (523) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 6057 (523) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 6058 (523) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 6059 (523) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 6060 (523) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 6061 (523) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 6062 (523) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10498   ,  /* 6063 (523) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9184    ,  /* 6064 (523) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  10324   ,  /* 6065 (523) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 6066 (523) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10642   ,  /* 6067 (523) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9040    ,  /* 6068 (523) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9352    ,  /* 6069 (523) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 6070 (523) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9022    ,  /* 6071 (523) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10660   ,  /* 6072 (523) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  8872    ,  /* 6073 (523) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 6074 (523) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9202    ,  /* 6075 (523) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  10480   ,  /* 6076 (523) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  12064   ,  /* 6077 (523) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 6078 (523) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16426   ,  /* 6079 (523) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  3256    ,  /* 6080 (523) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  7672    ,  /* 6081 (523) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 6082 (523) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3310    ,  /* 6083 (523) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  16372   ,  /* 6084 (523) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  7636    ,  /* 6085 (523) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 6086 (523) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3250    ,  /* 6087 (523) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  16432   ,  /* 6088 (523) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  11992   ,  /* 6089 (523) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 6090 (523) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16378   ,  /* 6091 (523) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  3304    ,  /* 6092 (523) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  1551758 ,  /* 6093 (523) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 6094 (523) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1543172 ,  /* 6095 (523) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1547058 ,  /* 6096 (523) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1538474 ,  /* 6097 (523) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 6098 (523) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1547060 ,  /* 6099 (523) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1543170 ,  /* 6100 (523) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1551594 ,  /* 6101 (523) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 6102 (523) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1542684 ,  /* 6103 (523) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1547546 ,  /* 6104 (523) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1538634 ,  /* 6105 (523) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 6106 (523) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1547544 ,  /* 6107 (523) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1542686 ,  /* 6108 (523) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  1546090 ,  /* 6109 (523) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544140 ,  /* 6110 (523) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545772 ,  /* 6111 (523) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 6, 6] */
  1544458 ,  /* 6112 (523) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 6] */
  1545598 ,  /* 6113 (523) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544632 ,  /* 6114 (523) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1545916 ,  /* 6115 (523) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 6, 6] */
  1544314 ,  /* 6116 (523) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 6, 6] */
  1544626 ,  /* 6117 (523) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545604 ,  /* 6118 (523) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544296 ,  /* 6119 (523) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  1545934 ,  /* 6120 (523) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  1544146 ,  /* 6121 (523) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1546084 ,  /* 6122 (523) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544476 ,  /* 6123 (523) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  1545754 ,  /* 6124 (523) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  1547338 ,  /* 6125 (523) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542892 ,  /* 6126 (523) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1551700 ,  /* 6127 (523) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 6] */
  1538530 ,  /* 6128 (523) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 6] */
  1542946 ,  /* 6129 (523) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547284 ,  /* 6130 (523) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1538584 ,  /* 6131 (523) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  1551646 ,  /* 6132 (523) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  1542910 ,  /* 6133 (523) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1547320 ,  /* 6134 (523) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1538524 ,  /* 6135 (523) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  1551706 ,  /* 6136 (523) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  1547266 ,  /* 6137 (523) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542964 ,  /* 6138 (523) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1551652 ,  /* 6139 (523) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 6] */
  1538578 ,  /* 6140 (523) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 6] */
  17140694,  /* 6141 (523) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 6142 (523) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17132108,  /* 6143 (523) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135994,  /* 6144 (523) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17127410,  /* 6145 (523) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 6146 (523) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17135996,  /* 6147 (523) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17132106,  /* 6148 (523) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17140530,  /* 6149 (523) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 6150 (523) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17131620,  /* 6151 (523) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17136482,  /* 6152 (523) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17127570,  /* 6153 (523) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 6154 (523) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17136480,  /* 6155 (523) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17131622,  /* 6156 (523) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  17135026,  /* 6157 (523) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133076,  /* 6158 (523) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134708,  /* 6159 (523) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 0, 6] */
  17133394,  /* 6160 (523) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 6] */
  17134534,  /* 6161 (523) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133568,  /* 6162 (523) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17134852,  /* 6163 (523) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 0, 6] */
  17133250,  /* 6164 (523) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 0, 6] */
  17133562,  /* 6165 (523) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134540,  /* 6166 (523) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133232,  /* 6167 (523) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 0, 6] */
  17134870,  /* 6168 (523) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 6] */
  17133082,  /* 6169 (523) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17135020,  /* 6170 (523) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133412,  /* 6171 (523) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 0, 6] */
  17134690,  /* 6172 (523) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 6] */
  17136274,  /* 6173 (523) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131828,  /* 6174 (523) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17140636,  /* 6175 (523) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 6] */
  17127466,  /* 6176 (523) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 6] */
  17131882,  /* 6177 (523) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136220,  /* 6178 (523) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17127520,  /* 6179 (523) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 6] */
  17140582,  /* 6180 (523) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 6] */
  17131846,  /* 6181 (523) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136256,  /* 6182 (523) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17127460,  /* 6183 (523) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 6] */
  17140642,  /* 6184 (523) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 6] */
  17136202,  /* 6185 (523) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17131900,  /* 6186 (523) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17140588,  /* 6187 (523) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 6] */
  17127514,  /* 6188 (523) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 6] */
  18439772,  /* 6189 (523) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 6190 (523) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18431186,  /* 6191 (523) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18435072,  /* 6192 (523) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18426488,  /* 6193 (523) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 6194 (523) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18435074,  /* 6195 (523) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18431184,  /* 6196 (523) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18439608,  /* 6197 (523) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 6198 (523) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18430698,  /* 6199 (523) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 6200 (523) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18426648,  /* 6201 (523) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 6202 (523) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18435558,  /* 6203 (523) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18430700,  /* 6204 (523) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  18434104,  /* 6205 (523) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432154,  /* 6206 (523) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433786,  /* 6207 (523) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 0] */
  18432472,  /* 6208 (523) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 0] */
  18433612,  /* 6209 (523) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432646,  /* 6210 (523) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18433930,  /* 6211 (523) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 0] */
  18432328,  /* 6212 (523) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 0] */
  18432640,  /* 6213 (523) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433618,  /* 6214 (523) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432310,  /* 6215 (523) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 0] */
  18433948,  /* 6216 (523) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 0] */
  18432160,  /* 6217 (523) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18434098,  /* 6218 (523) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432490,  /* 6219 (523) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 0] */
  18433768,  /* 6220 (523) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 0] */
  18435352,  /* 6221 (523) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430906,  /* 6222 (523) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18439714,  /* 6223 (523) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 0] */
  18426544,  /* 6224 (523) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 0] */
  18430960,  /* 6225 (523) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18435298,  /* 6226 (523) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18426598,  /* 6227 (523) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 0] */
  18439660,  /* 6228 (523) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 0] */
  18430924,  /* 6229 (523) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435334,  /* 6230 (523) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18426538,  /* 6231 (523) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 0] */
  18439720,  /* 6232 (523) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 0] */
  18435280,  /* 6233 (523) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430978,  /* 6234 (523) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18439666,  /* 6235 (523) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 0] */
  18426592,  /* 6236 (523) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 0] */
  16484   ,  /* 6237 (524) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 6238 (524) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  125996  ,  /* 6239 (524) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  129882  ,  /* 6240 (524) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  3200    ,  /* 6241 (524) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 6242 (524) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  129884  ,  /* 6243 (524) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  125994  ,  /* 6244 (524) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  16320   ,  /* 6245 (524) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 6246 (524) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  125508  ,  /* 6247 (524) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 6248 (524) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  3360    ,  /* 6249 (524) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 6250 (524) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  130368  ,  /* 6251 (524) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  125510  ,  /* 6252 (524) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  10816   ,  /* 6253 (524) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 6254 (524) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  17016610,  /* 6255 (524) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 0, 0] */
  17015296,  /* 6256 (524) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 0] */
  10324   ,  /* 6257 (524) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 6258 (524) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  17016754,  /* 6259 (524) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 0, 0] */
  17015152,  /* 6260 (524) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 0, 0] */
  9352    ,  /* 6261 (524) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 6262 (524) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  17015134,  /* 6263 (524) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 0, 0] */
  17016772,  /* 6264 (524) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 0] */
  8872    ,  /* 6265 (524) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 6266 (524) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  17015314,  /* 6267 (524) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 0, 0] */
  17016592,  /* 6268 (524) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 0] */
  12064   ,  /* 6269 (524) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 6270 (524) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  1433602 ,  /* 6271 (524) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 0] */
  1420432 ,  /* 6272 (524) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 0] */
  7672    ,  /* 6273 (524) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 6274 (524) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  1420486 ,  /* 6275 (524) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 0] */
  1433548 ,  /* 6276 (524) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 0] */
  7636    ,  /* 6277 (524) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 6278 (524) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  1420426 ,  /* 6279 (524) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 0] */
  1433608 ,  /* 6280 (524) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 0] */
  11992   ,  /* 6281 (524) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 6282 (524) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  1433554 ,  /* 6283 (524) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 0] */
  1420480 ,  /* 6284 (524) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 0] */
  1551758 ,  /* 6285 (524) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 6286 (524) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1425074 ,  /* 6287 (524) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1428960 ,  /* 6288 (524) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1538474 ,  /* 6289 (524) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1551756 ,  /* 6290 (524) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1428962 ,  /* 6291 (524) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1425072 ,  /* 6292 (524) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1551594 ,  /* 6293 (524) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  1538636 ,  /* 6294 (524) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1424586 ,  /* 6295 (524) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1429448 ,  /* 6296 (524) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1538634 ,  /* 6297 (524) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 6298 (524) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  1429446 ,  /* 6299 (524) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1424588 ,  /* 6300 (524) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  1546090 ,  /* 6301 (524) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544140 ,  /* 6302 (524) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  18551884,  /* 6303 (524) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18550570,  /* 6304 (524) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  1545598 ,  /* 6305 (524) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1544632 ,  /* 6306 (524) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  18552028,  /* 6307 (524) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18550426,  /* 6308 (524) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  1544626 ,  /* 6309 (524) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545604 ,  /* 6310 (524) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  18550408,  /* 6311 (524) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18552046,  /* 6312 (524) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  1544146 ,  /* 6313 (524) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1546084 ,  /* 6314 (524) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  18550588,  /* 6315 (524) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18551866,  /* 6316 (524) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  1547338 ,  /* 6317 (524) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542892 ,  /* 6318 (524) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  134524  ,  /* 6319 (524) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 6] */
  121354  ,  /* 6320 (524) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 6] */
  1542946 ,  /* 6321 (524) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547284 ,  /* 6322 (524) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  121408  ,  /* 6323 (524) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 6] */
  134470  ,  /* 6324 (524) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 6] */
  1542910 ,  /* 6325 (524) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  1547320 ,  /* 6326 (524) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  121348  ,  /* 6327 (524) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 6] */
  134530  ,  /* 6328 (524) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 6] */
  1547266 ,  /* 6329 (524) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  1542964 ,  /* 6330 (524) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  134476  ,  /* 6331 (524) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 6] */
  121402  ,  /* 6332 (524) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 6] */
  17140694,  /* 6333 (524) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 6334 (524) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17014010,  /* 6335 (524) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17017896,  /* 6336 (524) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17127410,  /* 6337 (524) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 6338 (524) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17017898,  /* 6339 (524) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17014008,  /* 6340 (524) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17140530,  /* 6341 (524) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17127572,  /* 6342 (524) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17013522,  /* 6343 (524) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17018384,  /* 6344 (524) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17127570,  /* 6345 (524) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  17140532,  /* 6346 (524) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17018382,  /* 6347 (524) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17013524,  /* 6348 (524) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 0] */
  17135026,  /* 6349 (524) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133076,  /* 6350 (524) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  128596  ,  /* 6351 (524) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 6] */
  127282  ,  /* 6352 (524) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 6] */
  17134534,  /* 6353 (524) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17133568,  /* 6354 (524) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  128740  ,  /* 6355 (524) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 6] */
  127138  ,  /* 6356 (524) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 6] */
  17133562,  /* 6357 (524) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  17134540,  /* 6358 (524) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  127120  ,  /* 6359 (524) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 6] */
  128758  ,  /* 6360 (524) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 6] */
  17133082,  /* 6361 (524) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17135020,  /* 6362 (524) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  127300  ,  /* 6363 (524) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 6] */
  128578  ,  /* 6364 (524) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 6] */
  17136274,  /* 6365 (524) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131828,  /* 6366 (524) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  18557812,  /* 6367 (524) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  18544642,  /* 6368 (524) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  17131882,  /* 6369 (524) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  17136220,  /* 6370 (524) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  18544696,  /* 6371 (524) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18557758,  /* 6372 (524) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  17131846,  /* 6373 (524) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136256,  /* 6374 (524) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  18544636,  /* 6375 (524) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18557818,  /* 6376 (524) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  17136202,  /* 6377 (524) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  17131900,  /* 6378 (524) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  18557764,  /* 6379 (524) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18544690,  /* 6380 (524) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  18439772,  /* 6381 (524) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 6382 (524) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18549284,  /* 6383 (524) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 6384 (524) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18426488,  /* 6385 (524) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18439770,  /* 6386 (524) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18553172,  /* 6387 (524) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 6388 (524) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18439608,  /* 6389 (524) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 6390 (524) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18548796,  /* 6391 (524) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 6392 (524) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18426648,  /* 6393 (524) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  18439610,  /* 6394 (524) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  18553656,  /* 6395 (524) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 6396 (524) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18434104,  /* 6397 (524) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432154,  /* 6398 (524) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1427674 ,  /* 6399 (524) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 6, 0] */
  1426360 ,  /* 6400 (524) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 0] */
  18433612,  /* 6401 (524) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432646,  /* 6402 (524) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  1427818 ,  /* 6403 (524) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 6, 0] */
  1426216 ,  /* 6404 (524) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 6, 0] */
  18432640,  /* 6405 (524) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18433618,  /* 6406 (524) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  1426198 ,  /* 6407 (524) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 6, 0] */
  1427836 ,  /* 6408 (524) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 0] */
  18432160,  /* 6409 (524) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18434098,  /* 6410 (524) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  1426378 ,  /* 6411 (524) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 6, 0] */
  1427656 ,  /* 6412 (524) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 0] */
  18435352,  /* 6413 (524) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430906,  /* 6414 (524) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  17022538,  /* 6415 (524) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 0] */
  17009368,  /* 6416 (524) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 0] */
  18430960,  /* 6417 (524) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  18435298,  /* 6418 (524) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  17009422,  /* 6419 (524) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 0] */
  17022484,  /* 6420 (524) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 0] */
  18430924,  /* 6421 (524) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  18435334,  /* 6422 (524) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  17009362,  /* 6423 (524) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 0] */
  17022544,  /* 6424 (524) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 0] */
  18435280,  /* 6425 (524) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430978,  /* 6426 (524) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  17022490,  /* 6427 (524) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 0] */
  17009416,  /* 6428 (524) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 0] */
  16484   ,  /* 6429 (525) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  9278591 ,  /* 6430 (525) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 3] */
  1538474 ,  /* 6431 (525) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  26406689,  /* 6432 (525) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 9] */
  16320   ,  /* 6433 (525) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10813377,  /* 6434 (525) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 9] */
  1538634 ,  /* 6435 (525) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  27706251,  /* 6436 (525) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  10816   ,  /* 6437 (525) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9281191 ,  /* 6438 (525) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 3, 3] */
  17134534,  /* 6439 (525) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  27704623,  /* 6440 (525) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 9, 3] */
  9352    ,  /* 6441 (525) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  26403925,  /* 6442 (525) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 3, 9] */
  17133082,  /* 6443 (525) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  10815169,  /* 6444 (525) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 9, 9] */
  12064   ,  /* 6445 (525) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  9287119 ,  /* 6446 (525) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 3] */
  18430960,  /* 6447 (525) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  10809277,  /* 6448 (525) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 9, 9] */
  7636    ,  /* 6449 (525) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  27697231,  /* 6450 (525) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 9, 3] */
  18435280,  /* 6451 (525) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  26411281,  /* 6452 (525) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 9] */
  9273891 ,  /* 6453 (525) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 3, 3] */
  11784   ,  /* 6454 (525) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  10822449,  /* 6455 (525) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 9] */
  17132106,  /* 6456 (525) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 6] */
  9274055 ,  /* 6457 (525) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 3, 3] */
  1547546 ,  /* 6458 (525) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  10822289,  /* 6459 (525) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 9] */
  18430700,  /* 6460 (525) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 0] */
  9279559 ,  /* 6461 (525) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  3, 3, 3] */
  9184    ,  /* 6462 (525) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  26404261,  /* 6463 (525) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  9, 3, 9] */
  18432328,  /* 6464 (525) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 0] */
  9281023 ,  /* 6465 (525) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 3, 3] */
  17134870,  /* 6466 (525) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 6] */
  26405713,  /* 6467 (525) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 3, 9] */
  1545754 ,  /* 6468 (525) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  9278311 ,  /* 6469 (525) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  3, 3, 3] */
  3256    ,  /* 6470 (525) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  27705991,  /* 6471 (525) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 9, 3] */
  1551646 ,  /* 6472 (525) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  9282739 ,  /* 6473 (525) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 3, 3] */
  18439720,  /* 6474 (525) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 0] */
  27701671,  /* 6475 (525) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  9, 9, 3] */
  17127514,  /* 6476 (525) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 6] */
  1551758 ,  /* 6477 (525) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  10813865,  /* 6478 (525) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 9] */
  3200    ,  /* 6479 (525) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  27705767,  /* 6480 (525) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  1551594 ,  /* 6481 (525) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  9278103 ,  /* 6482 (525) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  3360    ,  /* 6483 (525) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  26407173,  /* 6484 (525) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 9] */
  1546090 ,  /* 6485 (525) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  10816465,  /* 6486 (525) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 9, 9] */
  18433612,  /* 6487 (525) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26405545,  /* 6488 (525) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 9] */
  1544626 ,  /* 6489 (525) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27703003,  /* 6490 (525) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 9, 3] */
  18432160,  /* 6491 (525) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  9279895 ,  /* 6492 (525) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 3, 3] */
  1547338 ,  /* 6493 (525) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  10822393,  /* 6494 (525) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 9] */
  17131882,  /* 6495 (525) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  9274003 ,  /* 6496 (525) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 3, 3] */
  1542910 ,  /* 6497 (525) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  26398153,  /* 6498 (525) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 3, 9] */
  17136202,  /* 6499 (525) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  27710359,  /* 6500 (525) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 3] */
  10809165,  /* 6501 (525) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 9, 9] */
  1547058 ,  /* 6502 (525) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  9287175 ,  /* 6503 (525) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 3] */
  18431184,  /* 6504 (525) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 0] */
  10809329,  /* 6505 (525) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 9, 9] */
  12272   ,  /* 6506 (525) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  9287015 ,  /* 6507 (525) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 3] */
  17131622,  /* 6508 (525) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 6] */
  10814833,  /* 6509 (525) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  3, 9, 9] */
  1544458 ,  /* 6510 (525) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 6] */
  27703339,  /* 6511 (525) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  9, 9, 3] */
  17133250,  /* 6512 (525) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 0, 6] */
  10816297,  /* 6513 (525) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 9, 9] */
  18433948,  /* 6514 (525) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 0] */
  27704791,  /* 6515 (525) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 9, 3] */
  10480   ,  /* 6516 (525) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10813585,  /* 6517 (525) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  3, 9, 9] */
  1538530 ,  /* 6518 (525) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 6] */
  26406913,  /* 6519 (525) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 3, 9] */
  16372   ,  /* 6520 (525) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  10818013,  /* 6521 (525) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 9, 9] */
  17140642,  /* 6522 (525) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 6] */
  26402593,  /* 6523 (525) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  9, 3, 9] */
  18426592,  /* 6524 (525) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 0] */
  17140694,  /* 6525 (525) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  26402801,  /* 6526 (525) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  18426488,  /* 6527 (525) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  9282479 ,  /* 6528 (525) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 3] */
  17140530,  /* 6529 (525) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  27701391,  /* 6530 (525) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 3] */
  18426648,  /* 6531 (525) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  10818237,  /* 6532 (525) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  17135026,  /* 6533 (525) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  26405401,  /* 6534 (525) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 3, 9] */
  10324   ,  /* 6535 (525) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10816609,  /* 6536 (525) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 9] */
  17133562,  /* 6537 (525) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  9279715 ,  /* 6538 (525) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 3, 3] */
  8872    ,  /* 6539 (525) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  27703183,  /* 6540 (525) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 9, 3] */
  17136274,  /* 6541 (525) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  26411329,  /* 6542 (525) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 9] */
  1542946 ,  /* 6543 (525) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  27697291,  /* 6544 (525) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 9, 3] */
  17131846,  /* 6545 (525) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  10809217,  /* 6546 (525) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 9, 9] */
  1547266 ,  /* 6547 (525) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  9287071 ,  /* 6548 (525) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 3] */
  26398101,  /* 6549 (525) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 3, 9] */
  17135994,  /* 6550 (525) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  27710463,  /* 6551 (525) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 3] */
  7896    ,  /* 6552 (525) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  26398265,  /* 6553 (525) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 3, 9] */
  18435560,  /* 6554 (525) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  27710303,  /* 6555 (525) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 3] */
  1542686 ,  /* 6556 (525) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 6] */
  26403769,  /* 6557 (525) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  9, 3, 9] */
  17133394,  /* 6558 (525) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 6] */
  9280051 ,  /* 6559 (525) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  3, 3, 3] */
  1544314 ,  /* 6560 (525) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 6, 6] */
  26405233,  /* 6561 (525) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 3, 9] */
  10660   ,  /* 6562 (525) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  9281503 ,  /* 6563 (525) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 3, 3] */
  18433768,  /* 6564 (525) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 0] */
  26402521,  /* 6565 (525) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  9, 3, 9] */
  17127466,  /* 6566 (525) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 6] */
  10817977,  /* 6567 (525) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 9, 9] */
  18439660,  /* 6568 (525) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 0] */
  26406949,  /* 6569 (525) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 3, 9] */
  1551706 ,  /* 6570 (525) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  10813657,  /* 6571 (525) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  3, 9, 9] */
  3304    ,  /* 6572 (525) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  18439772,  /* 6573 (525) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  27701879,  /* 6574 (525) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 3] */
  17127410,  /* 6575 (525) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  10817753,  /* 6576 (525) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 9] */
  18439608,  /* 6577 (525) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  26402313,  /* 6578 (525) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  17127570,  /* 6579 (525) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  9282963 ,  /* 6580 (525) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 3] */
  18434104,  /* 6581 (525) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  27704479,  /* 6582 (525) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 9, 3] */
  1545598 ,  /* 6583 (525) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  9281335 ,  /* 6584 (525) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 3, 3] */
  18432640,  /* 6585 (525) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  10814989,  /* 6586 (525) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 9, 9] */
  1544146 ,  /* 6587 (525) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  26404105,  /* 6588 (525) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 3, 9] */
  18435352,  /* 6589 (525) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  27710407,  /* 6590 (525) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 3] */
  7672    ,  /* 6591 (525) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  26398213,  /* 6592 (525) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 3, 9] */
  18430924,  /* 6593 (525) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  9273943 ,  /* 6594 (525) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 3, 3] */
  11992   ,  /* 6595 (525) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  10822345,  /* 6596 (525) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 9] */
  27697179,  /* 6597 (525) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 9, 3] */
  18435072,  /* 6598 (525) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  26411385,  /* 6599 (525) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 9] */
  1543170 ,  /* 6600 (525) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 6] */
  27697343,  /* 6601 (525) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 9, 3] */
  17136482,  /* 6602 (525) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  26411225,  /* 6603 (525) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 9] */
  7412    ,  /* 6604 (525) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  27702847,  /* 6605 (525) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  9, 9, 3] */
  18432472,  /* 6606 (525) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 0] */
  10815325,  /* 6607 (525) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  3, 9, 9] */
  9040    ,  /* 6608 (525) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  27704311,  /* 6609 (525) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 9, 3] */
  1545934 ,  /* 6610 (525) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  10816777,  /* 6611 (525) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 9, 9] */
  17134690,  /* 6612 (525) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 6] */
  27701599,  /* 6613 (525) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  9, 9, 3] */
  18426544,  /* 6614 (525) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 0] */
  9282703 ,  /* 6615 (525) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 3, 3] */
  17140582,  /* 6616 (525) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 6] */
  27706027,  /* 6617 (525) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 9, 3] */
  16432   ,  /* 6618 (525) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  9278383 ,  /* 6619 (525) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  3, 3, 3] */
  1538578 ,  /* 6620 (525) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 6] */
  16484   ,  /* 6621 (526) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 6622 (526) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  775535  ,  /* 6623 (526) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 3, 3] */
  779421  ,  /* 6624 (526) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 3, 3] */
  26339054,  /* 6625 (526) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 6] */
  26352336,  /* 6626 (526) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 6] */
  25698101,  /* 6627 (526) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 0, 9] */
  25694211,  /* 6628 (526) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 0, 9] */
  783957  ,  /* 6629 (526) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 3, 3] */
  770999  ,  /* 6630 (526) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 3, 3] */
  7410    ,  /* 6631 (526) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 6632 (526) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  25689675,  /* 6633 (526) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 0, 9] */
  25702637,  /* 6634 (526) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 0, 9] */
  26348124,  /* 6635 (526) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 6] */
  26343266,  /* 6636 (526) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 6] */
  10816   ,  /* 6637 (526) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 6638 (526) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  8572603 ,  /* 6639 (526) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 0, 3] */
  8571289 ,  /* 6640 (526) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  3, 0, 3] */
  19201249,  /* 6641 (526) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 9, 3] */
  19200283,  /* 6642 (526) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 9, 3] */
  27645574,  /* 6643 (526) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 9, 0] */
  27643972,  /* 6644 (526) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  9, 9, 0] */
  8571457 ,  /* 6645 (526) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  3, 0, 3] */
  8572435 ,  /* 6646 (526) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 0, 3] */
  9022    ,  /* 6647 (526) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10660   ,  /* 6648 (526) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  27643804,  /* 6649 (526) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  9, 9, 0] */
  27645742,  /* 6650 (526) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 9, 0] */
  19200127,  /* 6651 (526) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 9, 3] */
  19201405,  /* 6652 (526) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 9, 3] */
  12064   ,  /* 6653 (526) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 6654 (526) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  9228070 ,  /* 6655 (526) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 0] */
  9214900 ,  /* 6656 (526) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 0] */
  10105051,  /* 6657 (526) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  3, 6, 9] */
  10109389,  /* 6658 (526) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 6, 9] */
  2306221 ,  /* 6659 (526) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 9, 9] */
  2319283 ,  /* 6660 (526) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 9, 9] */
  9219280 ,  /* 6661 (526) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  3, 3, 0] */
  9223690 ,  /* 6662 (526) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 3, 0] */
  3250    ,  /* 6663 (526) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  16432   ,  /* 6664 (526) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  2314903 ,  /* 6665 (526) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 9, 9] */
  2310601 ,  /* 6666 (526) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 9, 9] */
  10113757,  /* 6667 (526) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 6, 9] */
  10100683,  /* 6668 (526) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 6, 9] */
  1551758 ,  /* 6669 (526) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 6670 (526) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  2310809 ,  /* 6671 (526) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 9, 9] */
  2314695 ,  /* 6672 (526) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 9, 9] */
  27638132,  /* 6673 (526) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 0] */
  27651414,  /* 6674 (526) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 0] */
  26997179,  /* 6675 (526) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 6, 3] */
  26993289,  /* 6676 (526) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 6, 3] */
  2319231 ,  /* 6677 (526) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 9, 9] */
  2306273 ,  /* 6678 (526) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 9, 9] */
  1542684 ,  /* 6679 (526) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 6] */
  1547546 ,  /* 6680 (526) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 6] */
  26988753,  /* 6681 (526) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 6, 3] */
  27001715,  /* 6682 (526) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 6, 3] */
  27647202,  /* 6683 (526) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 0] */
  27642344,  /* 6684 (526) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 0] */
  1546090 ,  /* 6685 (526) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544140 ,  /* 6686 (526) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  10107877,  /* 6687 (526) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 6, 9] */
  10106563,  /* 6688 (526) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  3, 6, 9] */
  17902171,  /* 6689 (526) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 3, 9] */
  17901205,  /* 6690 (526) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 3, 9] */
  26346496,  /* 6691 (526) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 6] */
  26344894,  /* 6692 (526) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  9, 3, 6] */
  10106731,  /* 6693 (526) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  3, 6, 9] */
  10107709,  /* 6694 (526) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 6, 9] */
  1544296 ,  /* 6695 (526) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 6, 6] */
  1545934 ,  /* 6696 (526) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 6] */
  26344726,  /* 6697 (526) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  9, 3, 6] */
  26346664,  /* 6698 (526) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 3, 6] */
  17901049,  /* 6699 (526) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 3, 9] */
  17902327,  /* 6700 (526) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 3, 9] */
  1547338 ,  /* 6701 (526) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542892 ,  /* 6702 (526) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  10763344,  /* 6703 (526) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 6] */
  10750174,  /* 6704 (526) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 6] */
  8569777 ,  /* 6705 (526) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  3, 0, 3] */
  8574115 ,  /* 6706 (526) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 0, 3] */
  770947  ,  /* 6707 (526) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 3, 3] */
  784009  ,  /* 6708 (526) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 3, 3] */
  10754554,  /* 6709 (526) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  3, 9, 6] */
  10758964,  /* 6710 (526) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 9, 6] */
  1538524 ,  /* 6711 (526) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 6] */
  1551706 ,  /* 6712 (526) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 6] */
  779629  ,  /* 6713 (526) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 3, 3] */
  775327  ,  /* 6714 (526) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 3, 3] */
  8578483 ,  /* 6715 (526) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 0, 3] */
  8565409 ,  /* 6716 (526) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 0, 3] */
  17140694,  /* 6717 (526) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 6718 (526) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17899745,  /* 6719 (526) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 3, 9] */
  17903631,  /* 6720 (526) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 3, 9] */
  9214844 ,  /* 6721 (526) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 0] */
  9228126 ,  /* 6722 (526) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 0] */
  8573891 ,  /* 6723 (526) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 0, 3] */
  8570001 ,  /* 6724 (526) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 0, 3] */
  17908167,  /* 6725 (526) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 3, 9] */
  17895209,  /* 6726 (526) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 3, 9] */
  17131620,  /* 6727 (526) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 6] */
  17136482,  /* 6728 (526) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 6] */
  8565465 ,  /* 6729 (526) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 0, 3] */
  8578427 ,  /* 6730 (526) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 0, 3] */
  9223914 ,  /* 6731 (526) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 0] */
  9219056 ,  /* 6732 (526) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 0] */
  17135026,  /* 6733 (526) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133076,  /* 6734 (526) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  25696813,  /* 6735 (526) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 0, 9] */
  25695499,  /* 6736 (526) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  9, 0, 9] */
  2313235 ,  /* 6737 (526) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 9, 9] */
  2312269 ,  /* 6738 (526) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 9, 9] */
  10757560,  /* 6739 (526) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 6] */
  10755958,  /* 6740 (526) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  3, 9, 6] */
  25695667,  /* 6741 (526) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  9, 0, 9] */
  25696645,  /* 6742 (526) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 0, 9] */
  17133232,  /* 6743 (526) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 0, 6] */
  17134870,  /* 6744 (526) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 6] */
  10755790,  /* 6745 (526) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  3, 9, 6] */
  10757728,  /* 6746 (526) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 9, 6] */
  2312113 ,  /* 6747 (526) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 9, 9] */
  2313391 ,  /* 6748 (526) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 9, 9] */
  17136274,  /* 6749 (526) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131828,  /* 6750 (526) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  26352280,  /* 6751 (526) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 6] */
  26339110,  /* 6752 (526) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 6] */
  26993065,  /* 6753 (526) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  9, 6, 3] */
  26997403,  /* 6754 (526) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 6, 3] */
  19194235,  /* 6755 (526) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 9, 3] */
  19207297,  /* 6756 (526) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 9, 3] */
  26343490,  /* 6757 (526) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  9, 3, 6] */
  26347900,  /* 6758 (526) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 3, 6] */
  17127460,  /* 6759 (526) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 6] */
  17140642,  /* 6760 (526) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 6] */
  19202917,  /* 6761 (526) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 9, 3] */
  19198615,  /* 6762 (526) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 9, 3] */
  27001771,  /* 6763 (526) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 6, 3] */
  26988697,  /* 6764 (526) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 6, 3] */
  18439772,  /* 6765 (526) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 6766 (526) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  19198823,  /* 6767 (526) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 9, 3] */
  19202709,  /* 6768 (526) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 9, 3] */
  10750118,  /* 6769 (526) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 6] */
  10763400,  /* 6770 (526) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 6] */
  10109165,  /* 6771 (526) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 6, 9] */
  10105275,  /* 6772 (526) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 6, 9] */
  19207245,  /* 6773 (526) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 9, 3] */
  19194287,  /* 6774 (526) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 9, 3] */
  18430698,  /* 6775 (526) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 0] */
  18435560,  /* 6776 (526) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 0] */
  10100739,  /* 6777 (526) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 6, 9] */
  10113701,  /* 6778 (526) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 6, 9] */
  10759188,  /* 6779 (526) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 6] */
  10754330,  /* 6780 (526) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 6] */
  18434104,  /* 6781 (526) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432154,  /* 6782 (526) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26995891,  /* 6783 (526) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 6, 3] */
  26994577,  /* 6784 (526) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  9, 6, 3] */
  777961  ,  /* 6785 (526) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 3, 3] */
  776995  ,  /* 6786 (526) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 3, 3] */
  9222286 ,  /* 6787 (526) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 3, 0] */
  9220684 ,  /* 6788 (526) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  3, 3, 0] */
  26994745,  /* 6789 (526) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  9, 6, 3] */
  26995723,  /* 6790 (526) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 6, 3] */
  18432310,  /* 6791 (526) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 0] */
  18433948,  /* 6792 (526) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 0] */
  9220516 ,  /* 6793 (526) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  3, 3, 0] */
  9222454 ,  /* 6794 (526) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 3, 0] */
  776839  ,  /* 6795 (526) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 3, 3] */
  778117  ,  /* 6796 (526) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 3, 3] */
  18435352,  /* 6797 (526) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430906,  /* 6798 (526) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  27651358,  /* 6799 (526) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 0] */
  27638188,  /* 6800 (526) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 0] */
  25693987,  /* 6801 (526) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  9, 0, 9] */
  25698325,  /* 6802 (526) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 0, 9] */
  17895157,  /* 6803 (526) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 3, 9] */
  17908219,  /* 6804 (526) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 3, 9] */
  27642568,  /* 6805 (526) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  9, 9, 0] */
  27646978,  /* 6806 (526) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 9, 0] */
  18426538,  /* 6807 (526) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 0] */
  18439720,  /* 6808 (526) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 0] */
  17903839,  /* 6809 (526) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 3, 9] */
  17899537,  /* 6810 (526) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 3, 9] */
  25702693,  /* 6811 (526) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 0, 9] */
  25689619,  /* 6812 (526) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 0, 9] */
  16484   ,  /* 6813 (527) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  9278591 ,  /* 6814 (527) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 3] */
  1538474 ,  /* 6815 (527) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  26406689,  /* 6816 (527) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 9] */
  16320   ,  /* 6817 (527) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10813377,  /* 6818 (527) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 9] */
  1538634 ,  /* 6819 (527) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  27706251,  /* 6820 (527) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  10816   ,  /* 6821 (527) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9281191 ,  /* 6822 (527) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 3, 3] */
  17134534,  /* 6823 (527) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  27704623,  /* 6824 (527) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 9, 3] */
  9352    ,  /* 6825 (527) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  26403925,  /* 6826 (527) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 3, 9] */
  17133082,  /* 6827 (527) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  10815169,  /* 6828 (527) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 9, 9] */
  12064   ,  /* 6829 (527) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  9287119 ,  /* 6830 (527) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 3] */
  18430960,  /* 6831 (527) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 0] */
  10809277,  /* 6832 (527) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 9, 9] */
  7636    ,  /* 6833 (527) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  27697231,  /* 6834 (527) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 9, 3] */
  18435280,  /* 6835 (527) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  26411281,  /* 6836 (527) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 9] */
  9391989 ,  /* 6837 (527) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 3, 9] */
  129882  ,  /* 6838 (527) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  10704351,  /* 6839 (527) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 3] */
  17014008,  /* 6840 (527) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 0, 0] */
  9392153 ,  /* 6841 (527) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 3, 9] */
  1429448 ,  /* 6842 (527) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 0] */
  10704191,  /* 6843 (527) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 3] */
  18548798,  /* 6844 (527) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  9397657 ,  /* 6845 (527) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  3, 3, 9] */
  127282  ,  /* 6846 (527) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 6] */
  26286163,  /* 6847 (527) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  9, 3, 3] */
  18550426,  /* 6848 (527) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  9399121 ,  /* 6849 (527) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 3, 9] */
  17016772,  /* 6850 (527) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 0] */
  26287615,  /* 6851 (527) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 3, 3] */
  1427656 ,  /* 6852 (527) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 6, 0] */
  9396409 ,  /* 6853 (527) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  3, 3, 9] */
  121354  ,  /* 6854 (527) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 6] */
  27824089,  /* 6855 (527) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 9, 9] */
  1433548 ,  /* 6856 (527) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 0] */
  9400837 ,  /* 6857 (527) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 3, 9] */
  18557818,  /* 6858 (527) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  27819769,  /* 6859 (527) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  9, 9, 9] */
  17009416,  /* 6860 (527) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 0, 0] */
  1551758 ,  /* 6861 (527) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  10813865,  /* 6862 (527) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 9] */
  3200    ,  /* 6863 (527) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  27705767,  /* 6864 (527) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  1551594 ,  /* 6865 (527) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  9278103 ,  /* 6866 (527) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  3360    ,  /* 6867 (527) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  26407173,  /* 6868 (527) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 9] */
  1546090 ,  /* 6869 (527) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  10816465,  /* 6870 (527) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 9, 9] */
  18433612,  /* 6871 (527) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  26405545,  /* 6872 (527) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 9] */
  1544626 ,  /* 6873 (527) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27703003,  /* 6874 (527) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 9, 3] */
  18432160,  /* 6875 (527) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 0] */
  9279895 ,  /* 6876 (527) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 3, 3] */
  1547338 ,  /* 6877 (527) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  10822393,  /* 6878 (527) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 9] */
  17131882,  /* 6879 (527) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  9274003 ,  /* 6880 (527) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 3, 3] */
  1542910 ,  /* 6881 (527) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  26398153,  /* 6882 (527) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 3, 9] */
  17136202,  /* 6883 (527) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  27710359,  /* 6884 (527) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 3] */
  10691067,  /* 6885 (527) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  3, 9, 3] */
  1428960 ,  /* 6886 (527) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 0] */
  9405273 ,  /* 6887 (527) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 9] */
  18549282,  /* 6888 (527) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  10691231,  /* 6889 (527) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  3, 9, 3] */
  130370  ,  /* 6890 (527) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  9405113 ,  /* 6891 (527) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 9] */
  17013524,  /* 6892 (527) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 0, 0] */
  10696735,  /* 6893 (527) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  3, 9, 3] */
  1426360 ,  /* 6894 (527) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 6, 0] */
  27821437,  /* 6895 (527) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  9, 9, 9] */
  17015152,  /* 6896 (527) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 0, 0] */
  10698199,  /* 6897 (527) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 9, 3] */
  18552046,  /* 6898 (527) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  27822889,  /* 6899 (527) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 9, 9] */
  128578  ,  /* 6900 (527) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 6] */
  10695487,  /* 6901 (527) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  3, 9, 3] */
  1420432 ,  /* 6902 (527) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 6, 0] */
  26288815,  /* 6903 (527) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 3, 3] */
  134470  ,  /* 6904 (527) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 6] */
  10699915,  /* 6905 (527) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 9, 3] */
  17022544,  /* 6906 (527) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 0] */
  26284495,  /* 6907 (527) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  9, 3, 3] */
  18544690,  /* 6908 (527) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  17140694,  /* 6909 (527) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  26402801,  /* 6910 (527) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  18426488,  /* 6911 (527) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 0] */
  9282479 ,  /* 6912 (527) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 3] */
  17140530,  /* 6913 (527) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  27701391,  /* 6914 (527) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 3] */
  18426648,  /* 6915 (527) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 0] */
  10818237,  /* 6916 (527) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  17135026,  /* 6917 (527) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  26405401,  /* 6918 (527) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 3, 9] */
  10324   ,  /* 6919 (527) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10816609,  /* 6920 (527) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 9] */
  17133562,  /* 6921 (527) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  9279715 ,  /* 6922 (527) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 3, 3] */
  8872    ,  /* 6923 (527) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  27703183,  /* 6924 (527) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 9, 3] */
  17136274,  /* 6925 (527) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  26411329,  /* 6926 (527) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 9] */
  1542946 ,  /* 6927 (527) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  27697291,  /* 6928 (527) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 9, 3] */
  17131846,  /* 6929 (527) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  10809217,  /* 6930 (527) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 9, 9] */
  1547266 ,  /* 6931 (527) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  9287071 ,  /* 6932 (527) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 3] */
  26280003,  /* 6933 (527) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 3, 3] */
  17017896,  /* 6934 (527) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  27828561,  /* 6935 (527) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 9] */
  125994  ,  /* 6936 (527) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 6] */
  26280167,  /* 6937 (527) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 3, 3] */
  18553658,  /* 6938 (527) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  27828401,  /* 6939 (527) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 9] */
  1424588 ,  /* 6940 (527) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 6, 0] */
  26285671,  /* 6941 (527) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  9, 3, 3] */
  17015296,  /* 6942 (527) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 0, 0] */
  9398149 ,  /* 6943 (527) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  3, 3, 9] */
  1426216 ,  /* 6944 (527) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 6, 0] */
  26287135,  /* 6945 (527) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 3, 3] */
  128758  ,  /* 6946 (527) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 6] */
  9399601 ,  /* 6947 (527) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 3, 9] */
  18551866,  /* 6948 (527) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  26284423,  /* 6949 (527) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  9, 3, 3] */
  17009368,  /* 6950 (527) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 0, 0] */
  10699879,  /* 6951 (527) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 9, 3] */
  18557758,  /* 6952 (527) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  26288851,  /* 6953 (527) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 3, 3] */
  1433608 ,  /* 6954 (527) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 0] */
  10695559,  /* 6955 (527) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  3, 9, 3] */
  121402  ,  /* 6956 (527) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 6] */
  18439772,  /* 6957 (527) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  27701879,  /* 6958 (527) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 3] */
  17127410,  /* 6959 (527) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  10817753,  /* 6960 (527) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 9] */
  18439608,  /* 6961 (527) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  26402313,  /* 6962 (527) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  17127570,  /* 6963 (527) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  9282963 ,  /* 6964 (527) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 3] */
  18434104,  /* 6965 (527) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  27704479,  /* 6966 (527) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 9, 3] */
  1545598 ,  /* 6967 (527) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  9281335 ,  /* 6968 (527) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 3, 3] */
  18432640,  /* 6969 (527) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 0] */
  10814989,  /* 6970 (527) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 9, 9] */
  1544146 ,  /* 6971 (527) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  26404105,  /* 6972 (527) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 3, 9] */
  18435352,  /* 6973 (527) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  27710407,  /* 6974 (527) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 3] */
  7672    ,  /* 6975 (527) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  26398213,  /* 6976 (527) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 3, 9] */
  18430924,  /* 6977 (527) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 0] */
  9273943 ,  /* 6978 (527) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 3, 3] */
  11992   ,  /* 6979 (527) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  10822345,  /* 6980 (527) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 9] */
  27815277,  /* 6981 (527) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  9, 9, 9] */
  18553170,  /* 6982 (527) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  26293287,  /* 6983 (527) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 3] */
  1425072 ,  /* 6984 (527) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 6, 0] */
  27815441,  /* 6985 (527) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  9, 9, 9] */
  17018384,  /* 6986 (527) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  26293127,  /* 6987 (527) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 3] */
  125510  ,  /* 6988 (527) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 6] */
  27820945,  /* 6989 (527) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  9, 9, 9] */
  18550570,  /* 6990 (527) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  10697227,  /* 6991 (527) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  3, 9, 3] */
  127138  ,  /* 6992 (527) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 6] */
  27822409,  /* 6993 (527) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 9, 9] */
  1427836 ,  /* 6994 (527) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 0] */
  10698679,  /* 6995 (527) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 9, 3] */
  17016592,  /* 6996 (527) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 0, 0] */
  27819697,  /* 6997 (527) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  9, 9, 9] */
  18544642,  /* 6998 (527) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  9400801 ,  /* 6999 (527) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 3, 9] */
  17022484,  /* 7000 (527) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 0] */
  27824125,  /* 7001 (527) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 9, 9] */
  134530  ,  /* 7002 (527) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 6] */
  9396481 ,  /* 7003 (527) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  3, 3, 9] */
  1420480 ,  /* 7004 (527) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 6, 0] */
  16484   ,  /* 7005 (528) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 7006 (528) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  893633  ,  /* 7007 (528) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 3, 9] */
  897519  ,  /* 7008 (528) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 3, 9] */
  26339054,  /* 7009 (528) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 3, 6] */
  26352336,  /* 7010 (528) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 3, 6] */
  25580003,  /* 7011 (528) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 0, 3] */
  25576113,  /* 7012 (528) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 0, 3] */
  783957  ,  /* 7013 (528) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 3, 3] */
  770999  ,  /* 7014 (528) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 3, 3] */
  125508  ,  /* 7015 (528) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 6] */
  130370  ,  /* 7016 (528) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 6] */
  25689675,  /* 7017 (528) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 0, 9] */
  25702637,  /* 7018 (528) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 0, 9] */
  26230026,  /* 7019 (528) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 0] */
  26225168,  /* 7020 (528) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 0] */
  10816   ,  /* 7021 (528) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 7022 (528) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  25578715,  /* 7023 (528) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 0, 3] */
  25577401,  /* 7024 (528) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  9, 0, 3] */
  19201249,  /* 7025 (528) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 9, 3] */
  19200283,  /* 7026 (528) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 9, 3] */
  10639462,  /* 7027 (528) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 0] */
  10637860,  /* 7028 (528) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  3, 9, 0] */
  8571457 ,  /* 7029 (528) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  3, 0, 3] */
  8572435 ,  /* 7030 (528) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 0, 3] */
  17015134,  /* 7031 (528) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 0, 0] */
  17016772,  /* 7032 (528) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 0, 0] */
  27643804,  /* 7033 (528) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  9, 9, 0] */
  27645742,  /* 7034 (528) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 9, 0] */
  2194015 ,  /* 7035 (528) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 9, 3] */
  2195293 ,  /* 7036 (528) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 9, 3] */
  12064   ,  /* 7037 (528) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 7038 (528) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  10645246,  /* 7039 (528) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 0] */
  10632076,  /* 7040 (528) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 0] */
  10105051,  /* 7041 (528) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  3, 6, 9] */
  10109389,  /* 7042 (528) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 6, 9] */
  889045  ,  /* 7043 (528) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 3, 9] */
  902107  ,  /* 7044 (528) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 3, 9] */
  9219280 ,  /* 7045 (528) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  3, 3, 0] */
  9223690 ,  /* 7046 (528) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 3, 0] */
  1420426 ,  /* 7047 (528) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 6, 0] */
  1433608 ,  /* 7048 (528) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 6, 0] */
  2314903 ,  /* 7049 (528) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 9, 9] */
  2310601 ,  /* 7050 (528) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 9, 9] */
  8696581 ,  /* 7051 (528) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 0, 9] */
  8683507 ,  /* 7052 (528) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 0, 9] */
  1551758 ,  /* 7053 (528) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 6, 6] */
  1538472 ,  /* 7054 (528) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 6, 6] */
  2192711 ,  /* 7055 (528) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 9, 3] */
  2196597 ,  /* 7056 (528) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 9, 3] */
  27638132,  /* 7057 (528) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 9, 0] */
  27651414,  /* 7058 (528) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 9, 0] */
  27115277,  /* 7059 (528) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 6, 9] */
  27111387,  /* 7060 (528) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 6, 9] */
  2319231 ,  /* 7061 (528) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 9, 9] */
  2306273 ,  /* 7062 (528) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 9, 9] */
  1424586 ,  /* 7063 (528) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 6, 0] */
  1429448 ,  /* 7064 (528) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 6, 0] */
  26988753,  /* 7065 (528) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  9, 6, 3] */
  27001715,  /* 7066 (528) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  9, 6, 3] */
  27765300,  /* 7067 (528) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 6] */
  27760442,  /* 7068 (528) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 6] */
  1546090 ,  /* 7069 (528) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 6, 6] */
  1544140 ,  /* 7070 (528) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 6, 6] */
  27113989,  /* 7071 (528) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 6, 9] */
  27112675,  /* 7072 (528) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  9, 6, 9] */
  17902171,  /* 7073 (528) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 3, 9] */
  17901205,  /* 7074 (528) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 3, 9] */
  9340384 ,  /* 7075 (528) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 3, 6] */
  9338782 ,  /* 7076 (528) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  3, 3, 6] */
  10106731,  /* 7077 (528) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  3, 6, 9] */
  10107709,  /* 7078 (528) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  3, 6, 9] */
  18550408,  /* 7079 (528) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18552046,  /* 7080 (528) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  26344726,  /* 7081 (528) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  9, 3, 6] */
  26346664,  /* 7082 (528) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  9, 3, 6] */
  894937  ,  /* 7083 (528) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 3, 9] */
  896215  ,  /* 7084 (528) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 3, 9] */
  1547338 ,  /* 7085 (528) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 6, 6] */
  1542892 ,  /* 7086 (528) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 6, 6] */
  9346168 ,  /* 7087 (528) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 6] */
  9332998 ,  /* 7088 (528) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 6] */
  8569777 ,  /* 7089 (528) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  3, 0, 3] */
  8574115 ,  /* 7090 (528) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  3, 0, 3] */
  2188123 ,  /* 7091 (528) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 9, 3] */
  2201185 ,  /* 7092 (528) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 9, 3] */
  10754554,  /* 7093 (528) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  3, 9, 6] */
  10758964,  /* 7094 (528) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  3, 9, 6] */
  121348  ,  /* 7095 (528) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 6] */
  134530  ,  /* 7096 (528) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 6] */
  779629  ,  /* 7097 (528) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 3, 3] */
  775327  ,  /* 7098 (528) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 3, 3] */
  9995659 ,  /* 7099 (528) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 6, 3] */
  9982585 ,  /* 7100 (528) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 6, 3] */
  17140694,  /* 7101 (528) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 0, 6] */
  17127408,  /* 7102 (528) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 0, 6] */
  17781647,  /* 7103 (528) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 3, 3] */
  17785533,  /* 7104 (528) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 3, 3] */
  9214844 ,  /* 7105 (528) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 3, 0] */
  9228126 ,  /* 7106 (528) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 3, 0] */
  8691989 ,  /* 7107 (528) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 0, 9] */
  8688099 ,  /* 7108 (528) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 0, 9] */
  17908167,  /* 7109 (528) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 3, 9] */
  17895209,  /* 7110 (528) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 3, 9] */
  17013522,  /* 7111 (528) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 0, 0] */
  17018384,  /* 7112 (528) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 0, 0] */
  8565465 ,  /* 7113 (528) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 0, 3] */
  8578427 ,  /* 7114 (528) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 0, 3] */
  9342012 ,  /* 7115 (528) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 6] */
  9337154 ,  /* 7116 (528) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 6] */
  17135026,  /* 7117 (528) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17133076,  /* 7118 (528) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 0, 6] */
  8690701 ,  /* 7119 (528) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 0, 9] */
  8689387 ,  /* 7120 (528) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  3, 0, 9] */
  2313235 ,  /* 7121 (528) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 9, 9] */
  2312269 ,  /* 7122 (528) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 9, 9] */
  27763672,  /* 7123 (528) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 9, 6] */
  27762070,  /* 7124 (528) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  9, 9, 6] */
  25695667,  /* 7125 (528) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  9, 0, 9] */
  25696645,  /* 7126 (528) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 0, 9] */
  127120  ,  /* 7127 (528) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 6] */
  128758  ,  /* 7128 (528) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 6] */
  10755790,  /* 7129 (528) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  3, 9, 6] */
  10757728,  /* 7130 (528) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 9, 6] */
  19318225,  /* 7131 (528) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 9, 9] */
  19319503,  /* 7132 (528) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 9, 9] */
  17136274,  /* 7133 (528) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 0, 6] */
  17131828,  /* 7134 (528) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 0, 6] */
  27769456,  /* 7135 (528) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 6] */
  27756286,  /* 7136 (528) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 6] */
  26993065,  /* 7137 (528) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  9, 6, 3] */
  26997403,  /* 7138 (528) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 6, 3] */
  17777059,  /* 7139 (528) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 3, 3] */
  17790121,  /* 7140 (528) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 3, 3] */
  26343490,  /* 7141 (528) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  9, 3, 6] */
  26347900,  /* 7142 (528) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 3, 6] */
  18544636,  /* 7143 (528) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18557818,  /* 7144 (528) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  19202917,  /* 7145 (528) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 9, 3] */
  19198615,  /* 7146 (528) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 9, 3] */
  25584595,  /* 7147 (528) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 0, 3] */
  25571521,  /* 7148 (528) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 0, 3] */
  18439772,  /* 7149 (528) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  18426486,  /* 7150 (528) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  19316921,  /* 7151 (528) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 9, 9] */
  19320807,  /* 7152 (528) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 9, 9] */
  10750118,  /* 7153 (528) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 9, 6] */
  10763400,  /* 7154 (528) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 9, 6] */
  9991067 ,  /* 7155 (528) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 6, 3] */
  9987177 ,  /* 7156 (528) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 6, 3] */
  19207245,  /* 7157 (528) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 9, 3] */
  19194287,  /* 7158 (528) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 9, 3] */
  18548796,  /* 7159 (528) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 7160 (528) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  10100739,  /* 7161 (528) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  3, 6, 9] */
  10113701,  /* 7162 (528) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  3, 6, 9] */
  10641090,  /* 7163 (528) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 0] */
  10636232,  /* 7164 (528) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 0] */
  18434104,  /* 7165 (528) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  18432154,  /* 7166 (528) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  9989779 ,  /* 7167 (528) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 6, 3] */
  9988465 ,  /* 7168 (528) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  3, 6, 3] */
  777961  ,  /* 7169 (528) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 3, 3] */
  776995  ,  /* 7170 (528) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 3, 3] */
  26228398,  /* 7171 (528) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 0] */
  26226796,  /* 7172 (528) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  9, 3, 0] */
  26994745,  /* 7173 (528) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  9, 6, 3] */
  26995723,  /* 7174 (528) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  9, 6, 3] */
  1426198 ,  /* 7175 (528) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 6, 0] */
  1427836 ,  /* 7176 (528) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 6, 0] */
  9220516 ,  /* 7177 (528) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  3, 3, 0] */
  9222454 ,  /* 7178 (528) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  3, 3, 0] */
  17782951,  /* 7179 (528) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 3, 3] */
  17784229,  /* 7180 (528) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 3, 3] */
  18435352,  /* 7181 (528) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  18430906,  /* 7182 (528) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  26234182,  /* 7183 (528) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 0] */
  26221012,  /* 7184 (528) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 0] */
  25693987,  /* 7185 (528) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  9, 0, 9] */
  25698325,  /* 7186 (528) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  9, 0, 9] */
  19312333,  /* 7187 (528) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 9, 9] */
  19325395,  /* 7188 (528) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 9, 9] */
  27642568,  /* 7189 (528) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  9, 9, 0] */
  27646978,  /* 7190 (528) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  9, 9, 0] */
  17009362,  /* 7191 (528) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 0, 0] */
  17022544,  /* 7192 (528) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 0, 0] */
  17903839,  /* 7193 (528) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 3, 9] */
  17899537,  /* 7194 (528) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 3, 9] */
  27119869,  /* 7195 (528) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 6, 9] */
  27106795,  /* 7196 (528) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 6, 9] */
  16484   ,  /* 7197 (529) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 7198 (529) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  7898    ,  /* 7199 (529) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  11784   ,  /* 7200 (529) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  3200    ,  /* 7201 (529) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  16482   ,  /* 7202 (529) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  11786   ,  /* 7203 (529) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  7896    ,  /* 7204 (529) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  16320   ,  /* 7205 (529) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  3362    ,  /* 7206 (529) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  7410    ,  /* 7207 (529) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  0, 0, 0] */
  12272   ,  /* 7208 (529) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0] */
  3360    ,  /* 7209 (529) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 0, 0] */
  16322   ,  /* 7210 (529) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 0, 0] */
  12270   ,  /* 7211 (529) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  0, 0, 0] */
  7412    ,  /* 7212 (529) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  0, 0, 0] */
  10816   ,  /* 7213 (529) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 7214 (529) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10498   ,  /* 7215 (529) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9184    ,  /* 7216 (529) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  10324   ,  /* 7217 (529) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9358    ,  /* 7218 (529) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10642   ,  /* 7219 (529) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  0, 0, 0] */
  9040    ,  /* 7220 (529) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  0, 0, 0] */
  9352    ,  /* 7221 (529) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  10330   ,  /* 7222 (529) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  9022    ,  /* 7223 (529) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  10660   ,  /* 7224 (529) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  8872    ,  /* 7225 (529) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 0, 0] */
  10810   ,  /* 7226 (529) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9202    ,  /* 7227 (529) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  0, 0, 0] */
  10480   ,  /* 7228 (529) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  0, 0, 0] */
  12064   ,  /* 7229 (529) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 7230 (529) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  16426   ,  /* 7231 (529) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  3256    ,  /* 7232 (529) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  7672    ,  /* 7233 (529) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  12010   ,  /* 7234 (529) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  3310    ,  /* 7235 (529) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  16372   ,  /* 7236 (529) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  7636    ,  /* 7237 (529) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 0, 0] */
  12046   ,  /* 7238 (529) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 0, 0] */
  3250    ,  /* 7239 (529) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  0, 0, 0] */
  16432   ,  /* 7240 (529) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0] */
  11992   ,  /* 7241 (529) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  7690    ,  /* 7242 (529) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  16378   ,  /* 7243 (529) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  0, 0, 0] */
  3304    ,  /* 7244 (529) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  0, 0, 0] */
  18557870,  /* 7245 (529) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 7246 (529) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18549284,  /* 7247 (529) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18553170,  /* 7248 (529) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18544586,  /* 7249 (529) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18557868,  /* 7250 (529) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18553172,  /* 7251 (529) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18549282,  /* 7252 (529) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18557706,  /* 7253 (529) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  18544748,  /* 7254 (529) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18548796,  /* 7255 (529) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18553658,  /* 7256 (529) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18544746,  /* 7257 (529) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 6, 6] */
  18557708,  /* 7258 (529) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 6, 6] */
  18553656,  /* 7259 (529) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  6, 6, 6] */
  18548798,  /* 7260 (529) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  6, 6, 6] */
  18552202,  /* 7261 (529) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550252,  /* 7262 (529) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18551884,  /* 7263 (529) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18550570,  /* 7264 (529) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18551710,  /* 7265 (529) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550744,  /* 7266 (529) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18552028,  /* 7267 (529) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  6, 6, 6] */
  18550426,  /* 7268 (529) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  6, 6, 6] */
  18550738,  /* 7269 (529) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18551716,  /* 7270 (529) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550408,  /* 7271 (529) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18552046,  /* 7272 (529) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18550258,  /* 7273 (529) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18552196,  /* 7274 (529) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 6, 6] */
  18550588,  /* 7275 (529) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  6, 6, 6] */
  18551866,  /* 7276 (529) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  6, 6, 6] */
  18553450,  /* 7277 (529) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18549004,  /* 7278 (529) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18557812,  /* 7279 (529) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  18544642,  /* 7280 (529) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18549058,  /* 7281 (529) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18553396,  /* 7282 (529) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18544696,  /* 7283 (529) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18557758,  /* 7284 (529) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18549022,  /* 7285 (529) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 6, 6] */
  18553432,  /* 7286 (529) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 6, 6] */
  18544636,  /* 7287 (529) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  6, 6, 6] */
  18557818,  /* 7288 (529) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  6, 6, 6] */
  18553378,  /* 7289 (529) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  18549076,  /* 7290 (529) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18557764,  /* 7291 (529) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  6, 6, 6] */
  18544690,  /* 7292 (529) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  6, 6, 6] */
  16484   ,  /* 7293 (530) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0] */
  3198    ,  /* 7294 (530) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 0] */
  10695767,  /* 7295 (530) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  3, 9, 3] */
  10699653,  /* 7296 (530) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  3, 9, 3] */
  17127410,  /* 7297 (530) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 6] */
  17140692,  /* 7298 (530) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 6] */
  9400577 ,  /* 7299 (530) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  3, 3, 9] */
  9396687 ,  /* 7300 (530) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  3, 3, 9] */
  134418  ,  /* 7301 (530) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  0, 0, 6] */
  121460  ,  /* 7302 (530) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 6] */
  9278103 ,  /* 7303 (530) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  3, 3, 3] */
  9282965 ,  /* 7304 (530) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  3, 3, 3] */
  17009472,  /* 7305 (530) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  6, 0, 0] */
  17022434,  /* 7306 (530) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  6, 0, 0] */
  10818237,  /* 7307 (530) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  3, 9, 9] */
  10813379,  /* 7308 (530) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  3, 9, 9] */
  10816   ,  /* 7309 (530) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  0, 0, 0] */
  8866    ,  /* 7310 (530) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  0, 0, 0] */
  9399289 ,  /* 7311 (530) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  3, 3, 9] */
  9397975 ,  /* 7312 (530) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  3, 3, 9] */
  18433612,  /* 7313 (530) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  6, 6, 0] */
  18432646,  /* 7314 (530) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  6, 6, 0] */
  26287447,  /* 7315 (530) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  9, 3, 3] */
  26285845,  /* 7316 (530) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  9, 3, 3] */
  17015464,  /* 7317 (530) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  6, 0, 0] */
  17016442,  /* 7318 (530) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  6, 0, 0] */
  9279715 ,  /* 7319 (530) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  3, 3, 3] */
  9281353 ,  /* 7320 (530) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  3, 3, 3] */
  1426048 ,  /* 7321 (530) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  0, 6, 0] */
  1427986 ,  /* 7322 (530) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  0, 6, 0] */
  26404105,  /* 7323 (530) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  9, 3, 9] */
  26405383,  /* 7324 (530) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  9, 3, 9] */
  12064   ,  /* 7325 (530) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 0] */
  7618    ,  /* 7326 (530) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 0] */
  26293231,  /* 7327 (530) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 3, 3] */
  26280061,  /* 7328 (530) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 3, 3] */
  1542946 ,  /* 7329 (530) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 6] */
  1547284 ,  /* 7330 (530) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 6] */
  10691179,  /* 7331 (530) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 9, 3] */
  10704241,  /* 7332 (530) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 9, 3] */
  1424812 ,  /* 7333 (530) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  0, 6, 0] */
  1429222 ,  /* 7334 (530) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  0, 6, 0] */
  9273943 ,  /* 7335 (530) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  3, 3, 3] */
  9287125 ,  /* 7336 (530) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  3, 3, 3] */
  130090  ,  /* 7337 (530) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  0, 0, 6] */
  125788  ,  /* 7338 (530) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  0, 0, 6] */
  27710359,  /* 7339 (530) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  9, 9, 3] */
  27697285,  /* 7340 (530) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  9, 9, 3] */
  18557870,  /* 7341 (530) [  1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 6] */
  18544584,  /* 7342 (530) [ -1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 6] */
  26402801,  /* 7343 (530) [  0,-1, 0, 1, 0, 0, 0, 0, 1,  9, 3, 9] */
  26406687,  /* 7344 (530) [  0, 1, 0,-1, 0, 0, 0, 0,-1,  9, 3, 9] */
  1420376 ,  /* 7345 (530) [ -1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 0] */
  1433658 ,  /* 7346 (530) [  1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 0] */
  27705767,  /* 7347 (530) [  0, 1, 0,-1, 0, 0, 0, 0, 1,  9, 9, 3] */
  27701877,  /* 7348 (530) [  0,-1, 0, 1, 0, 0, 0, 0,-1,  9, 9, 3] */
  18439608,  /* 7349 (530) [  1, 0, 0, 0,-1, 0, 0, 0,-1,  6, 6, 0] */
  18426650,  /* 7350 (530) [ -1, 0, 0, 0, 1, 0, 0, 0, 1,  6, 6, 0] */
  27819489,  /* 7351 (530) [  0,-1, 0,-1, 0, 0, 0, 0,-1,  9, 9, 9] */
  27824351,  /* 7352 (530) [  0, 1, 0, 1, 0, 0, 0, 0, 1,  9, 9, 9] */
  1538634 ,  /* 7353 (530) [ -1, 0, 0, 0, 1, 0, 0, 0,-1,  0, 6, 6] */
  1551596 ,  /* 7354 (530) [  1, 0, 0, 0,-1, 0, 0, 0, 1,  0, 6, 6] */
  26289075,  /* 7355 (530) [  0, 1, 0, 1, 0, 0, 0, 0,-1,  9, 3, 3] */
  26284217,  /* 7356 (530) [  0,-1, 0,-1, 0, 0, 0, 0, 1,  9, 3, 3] */
  18552202,  /* 7357 (530) [  0, 0, 1, 1, 0, 0, 0, 1, 0,  6, 6, 6] */
  18550252,  /* 7358 (530) [  0, 0,-1,-1, 0, 0, 0,-1, 0,  6, 6, 6] */
  27704479,  /* 7359 (530) [  0, 0, 1, 0,-1, 0, 1, 0, 0,  9, 9, 3] */
  27703165,  /* 7360 (530) [  0, 0,-1, 0, 1, 0,-1, 0, 0,  9, 9, 3] */
  128422  ,  /* 7361 (530) [  0, 0, 1,-1, 0, 0, 0,-1, 0,  0, 0, 6] */
  127456  ,  /* 7362 (530) [  0, 0,-1, 1, 0, 0, 0, 1, 0,  0, 0, 6] */
  10816609,  /* 7363 (530) [  0, 0, 1, 0, 1, 0,-1, 0, 0,  3, 9, 9] */
  10815007,  /* 7364 (530) [  0, 0,-1, 0,-1, 0, 1, 0, 0,  3, 9, 9] */
  1544626 ,  /* 7365 (530) [  0, 0,-1, 1, 0, 0, 0,-1, 0,  0, 6, 6] */
  1545604 ,  /* 7366 (530) [  0, 0, 1,-1, 0, 0, 0, 1, 0,  0, 6, 6] */
  27821101,  /* 7367 (530) [  0, 0,-1, 0,-1, 0,-1, 0, 0,  9, 9, 9] */
  27822739,  /* 7368 (530) [  0, 0, 1, 0, 1, 0, 1, 0, 0,  9, 9, 9] */
  17133082,  /* 7369 (530) [  0, 0,-1,-1, 0, 0, 0, 1, 0,  6, 0, 6] */
  17135020,  /* 7370 (530) [  0, 0, 1, 1, 0, 0, 0,-1, 0,  6, 0, 6] */
  10697071,  /* 7371 (530) [  0, 0,-1, 0, 1, 0, 1, 0, 0,  3, 9, 3] */
  10698349,  /* 7372 (530) [  0, 0, 1, 0,-1, 0,-1, 0, 0,  3, 9, 3] */
  18553450,  /* 7373 (530) [  0, 1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 6] */
  18549004,  /* 7374 (530) [  0,-1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 6] */
  10822393,  /* 7375 (530) [  1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 9, 9] */
  10809223,  /* 7376 (530) [ -1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 9, 9] */
  17013784,  /* 7377 (530) [  0,-1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 0] */
  17018122,  /* 7378 (530) [  0, 1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 0] */
  26398213,  /* 7379 (530) [ -1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 3, 9] */
  26411275,  /* 7380 (530) [  1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 3, 9] */
  17131846,  /* 7381 (530) [  0,-1, 0, 0, 0,-1, 1, 0, 0,  6, 0, 6] */
  17136256,  /* 7382 (530) [  0, 1, 0, 0, 0, 1,-1, 0, 0,  6, 0, 6] */
  27815329,  /* 7383 (530) [ -1, 0, 0, 0, 0,-1, 0,-1, 0,  9, 9, 9] */
  27828511,  /* 7384 (530) [  1, 0, 0, 0, 0, 1, 0, 1, 0,  9, 9, 9] */
  18435280,  /* 7385 (530) [  0, 1, 0, 0, 0,-1,-1, 0, 0,  6, 6, 0] */
  18430978,  /* 7386 (530) [  0,-1, 0, 0, 0, 1, 1, 0, 0,  6, 6, 0] */
  9405169 ,  /* 7387 (530) [  1, 0, 0, 0, 0,-1, 0, 1, 0,  3, 3, 9] */
  9392095 ,  /* 7388 (530) [ -1, 0, 0, 0, 0, 1, 0,-1, 0,  3, 3, 9] */
};

static const int symmetry_operation_index[][2] = {
  {   0,    0}, /*   0 */
  {   1,    1}, /*   1 */
  {   2,    2}, /*   2 */
  {   2,    4}, /*   3 */
  {   2,    6}, /*   4 */
  {   2,    8}, /*   5 */
  {   2,   10}, /*   6 */
  {   2,   12}, /*   7 */
  {   2,   14}, /*   8 */
  {   4,   16}, /*   9 */
  {   4,   20}, /*  10 */
  {   4,   24}, /*  11 */
  {   4,   28}, /*  12 */
  {   4,   32}, /*  13 */
  {   4,   36}, /*  14 */
  {   4,   40}, /*  15 */
  {   4,   44}, /*  16 */
  {   4,   48}, /*  17 */
  {   2,   52}, /*  18 */
  {   2,   54}, /*  19 */
  {   2,   56}, /*  20 */
  {   2,   58}, /*  21 */
  {   2,   60}, /*  22 */
  {   2,   62}, /*  23 */
  {   2,   64}, /*  24 */
  {   2,   66}, /*  25 */
  {   2,   68}, /*  26 */
  {   2,   70}, /*  27 */
  {   2,   72}, /*  28 */
  {   2,   74}, /*  29 */
  {   4,   76}, /*  30 */
  {   4,   80}, /*  31 */
  {   4,   84}, /*  32 */
  {   4,   88}, /*  33 */
  {   4,   92}, /*  34 */
  {   4,   96}, /*  35 */
  {   4,  100}, /*  36 */
  {   4,  104}, /*  37 */
  {   4,  108}, /*  38 */
  {   4,  112}, /*  39 */
  {   4,  116}, /*  40 */
  {   4,  120}, /*  41 */
  {   4,  124}, /*  42 */
  {   4,  128}, /*  43 */
  {   4,  132}, /*  44 */
  {   4,  136}, /*  45 */
  {   4,  140}, /*  46 */
  {   4,  144}, /*  47 */
  {   4,  148}, /*  48 */
  {   4,  152}, /*  49 */
  {   4,  156}, /*  50 */
  {   4,  160}, /*  51 */
  {   4,  164}, /*  52 */
  {   4,  168}, /*  53 */
  {   4,  172}, /*  54 */
  {   4,  176}, /*  55 */
  {   4,  180}, /*  56 */
  {   4,  184}, /*  57 */
  {   4,  188}, /*  58 */
  {   4,  192}, /*  59 */
  {   4,  196}, /*  60 */
  {   4,  200}, /*  61 */
  {   4,  204}, /*  62 */
  {   8,  208}, /*  63 */
  {   8,  216}, /*  64 */
  {   8,  224}, /*  65 */
  {   8,  232}, /*  66 */
  {   8,  240}, /*  67 */
  {   8,  248}, /*  68 */
  {   8,  256}, /*  69 */
  {   8,  264}, /*  70 */
  {   8,  272}, /*  71 */
  {   4,  280}, /*  72 */
  {   4,  284}, /*  73 */
  {   4,  288}, /*  74 */
  {   4,  292}, /*  75 */
  {   4,  296}, /*  76 */
  {   4,  300}, /*  77 */
  {   4,  304}, /*  78 */
  {   4,  308}, /*  79 */
  {   4,  312}, /*  80 */
  {   4,  316}, /*  81 */
  {   4,  320}, /*  82 */
  {   4,  324}, /*  83 */
  {   4,  328}, /*  84 */
  {   4,  332}, /*  85 */
  {   4,  336}, /*  86 */
  {   4,  340}, /*  87 */
  {   4,  344}, /*  88 */
  {   4,  348}, /*  89 */
  {   8,  352}, /*  90 */
  {   8,  360}, /*  91 */
  {   8,  368}, /*  92 */
  {   8,  376}, /*  93 */
  {   8,  384}, /*  94 */
  {   8,  392}, /*  95 */
  {   8,  400}, /*  96 */
  {   8,  408}, /*  97 */
  {   8,  416}, /*  98 */
  {   8,  424}, /*  99 */
  {   8,  432}, /* 100 */
  {   8,  440}, /* 101 */
  {   8,  448}, /* 102 */
  {   8,  456}, /* 103 */
  {   8,  464}, /* 104 */
  {   8,  472}, /* 105 */
  {   8,  480}, /* 106 */
  {   8,  488}, /* 107 */
  {   4,  496}, /* 108 */
  {   4,  500}, /* 109 */
  {   4,  504}, /* 110 */
  {   4,  508}, /* 111 */
  {   4,  512}, /* 112 */
  {   4,  516}, /* 113 */
  {   4,  520}, /* 114 */
  {   4,  524}, /* 115 */
  {   8,  528}, /* 116 */
  {   8,  536}, /* 117 */
  {   8,  544}, /* 118 */
  {   8,  552}, /* 119 */
  {   8,  560}, /* 120 */
  {   8,  568}, /* 121 */
  {  16,  576}, /* 122 */
  {   8,  592}, /* 123 */
  {   8,  600}, /* 124 */
  {   4,  608}, /* 125 */
  {   4,  612}, /* 126 */
  {   4,  616}, /* 127 */
  {   4,  620}, /* 128 */
  {   4,  624}, /* 129 */
  {   4,  628}, /* 130 */
  {   4,  632}, /* 131 */
  {   4,  636}, /* 132 */
  {   4,  640}, /* 133 */
  {   4,  644}, /* 134 */
  {   4,  648}, /* 135 */
  {   4,  652}, /* 136 */
  {   4,  656}, /* 137 */
  {   4,  660}, /* 138 */
  {   4,  664}, /* 139 */
  {   4,  668}, /* 140 */
  {   4,  672}, /* 141 */
  {   4,  676}, /* 142 */
  {   4,  680}, /* 143 */
  {   4,  684}, /* 144 */
  {   4,  688}, /* 145 */
  {   4,  692}, /* 146 */
  {   4,  696}, /* 147 */
  {   4,  700}, /* 148 */
  {   4,  704}, /* 149 */
  {   4,  708}, /* 150 */
  {   4,  712}, /* 151 */
  {   4,  716}, /* 152 */
  {   4,  720}, /* 153 */
  {   4,  724}, /* 154 */
  {   4,  728}, /* 155 */
  {   4,  732}, /* 156 */
  {   4,  736}, /* 157 */
  {   4,  740}, /* 158 */
  {   4,  744}, /* 159 */
  {   4,  748}, /* 160 */
  {   4,  752}, /* 161 */
  {   4,  756}, /* 162 */
  {   4,  760}, /* 163 */
  {   4,  764}, /* 164 */
  {   4,  768}, /* 165 */
  {   4,  772}, /* 166 */
  {   4,  776}, /* 167 */
  {   4,  780}, /* 168 */
  {   4,  784}, /* 169 */
  {   4,  788}, /* 170 */
  {   4,  792}, /* 171 */
  {   4,  796}, /* 172 */
  {   8,  800}, /* 173 */
  {   8,  808}, /* 174 */
  {   8,  816}, /* 175 */
  {   8,  824}, /* 176 */
  {   8,  832}, /* 177 */
  {   8,  840}, /* 178 */
  {   8,  848}, /* 179 */
  {   8,  856}, /* 180 */
  {   8,  864}, /* 181 */
  {   8,  872}, /* 182 */
  {   8,  880}, /* 183 */
  {   8,  888}, /* 184 */
  {   8,  896}, /* 185 */
  {   8,  904}, /* 186 */
  {   8,  912}, /* 187 */
  {   8,  920}, /* 188 */
  {   8,  928}, /* 189 */
  {   8,  936}, /* 190 */
  {   8,  944}, /* 191 */
  {   8,  952}, /* 192 */
  {   8,  960}, /* 193 */
  {   8,  968}, /* 194 */
  {   8,  976}, /* 195 */
  {   8,  984}, /* 196 */
  {   8,  992}, /* 197 */
  {   8, 1000}, /* 198 */
  {   8, 1008}, /* 199 */
  {   8, 1016}, /* 200 */
  {   8, 1024}, /* 201 */
  {   8, 1032}, /* 202 */
  {   8, 1040}, /* 203 */
  {   8, 1048}, /* 204 */
  {   8, 1056}, /* 205 */
  {   8, 1064}, /* 206 */
  {   8, 1072}, /* 207 */
  {   8, 1080}, /* 208 */
  {  16, 1088}, /* 209 */
  {  16, 1104}, /* 210 */
  {  16, 1120}, /* 211 */
  {  16, 1136}, /* 212 */
  {  16, 1152}, /* 213 */
  {  16, 1168}, /* 214 */
  {   8, 1184}, /* 215 */
  {   8, 1192}, /* 216 */
  {   8, 1200}, /* 217 */
  {   8, 1208}, /* 218 */
  {   8, 1216}, /* 219 */
  {   8, 1224}, /* 220 */
  {   8, 1232}, /* 221 */
  {   8, 1240}, /* 222 */
  {   8, 1248}, /* 223 */
  {   8, 1256}, /* 224 */
  {   8, 1264}, /* 225 */
  {   8, 1272}, /* 226 */
  {   8, 1280}, /* 227 */
  {   8, 1288}, /* 228 */
  {   8, 1296}, /* 229 */
  {   8, 1304}, /* 230 */
  {   8, 1312}, /* 231 */
  {   8, 1320}, /* 232 */
  {   8, 1328}, /* 233 */
  {   8, 1336}, /* 234 */
  {   8, 1344}, /* 235 */
  {   8, 1352}, /* 236 */
  {   8, 1360}, /* 237 */
  {   8, 1368}, /* 238 */
  {   8, 1376}, /* 239 */
  {   8, 1384}, /* 240 */
  {   8, 1392}, /* 241 */
  {   8, 1400}, /* 242 */
  {   8, 1408}, /* 243 */
  {   8, 1416}, /* 244 */
  {   8, 1424}, /* 245 */
  {   8, 1432}, /* 246 */
  {   8, 1440}, /* 247 */
  {   8, 1448}, /* 248 */
  {   8, 1456}, /* 249 */
  {   8, 1464}, /* 250 */
  {   8, 1472}, /* 251 */
  {   8, 1480}, /* 252 */
  {   8, 1488}, /* 253 */
  {   8, 1496}, /* 254 */
  {   8, 1504}, /* 255 */
  {   8, 1512}, /* 256 */
  {   8, 1520}, /* 257 */
  {   8, 1528}, /* 258 */
  {   8, 1536}, /* 259 */
  {   8, 1544}, /* 260 */
  {   8, 1552}, /* 261 */
  {   8, 1560}, /* 262 */
  {   8, 1568}, /* 263 */
  {   8, 1576}, /* 264 */
  {   8, 1584}, /* 265 */
  {   8, 1592}, /* 266 */
  {   8, 1600}, /* 267 */
  {   8, 1608}, /* 268 */
  {   8, 1616}, /* 269 */
  {   8, 1624}, /* 270 */
  {   8, 1632}, /* 271 */
  {   8, 1640}, /* 272 */
  {   8, 1648}, /* 273 */
  {   8, 1656}, /* 274 */
  {   8, 1664}, /* 275 */
  {   8, 1672}, /* 276 */
  {   8, 1680}, /* 277 */
  {   8, 1688}, /* 278 */
  {   8, 1696}, /* 279 */
  {   8, 1704}, /* 280 */
  {   8, 1712}, /* 281 */
  {   8, 1720}, /* 282 */
  {   8, 1728}, /* 283 */
  {   8, 1736}, /* 284 */
  {   8, 1744}, /* 285 */
  {   8, 1752}, /* 286 */
  {   8, 1760}, /* 287 */
  {   8, 1768}, /* 288 */
  {   8, 1776}, /* 289 */
  {   8, 1784}, /* 290 */
  {   8, 1792}, /* 291 */
  {   8, 1800}, /* 292 */
  {   8, 1808}, /* 293 */
  {   8, 1816}, /* 294 */
  {   8, 1824}, /* 295 */
  {   8, 1832}, /* 296 */
  {   8, 1840}, /* 297 */
  {  16, 1848}, /* 298 */
  {  16, 1864}, /* 299 */
  {  16, 1880}, /* 300 */
  {  16, 1896}, /* 301 */
  {  16, 1912}, /* 302 */
  {  16, 1928}, /* 303 */
  {  16, 1944}, /* 304 */
  {  16, 1960}, /* 305 */
  {  16, 1976}, /* 306 */
  {  16, 1992}, /* 307 */
  {  16, 2008}, /* 308 */
  {  16, 2024}, /* 309 */
  {  16, 2040}, /* 310 */
  {  16, 2056}, /* 311 */
  {  16, 2072}, /* 312 */
  {  16, 2088}, /* 313 */
  {  16, 2104}, /* 314 */
  {  16, 2120}, /* 315 */
  {  16, 2136}, /* 316 */
  {  16, 2152}, /* 317 */
  {  16, 2168}, /* 318 */
  {  16, 2184}, /* 319 */
  {  16, 2200}, /* 320 */
  {  16, 2216}, /* 321 */
  {  16, 2232}, /* 322 */
  {  16, 2248}, /* 323 */
  {  16, 2264}, /* 324 */
  {  16, 2280}, /* 325 */
  {  16, 2296}, /* 326 */
  {  16, 2312}, /* 327 */
  {  16, 2328}, /* 328 */
  {  16, 2344}, /* 329 */
  {  16, 2360}, /* 330 */
  {  16, 2376}, /* 331 */
  {  16, 2392}, /* 332 */
  {  16, 2408}, /* 333 */
  {  32, 2424}, /* 334 */
  {  32, 2456}, /* 335 */
  {  32, 2488}, /* 336 */
  {  16, 2520}, /* 337 */
  {  16, 2536}, /* 338 */
  {  16, 2552}, /* 339 */
  {  16, 2568}, /* 340 */
  {  16, 2584}, /* 341 */
  {  16, 2600}, /* 342 */
  {  16, 2616}, /* 343 */
  {  16, 2632}, /* 344 */
  {  16, 2648}, /* 345 */
  {  16, 2664}, /* 346 */
  {  16, 2680}, /* 347 */
  {  16, 2696}, /* 348 */
  {   4, 2712}, /* 349 */
  {   4, 2716}, /* 350 */
  {   4, 2720}, /* 351 */
  {   4, 2724}, /* 352 */
  {   8, 2728}, /* 353 */
  {   8, 2736}, /* 354 */
  {   4, 2744}, /* 355 */
  {   8, 2748}, /* 356 */
  {   8, 2756}, /* 357 */
  {   8, 2764}, /* 358 */
  {   8, 2772}, /* 359 */
  {   8, 2780}, /* 360 */
  {   8, 2788}, /* 361 */
  {   8, 2796}, /* 362 */
  {  16, 2804}, /* 363 */
  {  16, 2820}, /* 364 */
  {  16, 2836}, /* 365 */
  {   8, 2852}, /* 366 */
  {   8, 2860}, /* 367 */
  {   8, 2868}, /* 368 */
  {   8, 2876}, /* 369 */
  {   8, 2884}, /* 370 */
  {   8, 2892}, /* 371 */
  {   8, 2900}, /* 372 */
  {   8, 2908}, /* 373 */
  {  16, 2916}, /* 374 */
  {  16, 2932}, /* 375 */
  {   8, 2948}, /* 376 */
  {   8, 2956}, /* 377 */
  {   8, 2964}, /* 378 */
  {   8, 2972}, /* 379 */
  {   8, 2980}, /* 380 */
  {   8, 2988}, /* 381 */
  {   8, 2996}, /* 382 */
  {   8, 3004}, /* 383 */
  {  16, 3012}, /* 384 */
  {  16, 3028}, /* 385 */
  {  16, 3044}, /* 386 */
  {  16, 3060}, /* 387 */
  {   8, 3076}, /* 388 */
  {   8, 3084}, /* 389 */
  {   8, 3092}, /* 390 */
  {   8, 3100}, /* 391 */
  {   8, 3108}, /* 392 */
  {   8, 3116}, /* 393 */
  {   8, 3124}, /* 394 */
  {   8, 3132}, /* 395 */
  {  16, 3140}, /* 396 */
  {  16, 3156}, /* 397 */
  {  16, 3172}, /* 398 */
  {  16, 3188}, /* 399 */
  {  16, 3204}, /* 400 */
  {  16, 3220}, /* 401 */
  {  16, 3236}, /* 402 */
  {  16, 3252}, /* 403 */
  {  16, 3268}, /* 404 */
  {  16, 3284}, /* 405 */
  {  16, 3300}, /* 406 */
  {  16, 3316}, /* 407 */
  {  16, 3332}, /* 408 */
  {  16, 3348}, /* 409 */
  {  16, 3364}, /* 410 */
  {  16, 3380}, /* 411 */
  {  16, 3396}, /* 412 */
  {  16, 3412}, /* 413 */
  {  16, 3428}, /* 414 */
  {  16, 3444}, /* 415 */
  {  16, 3460}, /* 416 */
  {  16, 3476}, /* 417 */
  {  16, 3492}, /* 418 */
  {  16, 3508}, /* 419 */
  {  16, 3524}, /* 420 */
  {  16, 3540}, /* 421 */
  {  16, 3556}, /* 422 */
  {  16, 3572}, /* 423 */
  {  32, 3588}, /* 424 */
  {  32, 3620}, /* 425 */
  {  32, 3652}, /* 426 */
  {  32, 3684}, /* 427 */
  {  32, 3716}, /* 428 */
  {  32, 3748}, /* 429 */
  {   3, 3780}, /* 430 */
  {   3, 3783}, /* 431 */
  {   3, 3786}, /* 432 */
  {   9, 3789}, /* 433 */
  {   3, 3798}, /* 434 */
  {   6, 3801}, /* 435 */
  {  18, 3807}, /* 436 */
  {   6, 3825}, /* 437 */
  {   6, 3831}, /* 438 */
  {   6, 3837}, /* 439 */
  {   6, 3843}, /* 440 */
  {   6, 3849}, /* 441 */
  {   6, 3855}, /* 442 */
  {   6, 3861}, /* 443 */
  {  18, 3867}, /* 444 */
  {   6, 3885}, /* 445 */
  {   6, 3891}, /* 446 */
  {   6, 3897}, /* 447 */
  {   6, 3903}, /* 448 */
  {   6, 3909}, /* 449 */
  {  18, 3915}, /* 450 */
  {   6, 3933}, /* 451 */
  {  18, 3939}, /* 452 */
  {   6, 3957}, /* 453 */
  {  12, 3963}, /* 454 */
  {  12, 3975}, /* 455 */
  {  12, 3987}, /* 456 */
  {  12, 3999}, /* 457 */
  {  36, 4011}, /* 458 */
  {  12, 4047}, /* 459 */
  {  36, 4059}, /* 460 */
  {  12, 4095}, /* 461 */
  {   6, 4107}, /* 462 */
  {   6, 4113}, /* 463 */
  {   6, 4119}, /* 464 */
  {   6, 4125}, /* 465 */
  {   6, 4131}, /* 466 */
  {   6, 4137}, /* 467 */
  {   6, 4143}, /* 468 */
  {  12, 4149}, /* 469 */
  {  12, 4161}, /* 470 */
  {  12, 4173}, /* 471 */
  {  12, 4185}, /* 472 */
  {  12, 4197}, /* 473 */
  {  12, 4209}, /* 474 */
  {  12, 4221}, /* 475 */
  {  12, 4233}, /* 476 */
  {  12, 4245}, /* 477 */
  {  12, 4257}, /* 478 */
  {  12, 4269}, /* 479 */
  {  12, 4281}, /* 480 */
  {  12, 4293}, /* 481 */
  {  12, 4305}, /* 482 */
  {  12, 4317}, /* 483 */
  {  12, 4329}, /* 484 */
  {  24, 4341}, /* 485 */
  {  24, 4365}, /* 486 */
  {  24, 4389}, /* 487 */
  {  24, 4413}, /* 488 */
  {  12, 4437}, /* 489 */
  {  48, 4449}, /* 490 */
  {  24, 4497}, /* 491 */
  {  12, 4521}, /* 492 */
  {  24, 4533}, /* 493 */
  {  24, 4557}, /* 494 */
  {  24, 4581}, /* 495 */
  {  24, 4605}, /* 496 */
  {  96, 4629}, /* 497 */
  {  96, 4725}, /* 498 */
  {  96, 4821}, /* 499 */
  {  48, 4917}, /* 500 */
  {  24, 4965}, /* 501 */
  {  48, 4989}, /* 502 */
  {  24, 5037}, /* 503 */
  {  24, 5061}, /* 504 */
  {  96, 5085}, /* 505 */
  {  96, 5181}, /* 506 */
  {  48, 5277}, /* 507 */
  {  24, 5325}, /* 508 */
  {  24, 5349}, /* 509 */
  {  48, 5373}, /* 510 */
  {  24, 5421}, /* 511 */
  {  96, 5445}, /* 512 */
  {  48, 5541}, /* 513 */
  {  24, 5589}, /* 514 */
  {  96, 5613}, /* 515 */
  {  48, 5709}, /* 516 */
  {  48, 5757}, /* 517 */
  {  48, 5805}, /* 518 */
  {  48, 5853}, /* 519 */
  {  48, 5901}, /* 520 */
  {  48, 5949}, /* 521 */
  {  48, 5997}, /* 522 */
  { 192, 6045}, /* 523 */
  { 192, 6237}, /* 524 */
  { 192, 6429}, /* 525 */
  { 192, 6621}, /* 526 */
  { 192, 6813}, /* 527 */
  { 192, 7005}, /* 528 */
  {  96, 7197}, /* 529 */
  {  96, 7293}, /* 530 */
};

static int remove_space(char symbol[], const int num_char);
static void replace_equal_char(char symbol[], const int position);

int spgdb_get_operation(int rot[3][3], double trans[3], const int hall_number)
{
  int i, j, r, t, degit;

  /* A space group operation is compressed using ternary numerical system for */
  /* rotation and duodecimal system for translation. This is achieved because */
  /* each element of rotation matrix can have only one of {-1,0,1}, and */
  /* the translation can have one of {0,2,3,4,6,8,9,10} divided by */
  /* 12. Therefore 3^9 * 12^3 = 34012224 different values can map space */
  /* group operations. In principle, octal numerical system can be used */
  /* for translation, but duodecimal system is more convenient. */

  r = symmetry_operations[hall_number] % 19683; /* 19683 = 3**9 */
  degit = 6561; /* 6561 = 3**8 */
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      rot[i][j] = ( r % ( degit * 3 ) ) / degit - 1;
      degit /= 3;
    }
  }

  t = symmetry_operations[hall_number] / 19683; /* 19683 = 3**9 */
  degit = 144;
  for ( i = 0; i < 3; i++ ) {
    trans[i] = ( (double) ( ( t % ( degit * 12 ) ) / degit ) ) / 12;
    degit /= 12;
  }

  return 1;
}

void spgdb_get_operation_index(int indices[2], const int hall_number)
{
  indices[0] = symmetry_operation_index[hall_number][0];
  indices[1] = symmetry_operation_index[hall_number][1];
}

/* Return NULL if failed */
Symmetry * spgdb_get_spacegroup_operations(const int hall_number)
{
  int i;
  int operation_index[2];
  int rot[3][3];
  double trans[3];
  Symmetry *symmetry;

  symmetry = NULL;

  if (hall_number < 1 || 530 < hall_number) {
    return NULL;
  }

  spgdb_get_operation_index(operation_index, hall_number);

  if ((symmetry = sym_alloc_symmetry(operation_index[0])) == NULL) {
    return NULL;
  }

  for (i = 0; i < operation_index[0]; i++) {
    /* rotation matrix matching and set difference of translations */
    spgdb_get_operation(rot, trans, operation_index[1] + i);
    mat_copy_matrix_i3(symmetry->rot[i], rot);
    mat_copy_vector_d3(symmetry->trans[i], trans);
  }

  return symmetry;
}

/* Return spgtype.number = 0 if hall_number is out of range. */
SpacegroupType spgdb_get_spacegroup_type(const int hall_number)
{
  int position; 
  SpacegroupType spgtype;

  spgtype.number = 0;
  
  if (0 < hall_number || hall_number < 531) {
    spgtype = spacegroup_types[hall_number];
  } else {
    spgtype = spacegroup_types[0];
  }

  remove_space(spgtype.schoenflies, 7);
  position = remove_space(spgtype.hall_symbol, 17);
  replace_equal_char(spgtype.hall_symbol, position);
  remove_space(spgtype.international, 32);
  remove_space(spgtype.international_full, 20);
  remove_space(spgtype.international_short, 11);
  
  return spgtype;
}

static int remove_space(char symbol[], const int num_char) {
  int i;

  for (i = num_char - 2; i > -1; i--) {
    if (symbol[i] == ' ') {
      symbol[i] = '\0';
    } else {
      return i;
    }
  }
  return i;
}

static void replace_equal_char(char symbol[], const int position) {
  int i;

  for (i = position; i > -1; i--) {
    if (symbol[i] == '=') { symbol[i] = '\"'; }
  }
}
 
