#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
#AC_INIT([spglib], 1.9.1, [atz.togo@gmail.com])
AC_INIT([spglib],
        m4_esyscmd_s(grep VERSION src/version.h|awk '{printf("%d.", $3)}'|sed s/.$//),
        [atz.togo@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/spglib.c])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
# AC_PROG_RANLIB
AM_PROG_CC_C_O
AM_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_FUNCS([sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
