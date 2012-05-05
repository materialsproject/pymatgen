#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([spglib], 1.2, [atz.togo@gmail.com])
AM_INIT_AUTOMAKE([spglib], 1.2)
AC_CONFIG_SRCDIR([src/spglib.c])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
# AC_PROG_RANLIB
AM_PROG_CC_C_O
AM_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h math.h])

#
# For Python module
#
# if test -z $PYTHON; then
# 	PYTHON="python"
# fi
# PYTHON_INCDIR=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_inc()" 2>/dev/null`
# NUMPY_INCDIR=`$PYTHON -c "import numpy; print numpy.get_include()" 2>/dev/null`
# CPPFLAGS="$CPPFLAGS -I$PYTHON_INCDIR -I$NUMPY_INCDIR"
# AC_CHECK_HEADERS([Python.h])
# AC_CHECK_HEADERS([numpy/arrayobject.h])
# AC_CHECK_HEADERS([spglib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile
                 src/Makefile])
# AC_CONFIG_FILES([Makefile
#                  python/Makefile
#                  src/Makefile])
AC_OUTPUT
