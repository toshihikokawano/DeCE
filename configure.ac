#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([dece], [1.2.2-Pyrite], [kawano@lanl.gov])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt strchr getopt])

AC_CONFIG_FILES([Makefile
                 source/Makefile])

AC_OUTPUT

