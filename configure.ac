#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([calfuse], [3.2.3], [archive@stsci.edu])
AM_INIT_AUTOMAKE([foreign])
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([include/calfitsio.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_PROG_LN_S

# Checks for libraries.

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h strings.h unistd.h])
PKG_CHECK_MODULES([cfitsio], [cfitsio >= 3])
# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memchr memset modf pow sqrt strcasecmp strerror strncasecmp strrchr strtol])

AC_OUTPUT([Makefile src/Makefile src/libcf/Makefile src/slalib/Makefile \
		   src/fes/Makefile src/fuv/Makefile src/analysis/Makefile \
		   src/cal/Makefile src/cal/jitter/Makefile])
