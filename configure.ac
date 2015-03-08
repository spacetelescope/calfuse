#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.6])
AC_INIT([calfuse], [3.2.3], [archive@stsci.edu])
AM_INIT_AUTOMAKE([foreign])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([src/include/config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_PROG_LN_S
AC_CHECK_PROG(HAVE_JAVAC, javac, [`which javac`], [javac not found])

# Checks for libraries.

# Checks for header files.
AC_FUNC_ALLOCA
PKG_CHECK_MODULES([cfitsio], [cfitsio >= 2.7])
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h strings.h unistd.h fitsio.h])
# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memchr memset modf pow sqrt strcasecmp strerror strncasecmp strrchr strtol])

AC_OUTPUT([Makefile \
		   src/Makefile \
		   src/libcf/Makefile \
		   src/slalib/Makefile \
		   src/fes/Makefile \
		   src/fuv/Makefile \
		   src/analysis/Makefile \
		   src/cal/Makefile \
		   src/cal/jitter/Makefile \
		   src/cal/parm/Makefile \
		   src/cal/saa/Makefile \
		   src/cal/screen/Makefile \
		   src/cal/get_tle/Makefile \
		   src/cal/geom/Makefile \
		   scripts/Makefile \
		   doc/Makefile])

