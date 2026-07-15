AC_INIT([TpcDstV0Finder], [1.0])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([TpcDstV0Finder.cc])
AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
LT_INIT
PKG_PROG_PKG_CONFIG
AC_CHECK_PROG(ROOT_CONFIG, root-config, root-config, no)
AS_IF([test "x$ROOT_CONFIG" = "xno"], [AC_MSG_ERROR([root-config not found])])
ROOTCFLAGS=`$ROOT_CONFIG --cflags`
ROOTLIBS=`$ROOT_CONFIG --libs`
AC_SUBST(ROOTCFLAGS)
AC_SUBST(ROOTLIBS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
