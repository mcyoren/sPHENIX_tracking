AC_INIT([PHGarfieldDriftVelocityQA], [1.00])
AC_CONFIG_SRCDIR([PHGarfieldDriftVelocityQA.cc])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PROG_CXX
LT_INIT([disable-static])
AS_IF([test "x$GXX" = xyes], [CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"])
AC_PATH_PROG([ROOT_CONFIG], [root-config])
AS_IF([test -z "$ROOT_CONFIG"], [AC_MSG_ERROR([root-config not found; source your sPHENIX environment first])])
ROOT_CFLAGS=`$ROOT_CONFIG --cflags`
ROOT_LIBS=`$ROOT_CONFIG --libs`
AC_SUBST([ROOT_CFLAGS])
AC_SUBST([ROOT_LIBS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
