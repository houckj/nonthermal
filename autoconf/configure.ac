dnl -*- sh -*-
AC_INIT(src/nonthermal-module.c)
AC_PREFIX_DEFAULT(/usr/local)

AC_CONFIG_AUX_DIR(autoconf)
AC_CANONICAL_HOST
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Note that Fortran flags may depend on the C compiler (if gcc-3).
AC_PROG_F77

JD_INIT
JD_ANSI_CC
JD_ELF_COMPILER
JD_IEEE_CFLAGS

AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN(,
  AC_DEFINE(F77_DUMMY_MAIN,MAIN__)
)
JH_F77_DEFS

AC_PATH_XTRA

JD_WITH_LIBRARY(slang)

dnl# Add libraries here
JD_WITH_LIBRARY(gsl)
JD_WITH_LIBRARY(isis)

dnl# This macro inits the module installation dir
JD_SLANG_MODULE_INSTALL_DIR

JD_GET_MODULES(src/modules.lis)
JH_GET_F77_MODULES(src/f77_modules.lis)

dnl Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)

JD_SET_RPATH($libdir)
ELF_CFLAGS="$ELF_CFLAGS $IEEE_CFLAGS"
CFLAGS="$CFLAGS $IEEE_CFLAGS"

AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_OUTPUT(Makefile:autoconf/Makefile.in src/Makefile)

JD_CREATE_MODULE_ORULES(src)
JD_APPEND_RULES(src/Makefile)

JH_CREATE_F77_MODULE_ORULES(src)
JH_APPEND_F77_RULES(src/Makefile)
dnl JH_APPEND_F77_ELFRULES(src/Makefile)

echo ""
echo "You are compiling with the following compiler configuration:"
echo "        CC =" "$CC"
echo " CC_SHARED =" "$CC_SHARED"
echo "    CFLAGS =" "$CFLAGS"
echo "   LDFLAGS =" "$LDFLAGS" "$DYNAMIC_LINK_FLAGS"
echo "      F77 =" "$F77"
echo ""
echo "The modules will be installed in $MODULE_INSTALL_DIR."
echo "Any associated .sl files will be install in $SL_FILES_INSTALL_DIR"
echo ""
echo "If any of these quantities are incorrect, edit src/Makefile accordingly."
echo ""
