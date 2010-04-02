dnl -*- sh -*-
AC_INIT(src/nonthermal-module.c)
AC_PREFIX_DEFAULT(/usr/local)

AC_CONFIG_AUX_DIR(autoconf)
AC_CANONICAL_HOST
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Note that Fortran flags may depend on the C compiler (if gcc-3).
AC_PROG_FC

JD_INIT
JD_ANSI_CC
JD_ELF_COMPILER
JD_IEEE_CFLAGS

AC_FC_LIBRARY_LDFLAGS
AC_FC_DUMMY_MAIN(,
  AC_DEFINE(FC_DUMMY_MAIN,MAIN__)
)
JH_FC_DEFS

AC_PATH_XTRA

JD_WITH_LIBRARY(slang)

dnl# Add libraries here
JD_WITH_LIBRARY(gsl,gsl/gsl_const_cgsm.h)
JD_WITH_LIBRARY(isis)

dnl# This macro inits the module installation dir
JD_SLANG_MODULE_INSTALL_DIR

JD_GET_MODULES(src/modules.lis)
JH_GET_FC_MODULES(src/f77_modules.lis)

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

AC_ARG_WITH(cparamlib,
  [  --with-cparamlib[=DIR]  Use Kamae et al neutral pion-decay cross-sections ],
  [  jch_use_cparamlib=$withval ], [jch_use_cparamlib=no])
CPARAMLIB_INC=#
CPARAMLIB_LIB=#
if test "x$jch_use_cparamlib" != "xno" ; then
   CPARAMLIB_INC="-I${jch_use_cparamlib}/include"
   CPARAMLIB_LIB="-L${jch_use_cparamlib}/lib -lcparamlib"
   AC_DEFINE(USE_CPARAMLIB)
fi
AC_SUBST(CPARAMLIB_INC)
AC_SUBST(CPARAMLIB_LIB)

JD_SET_RPATH($libdir)
ELF_CFLAGS="$ELF_CFLAGS $IEEE_CFLAGS"
CFLAGS="$CFLAGS $IEEE_CFLAGS"

AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_OUTPUT(Makefile:autoconf/Makefile.in src/Makefile)

JD_CREATE_MODULE_ORULES(src)
JD_APPEND_ELFRULES(src/Makefile)

JH_CREATE_FC_MODULE_ORULES(src)
dnl JH_APPEND_FC_RULES(src/Makefile)
JH_APPEND_FC_ELFRULES(src/Makefile)

echo ""
echo "You are compiling with the following compiler configuration:"
echo "        CC =" "$CC"
echo " CC_SHARED =" "$CC_SHARED"
echo "    CFLAGS =" "$CFLAGS"
echo "   LDFLAGS =" "$LDFLAGS" "$DYNAMIC_LINK_FLAGS"
echo "      FC =" "$FC"
echo ""
echo "The modules will be installed in $MODULE_INSTALL_DIR."
echo "Any associated .sl files will be install in $SL_FILES_INSTALL_DIR"
echo ""
echo "If any of these quantities are incorrect, edit src/Makefile accordingly."
echo ""
