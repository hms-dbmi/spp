AC_INIT([SPP], 1.7)

AC_CHECK_LIB(bz2, BZ2_bzDecompressInit)
AC_SUBST(HAVE_LIBBZ2)
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/tokenizer.hpp], [],
    [AC_MSG_ERROR(Unable to find required Boost libraries.)])
AC_LANG_POP([C++])
AC_CONFIG_FILES([src/Makevars])
cp confdefs.h src/config.h
AC_OUTPUT
